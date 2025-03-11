import csv
import gzip
import json
import lzma
import math
import os
import shutil
import sqlite3
import sys
from datetime import datetime
from enum import Enum
from functools import partial
from pathlib import Path
from typing import Any, Callable

import typer
from bitarray import bitarray
from bitarray.util import sc_encode, serialize
from typer import Argument, Option
from typing_extensions import Annotated

from hclink.search import imap_search, infer_hiercc_code
from hclink.store import connect_db, lookup_st, read_bitmaps, read_st_info, write_bitmap_to_filehandle
from hclink.build import convert_to_profile, create_allele_db, download_alleles, download_hiercc_profiles, \
    download_profiles, get_species_scheme, read_raw_hiercc_profiles, st_info

app = typer.Typer()


@app.command()
def assign(
        query: Annotated[
            str,
            Argument(
                help="file path of cgMLST JSON or '-' for JSON on STDIN"
            )],
        db: Annotated[
            Path,
            Option(
                "-db",
                "--reference-db",
                help="Reference profiles file. CSV of 'sample,ST,HierCC...,encoded profile'",
                file_okay=False, dir_okay=True, writable=False,
                readable=True
            )] = "db",
        num_threads: Annotated[
            int,
            Option(
                "-T",
                "--num-threads",
                help="Number of threads to use for multiprocessing"
            )] = os.cpu_count(),
        batch_size: Annotated[
            int,
            Option(
                "-B",
                "--batch-size",
                help="Number of profiles to process per batch on each thread"
            )
        ] = 1000,
) -> None:
    print(f"Starting ({datetime.now()})", file=sys.stderr)
    if query == '-':
        query_code: str = json.loads(sys.stdin.read())["code"]
    else:
        with open(query) as f:
            query_code: str = json.load(f)["code"]

    db_path = Path(db)
    scheme_metadata: Path = db_path / "metadata.json"
    profile_db: Path = db_path / "profiles.xz"
    gap_db: Path = db_path / "gap_profiles.xz"
    allele_db_path: Path = db_path / "alleles.db"
    st_db: Path = db_path / "ST.txt.xz"

    print(f"Reading metadata from {scheme_metadata} ({datetime.now()})", file=sys.stderr)
    with open(scheme_metadata, "r") as scheme_metadata_fh:
        metadata = json.load(scheme_metadata_fh)

    print(f"Converting profile ({datetime.now()})", file=sys.stderr)
    db = connect_db(str(allele_db_path))
    cursor: sqlite3.Cursor = db.cursor()
    lookup: Callable[[str, int], int] = partial(lookup_st, cursor)
    query_profile: tuple[bitarray, bitarray] = convert_to_profile(query_code,
                                                                  metadata["array_size"],
                                                                  metadata["family_sizes"],
                                                                  lookup)
    print(f"Starting search ({datetime.now()})", file=sys.stderr)
    best_hit: dict[str, Any] = imap_search(read_bitmaps(gap_db),
                           read_st_info(st_db),
                           read_bitmaps(profile_db),
                           metadata["max_gaps"],
                           query_profile,
                           num_threads,
                           batch_size)
    print(f"Finished search ({datetime.now()})", file=sys.stderr)


    print(f"Finished task ({datetime.now()})", file=sys.stderr)
    hiercc_code: list[tuple[str, str]] = infer_hiercc_code(best_hit["hiercc_distance"],
                                                           metadata["thresholds"],
                                                           best_hit["st"][1],
                                                           metadata["prepend"])
    print(json.dumps({
        "versions": {
            "hclink": metadata["version"],
            "library": metadata["datestamp"],
        },
        "closestST": best_hit["st"][0],
        "distance": best_hit["distance"],
        "hierccDistance": round(best_hit["hiercc_distance"], 2),
        "sharedGaps": best_hit["gaps_both"],
        "queryGaps": best_hit["gaps_a"],
        "referenceGaps": best_hit["gaps_b"],
        "hierCC": hiercc_code
    }), file=sys.stdout)


def extract_genes(profiles_csv):
    with gzip.open(profiles_csv, "rt") as profiles_fh:
        reader = csv.reader(profiles_fh, delimiter='\t')
        header = next(reader)
        return header[1:]


class Database(Enum):
    SENTERICA: str = "senterica"
    ECOLI: str = "ecoli"


@app.command()
def build(
        version: str,
        api_key: str,
        species: Annotated[Database, typer.Option("-s", "--species")],
        db_dir: Annotated[
            Path,
            Option("-d", "--downloads", help="Download directory", file_okay=False, dir_okay=True, writable=True,
                   readable=True)] = "db",
        clean: bool = False):
    scheme_info = get_species_scheme(species.value)
    if not db_dir.exists():
        print(f"Creating directory '{db_dir}'", file=sys.stderr)
        db_dir.mkdir()
    hiercc_download: Path = download_hiercc_profiles(scheme_info["scheme"], api_key, db_dir)
    print(f"Downloaded HierCC profiles from '{scheme_info["scheme"]}'", file=sys.stderr)
    profiles_csv: Path = download_profiles(scheme_info["downloads"], db_dir)
    print(f"Downloaded profiles from '{profiles_csv}'", file=sys.stderr)
    genes: list[str] = extract_genes(profiles_csv)
    alleles_dir: Path = download_alleles(scheme_info["downloads"], genes, db_dir)
    write_db(version, profiles_csv, hiercc_download, db_dir)

    if clean:
        profiles_csv.unlink()
        hiercc_download.unlink()
        shutil.rmtree(alleles_dir)


@app.command()
def write_db(
        version: str,
        profiles_csv: Annotated[
            Path, Option("-p", "--profiles-csv", help="Gzipped profiles CSV location (.gz)", exists=True, readable=True,
                         file_okay=True, dir_okay=False)] = "db/cgmlst_profiles.csv.gz",
        hiercc_profiles_json: Annotated[
            Path, Option("-c", "--hiercc-profiles-json", help="Compressed HierCC CSV location (.xz extension)",
                         exists=True, readable=True, file_okay=True, dir_okay=False)] = "db/hiercc_profiles.json.gz",
        db_dir: Annotated[Path, Option("-d", help="Data directory", readable=True, writable=True, file_okay=False,
                                       dir_okay=True)] = "db",
        metadata_json: str = "db/metadata.json"
):
    print(f"Writing database to '{db_dir}'", file=sys.stderr)
    path = Path(db_dir)
    if Path(metadata_json).exists():
        with open(metadata_json, 'r') as metadata_fh:
            metadata = json.load(metadata_fh)
            family_sizes = metadata["family_sizes"]
            array_size = metadata["array_size"]
            max_gaps = metadata["max_gaps"]
    else:
        with gzip.open(profiles_csv, 'rt') as in_fh:
            reader = csv.reader(in_fh, delimiter='\t')
            row = next(reader)
            family_sizes: list[int] = [0] * (len(row) - 1)
            for row in reader:
                for i in range(1, len(row)):
                    if int(row[i]) > family_sizes[i - 1]:
                        family_sizes[i - 1] = int(row[i])

            array_size = sum(family_sizes) + len(family_sizes)
            max_gaps = math.floor(len(family_sizes) * 0.1) + 1

    genes: list[str] = extract_genes(profiles_csv)

    create_allele_db(genes, db_dir / "alleles", db_dir / "alleles.db")

    print("Generating profiles", file=sys.stderr)
    hiercc_profiles_data: tuple[dict[str, list[str]], str, list[int]] = read_raw_hiercc_profiles(
        hiercc_profiles_json)

    with open(metadata_json, 'w') as metadata_fh:
        print(json.dumps(
            {
                "array_size": array_size,
                "family_sizes": family_sizes,
                "max_gaps": max_gaps,
                "thresholds": hiercc_profiles_data[2],
                "prepend": hiercc_profiles_data[1],
                "datestamp": str(datetime.now()),
                "version": version
            }), file=metadata_fh)

    with (gzip.open(profiles_csv, 'rt') as in_fh, lzma.open(path / "profiles.xz", 'wb') as profile_out,
          lzma.open(path / "gap_profiles.xz", 'wb') as gap_profile_out,
          lzma.open(path / "ST.txt.xz", "wt") as st_db_out):
        reader = csv.reader(in_fh, delimiter="\t")
        next(reader)  # Skip header
        store_profile = partial(write_bitmap_to_filehandle, profile_out, sc_encode)
        store_gaps = partial(write_bitmap_to_filehandle, gap_profile_out, serialize)
        for row in reader:
            code = "_".join([value if value != "0" else "" for value in row[1:]])
            profile, gap_profile = convert_to_profile(code, array_size, family_sizes, lambda x: None)
            if len(profile) != array_size:
                raise Exception(f"Profile bitarray length {len(profile)}!= {array_size}")
            if len(gap_profile) != len(family_sizes):
                raise Exception(f"Gap profile bytes length {len(gap_profile)}!= {len(family_sizes)}")
            store_profile(profile)
            store_gaps(gap_profile)
            st = row[0]
            st_db_out.write(f"{st_info(st, hiercc_profiles_data[0].get(st, []))}\n")


if __name__ == "__main__":
    app()
