import csv
import gzip
import json
import lzma
import math
import shutil
import sqlite3
import sys
import traceback
from datetime import datetime
from enum import Enum
from functools import partial
from pathlib import Path
from typing import Callable

import numpy as np
import typer
from typer import Argument, Option
from typing_extensions import Annotated
from usearch.index import Match, Matches

from hclink.build import convert_to_vector, create_allele_db, download_alleles, download_hiercc_profiles, \
    download_profiles, extract_genes, get_species_scheme, read_raw_hiercc_profiles
from hclink.search import count_shared_zeros, count_zeros, infer_hiercc_code, select_best
from hclink.store import batch_read_profiles, connect_db, create_index, lookup_st, read_st_info, restore_index

app = typer.Typer()


class Database(Enum):
    SENTERICA: str = "senterica"
    ECOLI: str = "ecoli"


def format_match(query_vector: np.ndarray, match: Match, get_matched_vector, get_st_info):
    if match is None:
        gaps: int = count_zeros(query_vector)
        return {
            "st": ("", []),
            "distance": query_vector.size,
            "hiercc_distance": query_vector.size,
            "gaps_a": gaps,
            "gaps_b": -1,
            "gaps_both": gaps,
        }
    else:
        matched_vector = get_matched_vector(match.key)
        gaps_a = count_zeros(query_vector)
        gaps_b = count_zeros(matched_vector)
        gaps_both = count_shared_zeros(query_vector, matched_vector)
        return {
            "st": get_st_info(str(match.key)),
            "distance": match.distance.item(),
            "hiercc_distance": match.distance.item(),
            "gaps_a": gaps_a,
            "gaps_b": gaps_b,
            "gaps_both": gaps_both,
        }


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
                "-d",
                "--reference-db",
                help="Reference profiles file. CSV of 'sample,ST,HierCC...,encoded profile'",
                file_okay=False, dir_okay=True, writable=False,
                readable=True
            )] = "db",
        num_matches: Annotated[int, Option("-m", "--max-matches", help="Maximum number of matches to return")] = 10,
) -> None:
    print(f"Starting ({datetime.now()})", file=sys.stderr)
    if query == '-':
        query_code: str = json.loads(sys.stdin.read())["code"]
    else:
        with open(query) as f:
            query_code: str = json.load(f)["code"]

    db_path = Path(db)
    scheme_metadata: Path = db_path / "metadata.json"
    allele_db_path: Path = db_path / "alleles.db"
    st_db: Path = db_path / "ST.txt.xz"

    # print(f"Reading metadata from {scheme_metadata} ({datetime.now()})", file=sys.stderr)
    with open(scheme_metadata, "r") as scheme_metadata_fh:
        metadata = json.load(scheme_metadata_fh)

    # print(f"Converting profile ({datetime.now()})", file=sys.stderr)
    db = connect_db(str(allele_db_path))
    cursor: sqlite3.Cursor = db.cursor()
    lookup: Callable[[str, int], int] = partial(lookup_st, cursor)
    query_vector = convert_to_vector(query_code, lookup)
    query_gaps = count_zeros(query_vector)
    if query_gaps >= metadata["max_gaps"]:
        best_hit = {
            "st": ("", []),
            "distance": query_vector.size,
            "hiercc_distance": query_vector.size,
            "gaps_a": query_gaps,
            "gaps_b": -1,
            "gaps_both": query_gaps,
        }
    else:
        print(f"Starting search ({datetime.now()})", file=sys.stderr)
        index = restore_index(db_path / "index.usearch", max_gaps=metadata["max_gaps"])
        matches: Matches = index.search(query_vector, num_matches)
        print(f"Finished search ({datetime.now()})", file=sys.stderr)
        st_info = partial(read_st_info, st_db)
        best_hit = format_match(query_vector, select_best(matches), index.get, st_info)
        # print(f"Finished formatting ({datetime.now()})", file=sys.stderr)

    hiercc_code: list[tuple[str, str]] = infer_hiercc_code(best_hit["hiercc_distance"],
                                                           metadata["thresholds"],
                                                           best_hit["st"][1],
                                                           metadata["prepend"])
    print(f"Finished task ({datetime.now()})", file=sys.stderr)
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
        db_dir: Annotated[
            Path, Option("-d", "--database-dir", help="Data directory", readable=True, writable=True, file_okay=False,
                         dir_okay=True)] = "db",
        index_batch_size: Annotated[int, Option("-B", "--index-batch-size", help="Batch size for indexing")] = 10000,
):
    print(f"Writing database to '{db_dir}'", file=sys.stderr)
    metadata_json = db_dir / "metadata.json"
    if metadata_json.exists():
        with open(metadata_json, 'r') as metadata_fh:
            try:
                metadata = json.load(metadata_fh)
            except json.JSONDecodeError as jde:
                print(f"Error: Invalid JSON in metadata file. {jde}", file=sys.stderr)

            num_families = metadata["num_families"]
            max_gaps = metadata["max_gaps"]
    else:
        with gzip.open(profiles_csv, 'rt') as in_fh:
            reader = csv.reader(in_fh, delimiter='\t')
            row = next(reader)
            num_families = len(row) - 1  # subtract 1 for the ST column
            max_gaps = math.floor(num_families * 0.1) + 1

    genes: list[str] = extract_genes(profiles_csv)
    print("Creating allele database", file=sys.stderr)

    allele_db = db_dir / "alleles.db"
    if not allele_db.exists():
        create_allele_db(genes, db_dir / "alleles", allele_db)

    print("Generating profiles", file=sys.stderr)
    hiercc_profiles_data: tuple[dict[str, list[str]], str, list[int]] = read_raw_hiercc_profiles(hiercc_profiles_json)
    with lzma.open(db_dir / "ST.txt.xz", "wt") as st_db_out:
        for st, hiercc_profile in hiercc_profiles_data[0].items():
            print(f"{st},{','.join(hiercc_profile)}", file=st_db_out)

    with open(metadata_json, 'w') as metadata_fh:
        print(json.dumps(
            {
                "num_families": num_families,
                "max_gaps": max_gaps,
                "thresholds": hiercc_profiles_data[2],
                "prepend": hiercc_profiles_data[1],
                "datestamp": str(datetime.now()),
                "version": version
            }), file=metadata_fh)

    try:
        print(f"Starting build function with db={db_dir}, batch_size={index_batch_size}", file=sys.stderr)
        print(f"Creating index with ndim={num_families}", file=sys.stderr)
        index = create_index(num_families, max_gaps=max_gaps)

        count = 0
        for keys, profiles in batch_read_profiles(db_dir / "cgmlst_profiles.csv.gz", max_gaps, batch_size=index_batch_size):
            index.add(keys, profiles)
            count += index_batch_size
            # if count % index_batch_size == 0:
            print(f"Indexed {count:,} profiles", file=sys.stderr)

        print("Saving index...", file=sys.stderr)
        index.save(db_dir / "index.usearch")
        print(f"Successfully indexed {count:,} profiles and saved to index.usearch", file=sys.stderr)

    except FileNotFoundError as fnf:
        print(f"Error: File not found. {fnf}", file=sys.stderr)
        sys.exit(1)
    except gzip.BadGzipFile as bgf:
        print(f"Error: Invalid gzip file. {bgf}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    app()
