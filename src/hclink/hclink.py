import csv
import gzip
import json
import lzma
import math
import os
import sys
from datetime import datetime
from pathlib import Path

import typer
from bitarray import bitarray
from bitarray.util import sc_encode, serialize
from typer import Argument, Option
from typing_extensions import Annotated

from hclink.build import Database, convert_to_profile, download_hiercc_profiles, download_profiles, get_species_scheme, \
    read_raw_hiercc_profiles, st_info
from hclink.search import calculate_hiercc_distance, imap_search, infer_hiercc_code, \
    read_gap_profiles, read_profiles, read_st_info

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
    lengths_db: Path = db_path / "profile_lengths.txt"
    gap_db: Path = db_path / "gap_profiles.xz"
    gap_lengths_db: Path = db_path / "gap_lengths.txt"
    st_db: Path = db_path / "ST.txt"

    print(f"Reading metadata from {scheme_metadata} ({datetime.now()})", file=sys.stderr)
    with open(scheme_metadata, "r") as scheme_metadata_fh:
        metadata = json.load(scheme_metadata_fh)

    print(f"Converting profile ({datetime.now()})", file=sys.stderr)
    query_profile: tuple[bitarray, bitarray] = convert_to_profile(query_code,
                                                                  metadata["array_size"],
                                                                  metadata["family_sizes"])
    print(f"Starting search ({datetime.now()})", file=sys.stderr)
    best_hit = imap_search(read_gap_profiles(gap_db, gap_lengths_db, len(metadata["family_sizes"])),
                           read_st_info(st_db),
                           read_profiles(lengths_db, profile_db),
                           metadata["max_gaps"],
                           query_profile,
                           num_threads,
                           batch_size)
    print(f"Finished search ({datetime.now()})", file=sys.stderr)
    # hiercc_distance: float = calculate_hiercc_distance(
    #     best_hit["distance"],
    #     best_hit["gaps_a"],
    #     best_hit["gaps_b"],
    #     best_hit["gaps_both"],
    #     len(metadata["family_sizes"]))

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


@app.command()
def build(
        version: str,
        api_key: str,
        species: Annotated[Database, typer.Option("-s", "--species")] = "ecoli",
        db_dir: Annotated[
            Path,
            Option("-d", "--downloads", help="Download directory", file_okay=False, dir_okay=True, writable=True,
                   readable=True)] = "db",
        clean: bool = False):
    scheme_info = get_species_scheme(species.value)
    # Need to fetch all profiles from enterobase
    # Write to sqlite database
    if not db_dir.exists():
        print(f"Creating directory '{db_dir}'", file=sys.stderr)
        db_dir.mkdir()
    hiercc_download: Path = download_hiercc_profiles(scheme_info["scheme"], api_key, db_dir)
    print(f"Downloaded HierCC profiles from '{scheme_info["scheme"]}'", file=sys.stderr)
    profiles_csv: Path = download_profiles(scheme_info["profiles"], db_dir)
    print(f"Downloaded profiles from '{profiles_csv}'", file=sys.stderr)
    write_db(version, profiles_csv, hiercc_download)

    if clean:
        profiles_csv.unlink()
        hiercc_download.unlink()


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
    path = Path(db_dir)
    if Path(metadata_json).exists():
        with open(metadata_json, 'r') as metadata_fh:
            metadata = json.load(metadata_fh)
            family_sizes = metadata["family_sizes"]
            array_size = metadata["array_size"]
            gap_slice = metadata["gap_slice"]
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
            gap_slice = len(bitarray(len(family_sizes)).tobytes())
            max_gaps = math.floor(len(family_sizes) * 0.1) + 1

    hiercc_profiles_data: tuple[dict[str, list[str]], str, list[int]] = read_raw_hiercc_profiles(
        hiercc_profiles_json)

    with open(metadata_json, 'w') as metadata_fh:
        print(json.dumps(
            {
                "array_size": array_size,
                "family_sizes": family_sizes,
                "gap_slice": gap_slice,
                "max_gaps": max_gaps,
                "thresholds": hiercc_profiles_data[2],
                "prepend": hiercc_profiles_data[1],
                "datestamp": str(datetime.now()),
                "version": version
            }), file=metadata_fh)

    with (gzip.open(profiles_csv, 'rt') as in_fh, lzma.open(path / "profiles.xz", 'wb') as profile_out,
          lzma.open(path / "gap_profiles.xz", 'wb') as gap_profile_out,
          open(path / "profile_lengths.txt", "w") as lengths_fh,
          open(path / "gap_lengths.txt", "w") as gap_lengths_fh,
          open(path / "ST.txt", "w") as st_db_out):
        reader = csv.reader(in_fh, delimiter="\t")
        next(reader)  # Skip header
        for row in reader:
            code = "_".join([value if value != "0" else "" for value in row[1:]])
            profile, gap_profile = convert_to_profile(code,
                                                      array_size,
                                                      family_sizes)
            if len(profile) != array_size:
                raise Exception(f"Profile bitarray length {len(profile)}!= {array_size}")
            if len(gap_profile.tobytes()) != gap_slice:
                raise Exception(f"Gap profile bytes length {len(gap_profile.tobytes())}!= {gap_slice}")
            encoded_profile = sc_encode(profile)
            profile_out.write(encoded_profile)
            lengths_fh.write(f"{len(encoded_profile)}\n")
            encoded_gap_profile = serialize(gap_profile)
            gap_profile_out.write(encoded_gap_profile)
            gap_lengths_fh.write(f"{len(encoded_gap_profile)}\n")
            st = row[0]
            st_db_out.write(f"{st_info(st, hiercc_profiles_data[0].get(st, []))}\n")


if __name__ == "__main__":
    app()
