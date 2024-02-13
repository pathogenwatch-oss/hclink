import csv
import gzip
import json
import lzma
import multiprocessing
import os
import sys
from datetime import datetime
from itertools import repeat
from pathlib import Path
from typing import Any, Iterable

import typer
from bitarray import bitarray
from bitarray.util import count_and, count_xor, sc_decode, sc_encode, serialize
from typer import Argument, Option
from typing_extensions import Annotated

from lib.db_utils import (download_hiercc_profiles, download_profiles,
                          read_gap_profiles, read_raw_hiercc_profiles, read_reference_profiles)

app = typer.Typer()


@app.command()
def build(
        version: str,
        api_key: str,
        db_dir: Annotated[
            Path,
            Option("-s", "--downloads", help="Download directory", file_okay=False, dir_okay=True, writable=True,
                   readable=True)] = "db",
        clean: bool = False):
    # Need to fetch all profiles from enterobase
    # Write to sqlite database
    if not db_dir.exists():
        db_dir.mkdir()
    profiles_csv: Path = download_profiles(db_dir)
    hiercc_download: Path = download_hiercc_profiles(api_key, db_dir)
    write_db(version, profiles_csv, hiercc_download)

    if clean:
        profiles_csv.unlink()
        hiercc_download.unlink()


def st_info(st: str, hiercc_profile: list[str]) -> str:
    if int(st) < 0 or len(hiercc_profile) == 0:
        return f'{st},{",".join([""] * 14)}'
    else:
        return f'{st},{",".join(hiercc_profile)}'


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

        with open(path / "metadata.json", 'w') as metadata_fh:
            print(json.dumps(
                {
                    "array_size": array_size,
                    "family_sizes": family_sizes,
                    "gap_slice": gap_slice,
                    "datestamp": str(datetime.now()),
                    "version": version
                }), file=metadata_fh)

    hiercc_profiles: dict[str, list[str]] = read_raw_hiercc_profiles(hiercc_profiles_json)

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
            st_db_out.write(f"{st_info(st, hiercc_profiles.get(st, []))}\n")


def convert_to_profile(code, array_size: int, family_sizes: list[int]) -> tuple[bitarray, bitarray]:
    """
    Converts a CGMLST code to a pair of bitarrays: one for the profile and one for the gap positions.
    The bitarray is the length of the total number of alleles across all loci +
    the size of the scheme. Effectively, the bitarray made up of one bitarray per family. The bit in the sub-arrays is
    set according the allele code, with the final bit set if it is a novel code. The sub-arrays are joined into a single
    bitarray.
    
    :param code: 
    :param array_size:
    :param family_sizes: list of highest allele ST for each position in the profile
    :return: (profile_array, gap_array)
    """
    code_arr = code.split("_")
    profile_array = bitarray(array_size)
    gap_array = bitarray(len(family_sizes))
    offset = 0
    for i in range(0, len(code_arr)):
        if code_arr[i].isnumeric():
            profile_array[offset + int(code_arr[i]) - 1] = 1
        elif code_arr[i] != "":
            profile_array[offset + family_sizes[i]] = 1
        else:
            gap_array[i] = 1
        offset += family_sizes[i] + 1
    if len(profile_array) != array_size:
        raise Exception(f"Profile bitarray length {len(profile_array)}!= {array_size}")
    return profile_array, gap_array


def compare_profiles(
        profile_a: tuple[bitarray, bitarray],
        profile_b: tuple[bitarray, bitarray]
) -> tuple[int, int, int, int]:
    """
    Compares two profiles and returns the distance between them.
    The distance is the number of bits that are different between the two profiles, excluding the gap positions.
    """
    shared_gaps: int = count_and(profile_a[1], profile_b[1])
    profileA_gaps: int = profile_a[1].count() - shared_gaps
    profileB_gaps: int = profile_b[1].count() - shared_gaps
    gap_adjust: int = profileA_gaps + profileB_gaps
    distance = count_xor(profile_a[0], profile_b[0])
    return int((distance - gap_adjust) / 2), profileA_gaps, profileB_gaps, shared_gaps


def comparison(
        query_profile: tuple[bitarray, bitarray],
        reference_profile,
        gap_profile,
        st: tuple[str, list[str]]
):
    profile = sc_decode(reference_profile)
    distance, gaps_a, gaps_b, gaps_both = compare_profiles(query_profile, (profile, gap_profile))
    return distance, gaps_a, gaps_b, gaps_both, st


def calculate_hiercc_distance(distance, query_gaps, reference_gaps, shared_gaps, profile_size) -> float:
    # s = shared genes - 0.03 * profile_size

    if distance == 0 and query_gaps == 0 and reference_gaps == 0:
        cc_distance = 0.0
    else:
        query_core = (profile_size - query_gaps - shared_gaps) - 0.03 * profile_size
        # if query core is > s use equation 1, else use equation 2
        common_core = profile_size - query_gaps - reference_gaps - shared_gaps
        if common_core >= query_core:
            cc_distance = (profile_size * distance) / common_core + 0.5
        else:
            cc_distance = ((profile_size * (distance + query_core - common_core)) / query_core) + 0.5
    return cc_distance


def infer_hiercc_code(hier_cc_distance: float, hier_cc_map: Iterable[tuple[int, str]]) -> Iterable[tuple[str, str]]:
    for threshold, code in hier_cc_map:
        formatted_threshold = f"d{threshold}"
        if hier_cc_distance <= threshold:
            yield formatted_threshold, code
        else:
            yield formatted_threshold, ""


@app.command()
def assign(
        query: Annotated[str, Argument(help="file path of cgMLST JSON or '-' for JSON on STDIN")],
        db: Annotated[
            Path, Option("-db", "--reference-db",
                         help="Reference profiles file. CSV of 'sample,ST,HierCC...,encoded profile'",
                         file_okay=False, dir_okay=True, writable=False,
                         readable=True)] = "db",
        num_cpu: Annotated[
            int, Option("-n", "--num-cpu",
                        help="Number of CPUs to use for multiprocessing")
        ] = os.cpu_count(),
        max_gaps: Annotated[
            int, Option("-g", "--max-gaps",
                        help="Pairs of profiles with equal or more combined positions are ignored")
        ] = 301):
    if query == '-':
        query_json: dict[str, Any] = json.loads(sys.stdin.read())
    else:
        with open(query) as f:
            query_json: dict[str, Any] = json.load(f)

    db_path = Path(db)
    scheme_metadata = db_path / "metadata.json"
    profile_db = db_path / "profiles.xz"
    lengths_db = db_path / "profile_lengths.txt"
    gap_db = db_path / "gap_profiles.xz"
    gap_lengths_db = db_path / "gap_lengths.txt"
    st_db = db_path / "ST.txt"

    with open(scheme_metadata, "r") as scheme_metadata_fh:
        metadata = json.load(scheme_metadata_fh)

    array_size = metadata["array_size"]
    query_profile: tuple[bitarray, bitarray] = convert_to_profile(query_json["code"], array_size,
                                                                  metadata["family_sizes"])

    gap_profiles: list[bitarray] = read_gap_profiles(
        gap_db,
        gap_lengths_db,
        len(metadata["family_sizes"])
    )

    reference_profiles = read_reference_profiles(lengths_db, profile_db)

    best_hits: list[dict[str, Any]] = []
    lowest_distance: int = sys.maxsize

    with (open(st_db, "r") as st_db_fh, multiprocessing.Pool(num_cpu) as pool):
        sts: list[tuple[str, list[str]]] = []
        for line in st_db_fh.readlines():
            info = line.strip().split(",")
            sts.append((info[0], info[1:]))
        for distance, gaps_a, gaps_b, gaps_both, st in pool.starmap(comparison,
                                                                    zip(repeat(query_profile), reference_profiles,
                                                                        gap_profiles, sts)):

            total_gaps = gaps_a + gaps_b + gaps_both
            if total_gaps >= max_gaps:
                continue
            # for (distance, st) in results:
            # for reference_profile, st, gap_profile in zip(reference_profiles, sts, gap_profiles):
            if distance < lowest_distance:
                best_hits = [
                    {"st": st, "distance": distance, "gaps_a": gaps_a,
                     "gaps_b": gaps_b, "gaps_both": gaps_both}]
                lowest_distance = distance
            elif distance == lowest_distance:
                for hit in best_hits:
                    if total_gaps < hit["gaps_both"] + hit["gaps_a"] + hit["gaps_b"]:
                        best_hits = [
                            {"st": st, "distance": distance, "gaps_a": gaps_a,
                             "gaps_b": gaps_b, "gaps_both": gaps_both}]
                        break
                    elif total_gaps == hit["gaps_both"] + hit["gaps_a"] + hit["gaps_b"]:
                        best_hits.append(
                            {"st": st, "distance": distance, "gaps_a": gaps_a,
                             "gaps_b": gaps_b, "gaps_both": gaps_both})
                        break

    # Select a single best match. For now keep it simple
    best_hit: dict[str, Any] = best_hits[0]

    # profile_index: int = sts.index(best_hit["st"])
    hier_cc_distance: float = calculate_hiercc_distance(
        best_hit["distance"],
        best_hit["gaps_a"],
        best_hit["gaps_b"],
        best_hit["gaps_both"],
        len(metadata["family_sizes"]))

    hiercc_code: list[tuple[str, str]] = list(infer_hiercc_code(
        hier_cc_distance,
        zip([0, 2, 5, 10, 20, 50, 100, 150, 200, 400, 900, 2000, 2600, 2850], best_hit["st"][1])))

    print(json.dumps({
        "versions": {
            "hclink": metadata["version"],
            "library": metadata["datestamp"],
        },
        "libraryVersion": metadata["datestamp"],
        "closestST": best_hit["st"][0],
        "distance": best_hit["distance"],
        "hierccDistance": round(hier_cc_distance, 2),
        "sharedGaps": best_hit["gaps_both"],
        "queryGaps": best_hit["gaps_a"],
        "referenceGaps": best_hit["gaps_b"],
        "hierCC": hiercc_code
    }), file=sys.stdout)


if __name__ == "__main__":
    app()
