import gzip
import json
import re
import sys
from datetime import datetime
from hashlib import sha1
from pathlib import Path
from typing import Any, Iterator

import requests
from bitarray import bitarray
from tenacity import retry, stop_after_attempt, wait_exponential

from hclink.store import finalise_db, initialise_db


def st_info(st: str, hiercc_profile: list[str]) -> str:
    if int(st) < 0 or len(hiercc_profile) == 0:
        return f'{st},{",".join([""] * 14)}'
    else:
        return f'{st},{",".join(hiercc_profile)}'


def convert_to_profile(code, array_size: int, family_sizes: list[int], lookup) -> tuple[bitarray, bitarray]:
    """
    Converts a CGMLST code to a pair of bitarrays: one for the profile and one for the gap positions.
    The bitarray is the length of the total number of alleles across all loci +
    the size of the scheme. Effectively, the bitarray made up of one bitarray per family. The bit in the sub-arrays is
    set according the allele code, with the final bit set if it is a novel code. The sub-arrays are joined into a single
    bitarray.

    :param code:
    :param array_size:
    :param family_sizes: list of highest allele ST for each position in the profile
    :lookup: a function that takes the checksum and locus index to return the ST code.
    :return: (profile_array, gap_array)
    """
    code_arr = code.split("_")
    profile_array = bitarray(array_size)
    gap_array = bitarray(len(family_sizes))
    offset = 0
    for i in range(0, len(code_arr)):
        if code_arr[i].isnumeric():
            index = int(code_arr[i])
            if index > 0:
                profile_array[offset + index - 1] = 1
            else:
                gap_array[i] = 1
        elif code_arr[i] != "":
            # Convert the code
            st = lookup(code_arr[i],i)
            if st is not None:
                profile_array[offset + st - 1] = 1
            else:
                profile_array[offset + family_sizes[i]] = 1
        else:
            gap_array[i] = 1
        offset += family_sizes[i] + 1
    if len(profile_array) != array_size:
        raise Exception(f"Profile bitarray length {len(profile_array)}!= {array_size}")
    return profile_array, gap_array


def get_species_scheme(species: str, filepath: Path = Path("schemes.json")) -> dict[str, str]:
    with open(filepath) as f:
        schemes = json.load(f)
    if species not in schemes["schemes"]:
        raise ValueError(f"Species '{species}' not found in schemes.json")
    return {"scheme": schemes["schemes"].get(species, ""), "downloads": schemes["downloads"].get(species, "")}


def read_raw_hiercc_profiles(hiercc_profiles_json: Path) -> tuple[dict[str, list[str]], str, list[int]]:
    with gzip.open(hiercc_profiles_json, "rt") as hiercc_profiles_fh:
        profiles = json.loads(hiercc_profiles_fh.read())
    processed: dict[str, list[str]] = {}
    first_profile = profiles[0]
    first_hiercc_key = next(iter(first_profile["info"]["hierCC"].keys()))
    prepend = re.sub(r'[0-9]+$', '', first_hiercc_key)
    thresholds: list[int] = sorted([int(key.replace(prepend, '')) for key in first_profile["info"]["hierCC"].keys()])
    for profile in profiles:
        if "info" not in profile or "hierCC" not in profile["info"]:
            continue
        st = profile["ST_id"]
        if int(st) < 1:
            continue
        processed[st] = [item[1] for item in (
            sorted(
                [(int(item[0].replace("d", "")), item[1]) for item in profile["info"]["hierCC"].items()],
                key=lambda x: x[0]))]
    return processed, prepend, thresholds


# @retry(
#     wait=wait_exponential(multiplier=1, min=4, max=240),
#     stop=stop_after_attempt(10),
# )
# def download_resource(url: str, out_file: Path) -> Path:
#     try:
#         with requests.get(url, stream=True, timeout=30) as response:
#             response.raise_for_status()
#             with gzip.open(out_file, 'wb') as f_out:
#                 for chunk in response.iter_content(chunk_size=8192):
#                     if chunk:  # filter out keep-alive new chunks
#                         f_out.write(chunk)
#     except requests.exceptions.RequestException as e:
#         raise Exception(f"Failed to download '{url}': {str(e)}")
#     return out_file


@retry(
    wait=wait_exponential(multiplier=1, min=4, max=240),
    stop=stop_after_attempt(2),
)
def download_resource(url: str, output_path: Path) -> Path:
    """
    Download a file from a URL to a specified path.

    Args:
    url (str): The URL of the file to download.
    output_path (Path): The path where the file should be saved.
    is_gzipped (bool): Whether the file is gzipped. Defaults to False.

    Returns:
    Path: The path where the file was saved.

    Raises:
    Exception: If there's an error during the download process.
    """
    try:
        with requests.get(url, stream=True, timeout=30) as response:
            response.raise_for_status()
            with open(output_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

        return output_path
    except requests.exceptions.RequestException as e:
        raise Exception(f"Failed to download file from '{url}': {str(e)}")


def download_profiles(scheme_url: str, data_dir: Path) -> Path:
    profiles_csv: Path = data_dir / "cgmlst_profiles.csv.gz"
    profiles_url = f"{scheme_url}/profiles.list.gz"
    return download_resource(profiles_url, profiles_csv)


@retry(wait=wait_exponential(multiplier=1, min=10, max=7200))
def fetch_hiercc_batch(url: str, api_key: str, offset: int, limit: int) -> tuple[int, dict[str, Any]]:
    r = requests.get(
        f"{url}&limit={limit}&offset={offset}",
        headers={"Authorization": f"Basic {api_key}"}
    )
    if r.status_code != 200:
        print(r, file=sys.stderr)
        r.raise_for_status()
    return r.status_code, r.json()


def download_hiercc_profiles(
        url: str,
        api_key: str,
        data_dir: Path = "db",
        limit: int = 10000,
        safety_valve: int = 1000000
) -> Path:
    out_file = data_dir / "hiercc_profiles.json.gz"
    print(f"Downloading HierCC profiles from '{url}' to {out_file}", file=sys.stderr)
    with gzip.open(out_file, 'wt') as out_fh:
        sts = []
        offset: int = 0
        while offset < safety_valve:
            status, batch = fetch_hiercc_batch(url, api_key, offset, limit)

            if batch is None or not batch or not batch["STs"]:
                break
            sts.extend(batch["STs"])
            offset += limit
            print(f"{datetime.now()},{offset},{len(sts)}", file=sys.stderr)
        out_fh.write(json.dumps(sts))

    return out_file


def download_alleles(url: str, genes: list[str], db_dir: Path = "db") -> Path:
    out_dir = db_dir / "alleles"
    out_dir.mkdir(parents=True, exist_ok=True)
    for gene in genes:
        gene_url = f"{url}/{gene}.fasta.gz"
        gene_fasta = out_dir / f"{gene}.fasta.gz"
        if not gene_fasta.exists():
            download_resource(gene_url, gene_fasta)
    return out_dir


def hash_alleles(filename: Path, gene: str, idx: int) -> Iterator[tuple[str, int, int]]:
        with gzip.open(filename, "rt", encoding='utf-8') as fasta_fh:
            first_line = fasta_fh.readline()
            code = int(first_line.strip().replace(f">{gene}_", ""))
            for line in fasta_fh.readlines():
                if line.startswith(">"):
                    code = int(line.strip().replace(f">{gene}_", ""))
                else:
                    yield sha1(line.strip().lower().encode()).hexdigest(), idx, code


def create_allele_db(genes, alleles_dir: Path, dbfile):
    db = initialise_db(dbfile)
    cursor = db.cursor()
    for idx, gene in enumerate(genes):
        filename = alleles_dir / f"{gene}.fasta.gz"
        gene = filename.stem.replace(".fasta", "")
        cursor.executemany("INSERT INTO alleles(checksum, position, code) VALUES(?,?,?)", hash_alleles(filename, gene, idx))
    db.commit()
    cursor.close()
    finalise_db(db)
