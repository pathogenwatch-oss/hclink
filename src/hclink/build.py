import csv
import gzip
import json
import re
import sys
from datetime import datetime
from hashlib import sha1
from pathlib import Path
from typing import Any, Callable, Iterator

import numpy as np
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

from hclink.store import finalise_db, initialise_db


def st_info(st: str, hiercc_profile: list[str]) -> str:
    if int(st) < 0 or len(hiercc_profile) == 0:
        raise ValueError(f"Invalid ST code: {st} or profile {"_".join(hiercc_profile)}")
    else:
        return f'{st},{",".join(hiercc_profile)}'


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
    first_profile: dict[str, Any] = profiles[0]
    first_hiercc_key: str = next(iter(first_profile["info"]["hierCC"].keys()))
    prepend: str = re.sub(r'[0-9]+$', '', first_hiercc_key)
    thresholds: list[int] = sorted([int(key.replace(prepend, '')) for key in first_profile["info"]["hierCC"].keys()])
    for profile in profiles:
        if "info" not in profile or "hierCC" not in profile["info"]:
            continue
        st: str = profile["ST_id"]
        if int(st) < 1:
            continue
        processed[st]: list[str] = [item[1] for item in (
            sorted(
                [(int(item[0].replace(prepend, "")), item[1]) for item in profile["info"]["hierCC"].items()],
                key=lambda x: x[0]))]
    return processed, prepend, thresholds


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
    profiles_url: str = f"{scheme_url}/profiles.list.gz"
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
    out_dir: Path = db_dir / "alleles"
    out_dir.mkdir(parents=True, exist_ok=True)
    for gene in genes:
        gene_url: str = f"{url}/{gene}.fasta.gz"
        gene_fasta: Path = out_dir / f"{gene}.fasta.gz"
        if not gene_fasta.exists():
            download_resource(gene_url, gene_fasta)
    return out_dir


def hash_alleles(filename: Path, gene: str, idx: int, hash_size: int = 20) -> Iterator[tuple[str, int, int]]:
    with gzip.open(filename, "rt", encoding='utf-8') as fasta_fh:
        first_line = fasta_fh.readline()
        code = int(first_line.strip().replace(f">{gene}_", ""))
        for line in fasta_fh.readlines():
            if line.startswith(">"):
                code = int(line.strip().replace(f">{gene}_", ""))
            else:
                yield sha1(line.strip().lower().encode()).hexdigest()[0:hash_size], idx, code


def create_allele_db(genes, alleles_dir: Path, dbfile, hash_size: int = 20):
    db = initialise_db(dbfile)
    cursor = db.cursor()
    for idx, gene in enumerate(genes):
        filename = alleles_dir / f"{gene}.fasta.gz"
        gene = filename.stem.replace(".fasta", "")
        cursor.executemany("INSERT INTO alleles(checksum, position, code) VALUES(?,?,?)",
                           hash_alleles(filename, gene, idx, hash_size))
    db.commit()
    cursor.close()
    finalise_db(db)


def extract_genes(profiles_csv):
    with gzip.open(profiles_csv, "rt") as profiles_fh:
        reader = csv.reader(profiles_fh, delimiter='\t')
        header = next(reader)
        return header[1:]


def convert_to_vector(code: str,
                      lookup: Callable[[str, int], int],
                      hash_size: int = 20
                      ):
    code_arr = code.split("_")
    profile_array: list[float] = [0.0] * len(code_arr)
    for i in range(0, len(code_arr)):
        if code_arr[i].isnumeric():
            profile_array[i] = float(code_arr[i])
        elif code_arr[i] == "":
            profile_array[i] = 0.0
        else:
            st = lookup(code_arr[i][0:hash_size], i)
            if st is not None:
                profile_array[i] = float(st)
            else:
                profile_array[i] = float(10000000)
    return np.array(profile_array, dtype=np.float32)
