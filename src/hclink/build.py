import gzip
import json
import math
import re
import shutil
import socket
import ssl
import sys
import urllib.request
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

import requests
from bitarray import bitarray
from tenacity import retry, wait_exponential


def st_info(st: str, hiercc_profile: list[str]) -> str:
    if int(st) < 0 or len(hiercc_profile) == 0:
        return f'{st},{",".join([""] * 14)}'
    else:
        return f'{st},{",".join(hiercc_profile)}'


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


def get_species_scheme(species: str, filepath: Path = Path("schemes.json")) -> dict[str, str]:
    with open(filepath) as f:
        schemes = json.load(f)
    if species not in schemes["schemes"]:
        raise ValueError(f"Species '{species}' not found in schemes.json")
    return {"scheme": schemes["schemes"].get(species, ""), "profiles": schemes["profiles"].get(species, "")}


class Database(Enum):
    SENTERICA: str = "senterica"
    ECOLI: str = "ecoli"


def read_raw_hiercc_profiles(hiercc_profiles_json: Path) -> tuple[dict[str, list[str]], str, list[int], int]:
    with gzip.open(hiercc_profiles_json, "rt") as hiercc_profiles_fh:
        profiles = json.loads(hiercc_profiles_fh.read())
    processed: dict[str, list[str]] = {}
    first_profile = profiles[0]
    first_hiercc_key = next(iter(first_profile["info"]["hierCC"].keys()))
    prepend = re.sub(r'[0-9]+$', '', first_hiercc_key)
    thresholds: list[int] = sorted([int(key.replace(prepend, '')) for key in first_profile["info"]["hierCC"].keys()])
    max_gaps = math.floor(len(first_profile["info"]["hierCC"].keys()) / 10) + 1
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
    return processed, prepend, thresholds, max_gaps


def download_profiles(profiles_url: str, data_dir: Path) -> Path:
    profiles_csv: Path = data_dir / "cgmlst_profiles.csv.gz"
    req = urllib.request.Request(
        profiles_url,
        data=None,
        headers={
            'User-Agent': 'mlst-downloader (https://gist.github.com/bewt85/16f2b7b9c3b331f751ce40273240a2eb)'
        }
    )
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        r = urllib.request.urlopen(req, context=ctx, timeout=10)
    except KeyboardInterrupt:
        raise
    except socket.timeout:
        raise Exception(f"GET '{profiles_url}' timed out after {10} seconds")
    if r.getcode() != 200:
        raise Exception(f"GET '{profiles_url}' returned {r.getcode()}")
    with open(profiles_csv, 'wb') as out_fh:
        shutil.copyfileobj(r, out_fh)
    return profiles_csv


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
