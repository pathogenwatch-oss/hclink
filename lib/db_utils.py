import gzip
import json
import lzma
import shutil
import socket
import ssl
import sys
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any

import requests
from bitarray import bitarray
from bitarray.util import deserialize
from tenacity import retry, wait_exponential


def read_raw_hiercc_profiles(hiercc_profiles_json: Path) -> dict[str, list[str]]:
    with gzip.open(hiercc_profiles_json, "rt") as hiercc_profiles_fh:
        profiles = json.loads(hiercc_profiles_fh.read())
        processed: dict[str, list[str]] = {}
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
    return processed



def download_profiles(profiles_url: str, data_dir: Path) -> Path:
    profiles_csv: Path = data_dir / "cgmlst_profiles.csv.gz"
    # if profiles_csv.exists():
    #     return profiles_csv
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


def read_reference_profiles(lengths_db, profile_db) -> list[bytes]:
    reference_profiles: list[bytes] = []
    with open(lengths_db, "r") as lengths_fh, lzma.open(profile_db, "rb") as reference_profiles_fh:
        for index, length_str in enumerate(lengths_fh.readlines()):
            encoded_profile = reference_profiles_fh.read(int(length_str.strip()))
            reference_profiles.append(encoded_profile)
    return reference_profiles


def read_gap_profiles(gap_db, gap_lengths_db, num_families) -> list[bitarray]:
    gap_profiles: list[bitarray] = []
    with open(gap_lengths_db, "r") as gap_lengths_fh, lzma.open(gap_db, "rb") as gap_profiles_fh:
        for index, length_str in enumerate(gap_lengths_fh.readlines()):
            encoded_profile = gap_profiles_fh.read(int(length_str.strip()))
            gap_profiles.append(deserialize(encoded_profile))
            if len(gap_profiles[index]) != num_families:
                raise Exception(
                    f"Profile {index} bitarray length {len(gap_profiles[index])} != {num_families}")
    return gap_profiles
