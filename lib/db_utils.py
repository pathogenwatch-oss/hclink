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
from retry import retry


def read_raw_hiercc_profiles(hiercc_profiles_json: Path) -> dict[str, list[str]]:
    with lzma.open(hiercc_profiles_json, "r") as hiercc_profiles_fh:
        profiles = json.loads(hiercc_profiles_fh.read().decode('utf-8'))
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


profiles_url = "https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/profiles.list.gz"


def download_profiles(data_dir: Path) -> Path:
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


@retry(tries=3, delay=1, backoff=5)
def fetch_hiercc_batch(key: str, offset: int, limit: int) -> tuple[int, dict[str, Any]]:
    r = requests.get(
        f"https://enterobase.warwick.ac.uk/api/v2.0/senterica/cgMLST_v2/sts?limit={limit}&offset={offset}&scheme=cgMLST_v2",
        headers={"Authorization": f"Basic {key}"}
    )
    if r.status_code != 200:
        print(r, file=sys.stderr)
        r.raise_for_status()
    return r.status_code, r.json()


def download_hiercc_profiles(key: str, data_dir: Path = "db", limit: int = 10000, safety_valve: int = 250000) -> Path:
    out_file = data_dir / "hiercc_profiles.json.xz"
    # curl -X GET --header "Accept: application/json" --header "Authorization: Basic ZXlKaGJHY2lPaUpJVXpJMU5pSXNJbWxoZENJNk1UY3dOalV4T0RZd05Td2laWGh3SWpveE56TTRNRFUwTmpBMWZRLmV5SnBaQ0k2TVRVNExDSmxiV0ZwYkNJNkltMHVhaTV6WlhKblpXRnVkREZBWjI5dloyeGxiV0ZwYkM1amIyMGlMQ0oxYzJWeWJtRnRaU0k2SW1waFkyc2lMQ0pqYVhSNUlqb2lUV1YwY205d2IyeHBjeUlzSW1OdmRXNTBjbmtpT2lKVmJtbDBaV1FnUzJsdVoyUnZiU0lzSW1acGNuTjBibUZ0WlNJNklrcGhZMnNpTENKc1lYTjBibUZ0WlNJNklsTnRhWFJvSWl3aVkyOXVabWx5YldWa0lqb3hMQ0pwYm5OMGFYUjFkR2x2YmlJNklrRkRUVVVnYkdGaWN5SXNJbVJsY0dGeWRHMWxiblFpT2lKWFpXbHlaQ0JEYUdWdGFXTmhiSE1nUkdWd1lYSjBiV1Z1ZENJc0ltRmtiV2x1YVhOMGNtRjBiM0lpT201MWJHd3NJbUZqZEdsMlpTSTZiblZzYkN3aWRtbGxkMTl6Y0dWamFXVnpJam9pVkhKMVpTSXNJbUZ3YVY5aFkyTmxjM05mYzJWdWRHVnlhV05oSWpvaVZISjFaU0lzSW5acFpYZGZjM0JsWTJsbGN6RWlPaUpVY25WbElpd2laR1ZzWlhSbFgzTjBjbUZwYm5NaU9pSlVjblZsSWl3aVkyaGhibWRsWDJGemMyVnRZbXg1WDNOMFlYUjFjeUk2SWxSeWRXVWlMQ0pqYUdGdVoyVmZZWE56WlcxaWJIbGZjbVZzWldGelpWOWtZWFJsSWpvaVZISjFaU0lzSW1Ob1lXNW5aVjl6ZEhKaGFXNWZiM2R1WlhJaU9pSlVjblZsSWl3aWRYQnNiMkZrWDNKbFlXUnpNU0k2SWxSeWRXVWlmUS5MY1hsaWFTdDFuN1pBZmxXOHFFaF81NW9JZEtyREt4Vl91dl8tUmFQUEJVOicn" "https://enterobase.warwick.ac.uk/api/v2.0/senterica/cgMLST_v2/sts?limit=50&offset=0&scheme=cgMLST_v2"

    with lzma.open(out_file, 'w') as out_fh:
        sts = []
        offset: int = 0
        while offset < safety_valve:
            status, batch = fetch_hiercc_batch(key, offset, limit)
            if batch is None or not batch or not batch["STs"]:
                break
            sts.extend(batch["STs"])
            offset += limit
            print(f"{datetime.now()},{offset},{len(sts)}", file=sys.stderr)
        out_fh.write(json.dumps(sts).encode('utf-8'))
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
