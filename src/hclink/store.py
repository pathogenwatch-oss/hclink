import csv
import gzip
import lzma
import sqlite3
import sys
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
from usearch.compiled import ScalarKind
from usearch.index import Index

from hclink.search import count_zeros, get_compiled_metric


def read_st_info(st_db: Path, st: str) -> tuple[str, list[str]]:
    with lzma.open(st_db, "rt") as st_db_fh:
        for line in st_db_fh.readlines():
            info: list[str] = line.strip().split(",")
            if st == info[0]:
                return info[0], info[1:]
        else:
            raise ValueError(f"No ST information found for ST code: {st}")


def initialise_db(dbfile: str) -> sqlite3.Connection:
    conn: sqlite3.Connection = sqlite3.connect(dbfile)
    cursor: sqlite3.Cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA synchronous=OFF")

    cursor.execute('CREATE TABLE IF NOT EXISTS alleles(checksum TEXT, position INTEGER, code INTEGER)')
    cursor.execute("DELETE FROM alleles")
    cursor.execute("DROP INDEX IF EXISTS idx_checksum")

    cursor.close()
    conn.commit()
    return conn


def connect_db(db_file: str) -> sqlite3.Connection:
    return sqlite3.connect(db_file)


def finalise_db(db: sqlite3.Connection) -> None:
    cursor: sqlite3.Cursor = db.cursor()
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_checksum ON alleles(checksum, position)')
    cursor.execute('VACUUM')
    cursor.execute('ANALYZE')
    cursor.execute('PRAGMA optimize')
    cursor.execute('PRAGMA query_only = ON')  # Set to read-only mode
    db.commit()
    cursor.close()


def lookup_st(cursor: sqlite3.Cursor, st: str, position: int) -> int | None:
    result: tuple[int] | None = next(
        cursor.execute("SELECT code FROM alleles WHERE checksum =? AND position =?", (st, position)), None)
    return result[0] if result is not None else None


def create_index(ndim: int, max_gaps) -> Index:
    return Index(
        ndim=ndim,
        dtype=ScalarKind.F32,
        metric=get_compiled_metric(max_gaps)
    )


def restore_index(index_file: Path, max_gaps) -> Index:
    index = Index.restore(index_file, view=True)
    index.metric = get_compiled_metric(max_gaps)
    return index


def batch_read_profiles(profiles_file: Path, max_gaps=5000, batch_size: int = 10000, limit: int = -1) -> Iterable[
    Tuple[np.ndarray, np.ndarray]]:
    with gzip.open(profiles_file, "rt") as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  # Skip header

        batch_keys = []
        batch_vectors = []
        print(f"Using max_gaps = {max_gaps}", file=sys.stderr)
        for row, line in enumerate(reader):
            if limit != -1 and row >= limit:
                break

            key = int(line[0])
            vector = np.array(list(map(float, line[1:])), dtype=np.float32)
            if count_zeros(vector) >= max_gaps:  # If max_gaps is set and profile has too many gaps, skip it
                # print(f"Skipping profile {key} with {count_zeros(vector)} gaps", file=sys.stderr)
                continue
            batch_keys.append(key)
            batch_vectors.append(vector)

            if len(batch_keys) == batch_size:
                # print(f"Read {len(batch_keys)} profiles", file=sys.stderr)
                yield np.array(batch_keys), np.array(batch_vectors)
                batch_keys = []
                batch_vectors = []

        # Yield any remaining profiles
        if batch_keys:
            yield np.array(batch_keys), np.array(batch_vectors)
