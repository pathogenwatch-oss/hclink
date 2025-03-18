import lzma
import sqlite3
import struct
from pathlib import Path
from typing import BinaryIO, Callable, Iterator

from pyroaring import BitMap64


def write_bitmap_to_filehandle(filehandle: BinaryIO, encoder: Callable[[BitMap64], bytes], bitmap: BitMap64) -> None:
    serialized: bytes = encoder(bitmap)
    filehandle.write(struct.pack('<I', len(serialized)))  # Write length of serialized data
    filehandle.write(serialized)


def read_bitmaps(db: Path) -> Iterator[bytes]:
    with lzma.open(db, "rb") as f:
        while True:
            try:
                size_data: bytes = f.read(4)
                if not size_data:
                    break  # End of file reached
                size: int = struct.unpack('<I', size_data)[0]
                serialized: bytes = f.read(size)
                yield serialized
            except EOFError:
                break  # End of file reached


def read_st_info(st_db: Path) -> Iterator[tuple[str, list[str]]]:
    with lzma.open(st_db, "rt") as st_db_fh:
        for line in st_db_fh.readlines():
            info: list[str] = line.strip().split(",")
            yield info[0], info[1:]


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
    result: tuple[int] | None = next(cursor.execute("SELECT code FROM alleles WHERE checksum =? AND position =?", (st, position)), None)
    return result[0] if result is not None else None