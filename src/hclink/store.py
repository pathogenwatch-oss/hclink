import lzma
import sqlite3
import struct
from pathlib import Path
from typing import Iterator

from bitarray import bitarray


def write_bitmap_to_filehandle(filehandle, encoder, bitmap: bitarray):
    serialized = encoder(bitmap)
    filehandle.write(struct.pack('<I', len(serialized)))  # Write length of serialized data
    filehandle.write(serialized)


def read_bitmaps(db: Path) -> Iterator[bytes]:
    with lzma.open(db, "rb") as f:
        while True:
            try:
                size_data = f.read(4)
                if not size_data:
                    break  # End of file reached
                size = struct.unpack('<I', size_data)[0]
                serialized = f.read(size)
                yield serialized
            except EOFError:
                break  # End of file reached


def read_st_info(st_db) -> Iterator[tuple[str, list[str]]]:
    with lzma.open(st_db, "rt") as st_db_fh:
        for line in st_db_fh.readlines():
            info = line.strip().split(",")
            yield info[0], info[1:]


def initialise_db(dbfile):
    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    
    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA synchronous=OFF")
    
    cursor.execute('CREATE TABLE IF NOT EXISTS alleles(checksum TEXT, position INTEGER, code INTEGER)')
    cursor.execute("DELETE FROM alleles")
    cursor.execute("DROP INDEX IF EXISTS idx_checksum")
    
    cursor.close()
    conn.commit()
    return conn


def connect_db(dbfile):
    return sqlite3.connect(dbfile)


def finalise_db(db):
    cursor = db.cursor()
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_checksum ON alleles(checksum, position)')
    cursor.execute('VACUUM')
    cursor.execute('ANALYZE')
    cursor.execute('PRAGMA optimize')
    cursor.execute('PRAGMA query_only = ON')  # Set to read-only mode
    db.commit()
    cursor.close()

def lookup_st(cursor, st, position) -> int:
    result = next(cursor.execute("SELECT code FROM alleles WHERE checksum =? AND position =?", (st, position)), None)
    return result[0] if result is not None else None
