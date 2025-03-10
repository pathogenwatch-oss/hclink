import lzma
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


# def read_bitmaps(db: Path) -> Iterator[BitMap64]:
#     with lzma.open(db, "rb") as f:
#         while True:
#             try:
#                 size_data = f.read(4)
#                 if not size_data:
#                     break  # End of file reached
#                 size = struct.unpack('<I', size_data)[0]
#                 serialized = f.read(size)
#                 yield BitMap64.deserialize(serialized)
#             except EOFError:
#                 break  # End of file reached


def read_st_info(st_db) -> Iterator[tuple[str, list[str]]]:
    with lzma.open(st_db, "rb") as st_db_fh:
        for line in st_db_fh.readlines():
            info = line.strip().decode("utf-8").split(",")
            yield info[0], info[1:]
