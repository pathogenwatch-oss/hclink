import ctypes
import json
import sqlite3
from functools import partial
from typing import Callable

import numpy as np

from hclink.search import get_distance_score
from hclink.store import connect_db, lookup_st
from hclink.build import convert_to_vector


def test_convert_to_profile():

    cursor: sqlite3.Cursor = connect_db(str("alleles.db")).cursor()
    lookup: Callable[[str, int], int] = partial(lookup_st, cursor)


    with open("foo.json", "r") as f:
        code = json.load(f)["code"]
        profile = convert_to_vector(code, lookup)
    with open("example_profile1.bin", "wb") as f:
        np.save(f, profile)



def test_hiercc_distance_identical_arrays():
    a = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    b = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    ndim = 3
    hiercc_distance = get_distance_score(0)
    result = hiercc_distance(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             ndim)
    assert result == 0.0, f"Expected distance 0.0 for identical arrays, but got {result}"

def test_hiercc_distance_ten_differences():
    a = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], dtype=np.float32)
    b = np.array([11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0], dtype=np.float32)
    ndim = 10
    hiercc_distance = get_distance_score(0)
    result = hiercc_distance(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             ndim)
    expected = 10.5  # 10 differences, no gaps, so (10 * 10) / 10 + 0.5
    assert abs(result - expected) < 1e-6, f"Expected distance {expected} for arrays with 10 differences, but got {result}"

def test_hiercc_distance_ten_differences_seven_common():
    a = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0], dtype=np.float32)
    b = np.array([11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0], dtype=np.float32)
    ndim = 17
    hiercc_distance = get_distance_score(0)
    result = hiercc_distance(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             ndim)
    expected = 10.5  # 10 differences, 7 common, no gaps, so (17 * 10) / 17 + 0.5
    assert abs(result - expected) < 1e-6, f"Expected distance {expected} for arrays with 10 differences and 7 common elements, but got {result}"

def test_hiercc_distance_with_gaps():
    a = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0, 10.0, 21.0, 22.0, 23.0, 0.0, 25.0, 26.0, 27.0], dtype=np.float32)
    b = np.array([11.0, 12.0, 13.0, 14.0, 0.0, 16.0, 17.0, 18.0, 0.0, 20.0, 21.0, 22.0, 23.0, 24.0, 0.0, 26.0, 27.0], dtype=np.float32)
    ndim = 17
    hiercc_distance = get_distance_score(0)
    result = hiercc_distance(a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                             ndim)
    
    # Calculate expected result
    distance = 8  # 8 differences (positions 1-4, 6-8, 10)
    gaps_a = 1  # One gap in 'a' not shared with 'b' (position 14)
    gaps_b = 2  # Two gaps in 'b' not shared with 'a' (positions 5 and 15)
    gaps_both = 1  # One shared gap (position 9)
    common_core = ndim - gaps_a - gaps_b - gaps_both
    expected = (float(ndim) * float(distance)) / common_core + 0.5
    print(f"Expected distance: {expected}")
    
    assert abs(result - expected) < 1e-6, f"Expected distance {expected} for arrays with gaps, but got {result}"
