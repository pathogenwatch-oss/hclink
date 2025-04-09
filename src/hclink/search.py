import math

import numpy as np
from numba import carray, cfunc, njit, types
from usearch.compiled import MetricKind, MetricSignature
from usearch.index import CompiledMetric, Match, Matches

signature = types.float32(
    types.CPointer(types.float32),
    types.CPointer(types.float32),
    types.uint64)


@njit
def debug_print(*args):
    print(*args)


def get_compiled_metric(max_gaps: int):
    return CompiledMetric(pointer=get_distance_score(max_gaps).address, kind=MetricKind.Unknown, signature=MetricSignature.ArrayArraySize)


def get_distance_score(max_gaps: int):

    @cfunc(signature)
    @njit
    def hiercc_distance(a_ptr, b_ptr, ndim):
        a_array = carray(a_ptr, ndim)
        b_array = carray(b_ptr, ndim)
        
        gaps_a = 0
        gaps_b = 0
        gaps_both = 0
        distance = 0
        
        for i in range(ndim):
            a_val = a_array[i]
            b_val = b_array[i]
            if a_val == b_val:
                gaps_both += (a_val == 0)
            else:
                distance += (a_val != 0 and b_val != 0)
                gaps_a += (a_val == 0)
                gaps_b += (b_val == 0)
        
        if gaps_both >= max_gaps or distance >= ndim:
            return np.float32(ndim)
        
        if distance == 0 and gaps_a == 0 and gaps_b == 0:
            return np.float32(0.0)
        
        ndim_f32 = np.float32(ndim)
        query_core = ndim_f32 - np.float32(gaps_a + gaps_both) - np.float32(0.03 * ndim)
        common_core = ndim_f32 - np.float32(gaps_a + gaps_b + gaps_both)
        
        if common_core >= query_core:
            if common_core == 0:
                return ndim_f32
            cc_distance = (ndim_f32 * np.float32(distance)) / common_core + np.float32(0.5)
        else:
            if query_core == 0:
                return ndim_f32
            cc_distance = ((ndim_f32 * np.float32(distance + query_core - common_core)) / query_core) + np.float32(0.5)
        
        return min(cc_distance, ndim_f32)

    return hiercc_distance

def infer_hiercc_code(hier_cc_distance: float,
                      hiercc_thresholds: list[int],
                      profile: list[str],
                      prepend: str) -> list[tuple[str, str]]:
    inferred: list[tuple[str, str]] = []
    if not profile:
        profile = [""] * len(hiercc_thresholds)
    if len(profile) != len(hiercc_thresholds):
        raise ValueError(
            f"The profile length is not the same as the thresholds list ({len(profile)} v{len(hiercc_thresholds)})")
    for index, threshold in enumerate(hiercc_thresholds):
        inferred.append((f"{prepend}{threshold}", (profile[index] if hier_cc_distance <= threshold else "")))
    return inferred


def select_best(matches: Matches) -> Match | None:
    best_distance = np.inf
    best_matches = []
    for index, distance in enumerate(matches.distances):
        if distance < best_distance:
            best_distance = distance
            best_matches = [index]
        elif distance == best_distance:
            best_matches.append(index)
        else:
            break
    if not best_matches:
        return None
    return matches[best_matches[0]]


def count_zeros(profile: np.ndarray) -> int:
    return profile.size - np.count_nonzero(profile)


def count_shared_zeros(a: np.ndarray, b: np.ndarray) -> int:
    return np.sum(np.logical_and(a == 0, b == 0)).item()
