import multiprocessing
import sys
from functools import partial
from typing import Any, Iterator

from numba import float32, int32, jit, types
from pyroaring import BitMap64

def comparison(
        query_profile: tuple[BitMap64, BitMap64],
        profile_size: int,
        reference: tuple[bytes, bytes, str],
) -> tuple[float, int, int, int, int, str]:
    reference_profile = BitMap64.deserialize(reference[0])
    reference_gaps = BitMap64.deserialize(reference[1])

    shared_gaps: int = query_profile[1].intersection_cardinality(reference_gaps)
    profile_a_gaps: int = query_profile[1].get_statistics()["cardinality"]
    profile_b_gaps: int = reference_gaps.get_statistics()["cardinality"]
    bit_distance: int = query_profile[0].symmetric_difference_cardinality(reference_profile)
    raw_distance, gaps_a, gaps_b = compiled_compare(bit_distance, shared_gaps, profile_a_gaps, profile_b_gaps)
    hiercc_distance = calculate_hiercc_distance(raw_distance, gaps_a, gaps_b, shared_gaps, profile_size)
    return hiercc_distance, raw_distance, gaps_a, gaps_b, shared_gaps, reference[2]

@jit(types.Tuple((int32, int32, int32))(int32, int32, int32, int32), nopython=True, nogil=True)
def compiled_compare(distance:int, shared_gaps: int, profile_a_gaps: int, profile_b_gaps: int) -> tuple[int, int, int]:
    gaps_a: int = profile_a_gaps - shared_gaps
    gaps_b: int = profile_b_gaps - shared_gaps
    gap_adjust: int = gaps_a + gaps_b
    return int((distance - gap_adjust) / 2), gaps_a, gaps_b

@jit(float32(int32, int32, int32, int32, int32), nopython=True, nogil=True)
def calculate_hiercc_distance(distance: int, query_gaps: int, reference_gaps: int, shared_gaps: int,
                              profile_size: int) -> float:
    # Distance bigger than profile size (i.e. failed comparison) then just return profile size
    if distance >= profile_size:
        return float(profile_size)
    if distance == 0 and query_gaps == 0 and reference_gaps == 0:
        cc_distance: float = 0.0
    else:
        query_core: float = float(profile_size - query_gaps - shared_gaps) - 0.03 * float(profile_size)
        # if query core is > s use equation 1, else use equation 2
        common_core: float = float(profile_size - query_gaps - reference_gaps - shared_gaps)
        if common_core >= query_core:
            cc_distance: float = (float(profile_size) * float(distance)) / common_core + 0.5
        else:
            cc_distance: float = ((float(profile_size) * float(distance + query_core - common_core)) / query_core) + 0.5
    return cc_distance

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

def imap_search(gap_profiles: Iterator[bytes],
                sts: Iterator[tuple[str, list[str]]],
                profiles: Iterator[bytes],
                max_gaps: int,
                query_profile: tuple[BitMap64, BitMap64],
                profile_size: int,
                num_cpu: int = 1,
                chunksize: int = 10000,
                threshold: int = sys.maxsize
                ) -> dict[str, Any]:
    if query_profile[1].get_statistics()["cardinality"] >= max_gaps:
        return {
            "st": ("", []),
            "distance": threshold,
            "hiercc_distance": threshold,
            "gaps_a": query_profile[1].get_statistics()["cardinality"],
            "gaps_b": -1,
            "gaps_both": query_profile[1].get_statistics()["cardinality"],
        }
    lowest_distance: int = threshold
    best_hits: list[tuple[str, int, int, int, int]] = []

    query_comparison = partial(comparison, query_profile, profile_size)

    print(f"Using {num_cpu} CPU cores", file=sys.stderr)
    with multiprocessing.Pool(num_cpu) as pool:
        for result in pool.imap_unordered(query_comparison,
                                          zip(profiles,
                                              gap_profiles,
                                              sts),
                                          chunksize=chunksize):
            total_gaps: int = result[2] + result[3] + result[4]
            if total_gaps >= max_gaps:
                continue
            if result[0] < lowest_distance:
                best_hits.append(result)
                lowest_distance = result[0]
            elif result[0] == lowest_distance:
                if total_gaps < best_hits[-1][2] + best_hits[-1][3] + best_hits[-1][4]:
                    best_hits.append(result)
    if not best_hits:
        return {
            "st": ("", []),
            "distance": threshold,
            "hiercc_distance": threshold,
            "gaps_a": query_profile[1].get_statistics()["cardinality"],
            "gaps_b": -1,
            "gaps_both": query_profile[1].get_statistics()["cardinality"],
        }
    return {
        "st": best_hits[-1][5],
        "hiercc_distance": best_hits[-1][0],
        "distance": best_hits[-1][1],
        "gaps_a": best_hits[-1][2],
        "gaps_b": best_hits[-1][3],
        "gaps_both": best_hits[-1][4]
    }
