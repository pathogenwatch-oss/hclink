import lzma
import multiprocessing
import sys
from functools import partial
from pathlib import Path
from typing import Any, Iterator

from bitarray import bitarray
from bitarray.util import count_and, count_xor, deserialize, sc_decode


def compare_profiles(
        profile_a: tuple[bitarray, bitarray],
        profile_b: tuple[bitarray, bitarray]
) -> tuple[int, int, int, int]:
    """
    Compares two profiles and returns the distance between them.
    The distance is the number of bits that are different between the two profiles, excluding the gap positions.
    """
    shared_gaps: int = count_and(profile_a[1], profile_b[1])
    profileA_gaps: int = profile_a[1].count() - shared_gaps
    profileB_gaps: int = profile_b[1].count() - shared_gaps
    gap_adjust: int = profileA_gaps + profileB_gaps
    distance = count_xor(profile_a[0], profile_b[0])
    return int((distance - gap_adjust) / 2), profileA_gaps, profileB_gaps, shared_gaps


def comparison(
        query_profile: tuple[bitarray, bitarray],
        reference,
):
    profile = sc_decode(reference[0])
    distance, gaps_a, gaps_b, gaps_both = compare_profiles(query_profile, (profile, reference[1]))
    return distance, gaps_a, gaps_b, gaps_both, reference[2]


def calculate_hiercc_distance(distance, query_gaps, reference_gaps, shared_gaps, profile_size) -> float:
    # Distance bigger than profile size (i.e. failed comparison) then just return profile size
    if distance >= profile_size:
        return profile_size
    if distance == 0 and query_gaps == 0 and reference_gaps == 0:
        cc_distance = 0.0
    else:
        query_core = (profile_size - query_gaps - shared_gaps) - 0.03 * profile_size
        # if query core is > s use equation 1, else use equation 2
        common_core = profile_size - query_gaps - reference_gaps - shared_gaps
        if common_core >= query_core:
            cc_distance = (profile_size * distance) / common_core + 0.5
        else:
            cc_distance = ((profile_size * (distance + query_core - common_core)) / query_core) + 0.5
    return cc_distance


def infer_hiercc_code(hier_cc_distance: float,
                      hiercc_thresholds: list[int],
                      profile: list[str],
                      prepend: str) -> list[tuple[str, str]]:
    inferred: list[tuple[str, str]] = []
    if not profile:
        profile = [""] * len(hiercc_thresholds)
    if len(profile) != len(hiercc_thresholds):
        raise ValueError(f"The profile length is not the same as the thresholds list ({len(profile)} v{len(hiercc_thresholds)})")
    for index, threshold in enumerate(hiercc_thresholds):
        inferred.append((f"{prepend}{threshold}", (profile[index] if hier_cc_distance <= threshold else "")))
    return inferred


def read_profiles(lengths_db: Path, profile_db: Path) -> Iterator[bytes]:
    with open(lengths_db, "r") as lengths_fh, lzma.open(profile_db, "rb") as reference_profiles_fh:
        for length_str in lengths_fh.readlines():
            encoded_profile = reference_profiles_fh.read(int(length_str.strip()))
            yield encoded_profile


def read_profiles_list(lengths_db, profile_db) -> list[bytes]:
    profiles_list: list[bytes] = []
    with open(lengths_db, "r") as lengths_fh, lzma.open(profile_db, "rb") as reference_profiles_fh:
        for index, length_str in enumerate(lengths_fh.readlines()):
            encoded_profile = reference_profiles_fh.read(int(length_str.strip()))
            profiles_list.append(encoded_profile)
    return profiles_list


def read_gap_profiles(gap_db: Path, gap_lengths_db, num_families) -> Iterator[bitarray]:
    with open(gap_lengths_db, "r") as gap_lengths_fh, lzma.open(gap_db, "rb") as gap_profiles_fh:
        for index, length_str in enumerate(gap_lengths_fh.readlines()):
            encoded_profile = gap_profiles_fh.read(int(length_str.strip()))
            yield deserialize(encoded_profile)
            if len(deserialize(encoded_profile)) != num_families:
                raise Exception(f"Profile {index} bitarray length {len(deserialize(encoded_profile))}!= {num_families}")


def read_gap_profiles_list(gap_db, gap_lengths_db, num_families) -> list[bitarray]:
    gap_profiles: list[bitarray] = []
    with open(gap_lengths_db, "r") as gap_lengths_fh, lzma.open(gap_db, "rb") as gap_profiles_fh:
        for index, length_str in enumerate(gap_lengths_fh.readlines()):
            encoded_profile = gap_profiles_fh.read(int(length_str.strip()))
            gap_profiles.append(deserialize(encoded_profile))
            if len(gap_profiles[index]) != num_families:
                raise Exception(
                    f"Profile {index} bitarray length {len(gap_profiles[index])} != {num_families}")
    return gap_profiles


def read_st_info(st_db) -> Iterator[tuple[str, list[str]]]:
    with open(st_db, "r") as st_db_fh:
        for line in st_db_fh.readlines():
            info = line.strip().split(",")
            yield info[0], info[1:]


def imap_search(gap_profiles: Iterator[bitarray],
                sts: Iterator[tuple[str, list[str]]],
                profiles: Iterator[bytes],
                max_gaps: int,
                query_profile,
                num_cpu: int = 1,
                chunksize: int = 10000,
                threshold: int = sys.maxsize
                ) -> dict[str, Any]:

    if query_profile[1].count() >= max_gaps:
        return {
            "st": ("",[]),
            "distance": threshold,
            "gaps_a": query_profile[1].count(),
            "gaps_b": -1,
            "gaps_both": query_profile[1].count(),
        }
    lowest_distance: int = threshold
    best_hits: list[tuple[str, int, int, int, int]] = []
    query_comparison = partial(comparison, query_profile)

    print(f"Using {num_cpu} CPU cores", file=sys.stderr)
    with multiprocessing.Pool(num_cpu) as pool:
        for result in pool.imap_unordered(query_comparison,
                                          zip(profiles,
                                              gap_profiles,
                                              sts),
                                          chunksize=chunksize):
            total_gaps = result[1] + result[2] + result[3]
            if total_gaps >= max_gaps:
                continue
            if result[0] < lowest_distance:
                best_hits.append(result)
                lowest_distance = result[0]
            elif result[0] == lowest_distance:
                if total_gaps < best_hits[-1][1] + best_hits[-1][2] + best_hits[-1][3]:
                    best_hits.append(result)
    if not best_hits:
        return {
            "st": ("",[]),
            "distance": threshold,
            "gaps_a": query_profile[1].count(),
            "gaps_b": -1,
            "gaps_both": query_profile[1].count(),
        }
    return {
        "st": best_hits[-1][4],
        "distance": best_hits[-1][0],
        "gaps_a": best_hits[-1][1],
        "gaps_b": best_hits[-1][2],
        "gaps_both": best_hits[-1][3]}
