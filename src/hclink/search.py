import lzma
from pathlib import Path
from typing import Iterable, Iterator

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
        reference_profile,
        gap_profile,
        st: tuple[str, list[str]]
):
    profile = sc_decode(reference_profile)
    distance, gaps_a, gaps_b, gaps_both = compare_profiles(query_profile, (profile, gap_profile))
    return distance, gaps_a, gaps_b, gaps_both, st


def comparison2(
        query_profile: tuple[bitarray, bitarray],
        reference,
        # gap_profile,
        # st: tuple[str, list[str]]
):
    profile = sc_decode(reference[0])
    distance, gaps_a, gaps_b, gaps_both = compare_profiles(query_profile, (profile, reference[1]))
    return distance, gaps_a, gaps_b, gaps_both, reference[2]

def calculate_hiercc_distance(distance, query_gaps, reference_gaps, shared_gaps, profile_size) -> float:
    # s = shared genes - 0.03 * profile_size

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


def infer_hiercc_code(hier_cc_distance: float, hier_cc_map: Iterable[tuple[int, str]]) -> Iterable[tuple[str, str]]:
    for threshold, code in hier_cc_map:
        formatted_threshold = f"d{threshold}"
        if hier_cc_distance <= threshold:
            yield formatted_threshold, code
        else:
            yield formatted_threshold, ""


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
            if len(deserialize(encoded_profile))!= num_families:
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
