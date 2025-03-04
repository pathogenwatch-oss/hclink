from typing import Iterable

from bitarray import bitarray
from bitarray.util import count_and, count_xor, sc_decode


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
