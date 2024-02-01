from unittest import TestCase

from bitarray import bitarray

from hclink import convert_to_profile


class Test(TestCase):
    def test_convert_to_profile(self):
        family_sizes = [3, 2, 1, 4, 2]
        profile = "_".join(["2", "", "x", "2", "x"])

        bits, gaps = convert_to_profile(
            profile,
            sum(family_sizes) + len(family_sizes),
            family_sizes
        )
        self.assertEqual(bitarray("01000000101000001"), bits)
        self.assertEqual(bitarray("01000"), gaps)
