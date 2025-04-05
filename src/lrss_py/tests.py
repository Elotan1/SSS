"""Module of Shamirs Sharing and Leakage Resilient Sharing Tests."""
import unittest
from shamir import ShamirSS, LRSS

# Define class to test the program
class TestShamir(unittest.TestCase):
    """Provides test cases for Shamirs Sharing and Leakage Resilient Sharing"""

    # Functions to test Shamir Secret Sharing Scheme:
    def test_one_digit_share_shamirs(self):
        """Test case for single digit shares"""
        field = ShamirSS(2, 3, rep="int")
        deconstruct = field.deconstruct(0x5, 2, 3)
        reconstruct = field.reconstruct(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_two_digit_share_shamirs(self):
        """Test case for two digit shares"""
        field = ShamirSS(2, 4, rep="int")
        deconstruct = field.deconstruct(0x5, 10, 12)
        reconstruct = field.reconstruct(deconstruct, 10, 12)
        self.assertEqual(5, reconstruct)

    def test_three_digit_share_shamirs(self):
        """Test case for three digit shares"""
        field = ShamirSS(2, 7, rep="int")
        deconstruct = field.deconstruct(0x5, 100, 120)
        reconstruct = field.reconstruct(deconstruct, 100, 120)
        self.assertEqual(5, reconstruct)

    def test_large_field_shamirs(self):
        """Test case for large finite fields"""
        field = ShamirSS(2, 100, rep="int")
        deconstruct = field.deconstruct(0x5, 2, 3)
        reconstruct = field.reconstruct(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_large_field_and_digit_shamirs(self):
        """Test case for large field and three digit shares"""
        field = ShamirSS(2, 100, rep="int")
        deconstruct = field.deconstruct(0x5, 290, 390)
        reconstruct = field.reconstruct(deconstruct, 290, 390)
        self.assertEqual(5, reconstruct)

# # Functions to test Leakage Resilient Secret Sharing Scheme:

    def test_one_digit_leakage_resilience(self):
        """Test case for single digit shares"""
        field = LRSS(2,5,rep="int")
        deconstruct = field.lr_share(0x5, 2, 3, 32, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_two_digit_leakage_resilience(self):
        """Test case for two digit shares"""
        field = LRSS(2,6,rep="int")
        deconstruct = field.lr_share(0x5, 10, 12, 32, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 10, 12)
        self.assertEqual(5, reconstruct)

    def test_three_digit_leakage_resilience(self):
        """Test case for three digit shares"""
        field = LRSS(2,7,rep="int")
        deconstruct = field.lr_share(0x5, 100, 120, 32, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 100, 120)
        self.assertEqual(5, reconstruct)

    def test_large_leakage(self):
        """Test case for large leakage budget"""
        field = LRSS(2,5,rep="int")
        deconstruct = field.lr_share(0x5, 2, 3, 300, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_large_error(self):
        """Test case for large leakage error"""
        field = LRSS(2,5,rep="int")
        deconstruct = field.lr_share(0x5, 2, 3, 32, 0.5)
        reconstruct = field.lr_share_rec(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_large_field_lrss(self):
        """Test case for large finite field"""
        field = LRSS(2,100,rep="int")
        deconstruct = field.lr_share(0x5, 2, 3, 32, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 2, 3)
        self.assertEqual(5, reconstruct)

    def test_large_field_and_digit_lrss(self):
        """Test case for large field and three digit shares"""
        field = LRSS(2,100,rep="int")
        deconstruct = field.lr_share(0x5, 100, 120, 32, 0.02)
        reconstruct = field.lr_share_rec(deconstruct, 100, 120)
        self.assertEqual(5, reconstruct)


if __name__ == '__main__':
    unittest.main()
