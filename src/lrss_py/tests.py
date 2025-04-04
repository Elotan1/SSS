import unittest
import time
from shamirv2 import *
 
# Define class to test the program
class TestShamir(unittest.TestCase):
        
    # Functions to test Shamir Secret Sharing Scheme:

        def test_one_digit_share_shamirs(self):
            S = ShamirSS(2, 3, repr="int")
            deconstruct = S.deconstruct(0x5, 2, 3)
            reconstruct = S.reconstruct(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)

        def test_two_digit_share_shamirs(self):
            S = ShamirSS(2, 4, repr="int")
            deconstruct = S.deconstruct(0x5, 10, 12)
            reconstruct = S.reconstruct(deconstruct, 10, 12)
            self.assertEqual(5, reconstruct)

        def test_three_digit_share_shamirs(self):
            S = ShamirSS(2, 7, repr="int")
            deconstruct = S.deconstruct(0x5, 100, 120)
            reconstruct = S.reconstruct(deconstruct, 100, 120)
            self.assertEqual(5, reconstruct)

        def test_large_field_shamirs(self):
            S = ShamirSS(2, 100, repr="int")
            deconstruct = S.deconstruct(0x5, 2, 3)
            reconstruct = S.reconstruct(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)

        def test_large_field_and_digit_shamirs(self):
            S = ShamirSS(2, 100, repr="int")
            deconstruct = S.deconstruct(0x5, 290, 390)
            reconstruct = S.reconstruct(deconstruct, 290, 390)
            self.assertEqual(5, reconstruct)

    # Functions to test Leakage Resilient Secret Sharing Scheme:

        def test_one_digit_leakage_resilience(self):
            L = LRSS(2,5,repr="int")
            deconstruct = L.lr_share(0x5, 2, 3, 32, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)
        
        def test_two_digit_leakage_resilience(self):
            L = LRSS(2,6,repr="int")
            deconstruct = L.lr_share(0x5, 10, 12, 32, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 10, 12)
            self.assertEqual(5, reconstruct)

        def test_three_digit_leakage_resilience(self):
            L = LRSS(2,7,repr="int")
            deconstruct = L.lr_share(0x5, 100, 120, 32, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 100, 120)
            self.assertEqual(5, reconstruct)
            
        def test_large_leakage(self):
            L = LRSS(2,5,repr="int")
            deconstruct = L.lr_share(0x5, 2, 3, 300, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)

        def test_large_error(self):
            L = LRSS(2,5,repr="int")
            deconstruct = L.lr_share(0x5, 2, 3, 32, 0.5)
            reconstruct = L.lr_share_rec(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)

        def test_large_field_lrss(self):
            L = LRSS(2,100,repr="int")
            deconstruct = L.lr_share(0x5, 2, 3, 32, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 2, 3)
            self.assertEqual(5, reconstruct)

        def test_large_field_and_digit_lrss(self):
            L = LRSS(2,100,repr="int")
            deconstruct = L.lr_share(0x5, 100, 120, 32, 0.02)
            reconstruct = L.lr_share_rec(deconstruct, 100, 120)
            self.assertEqual(5, reconstruct)
    
 
if __name__ == '__main__':
    unittest.main()
