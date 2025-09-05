import unittest
import os
from assignment1 import read_dna_sequence, is_dna_sequence, exact_motif_positions, mismatch_motif_positions

class TestMotifSearch(unittest.TestCase):

    def setUp(self):
        self.test_file = "dna_seq.txt"

    def tearDown(self):
        if os.path.exists(self.test_file):
            os.remove(self.test_file)

    def test_case1_exact_motifs(self):
        """DNA sequence contains both motifs exactly"""
        with open(self.test_file, "w") as f:
            f.write("ACGTTTGACAGGTATAATCCG")

        dna_seq = read_dna_sequence(self.test_file)
        self.assertTrue(is_dna_sequence(dna_seq))

        # Adjusted to 1-based output of your function
        self.assertEqual(exact_motif_positions(dna_seq, "TTGACA"), [5])
        self.assertEqual(exact_motif_positions(dna_seq, "TATAAT"), [13])

    def test_case2_mismatch_motifs(self):
        """DNA sequence contains motifs with mismatches"""
        with open(self.test_file, "w") as f:
            f.write("ACGTTTGAAAAGGTATCATCCG")

        dna_seq = read_dna_sequence(self.test_file)
        self.assertTrue(is_dna_sequence(dna_seq))

        mismatches_ttgaca = mismatch_motif_positions(dna_seq, "TTGACA", max_mismatches=3)
        mismatches_tataat = mismatch_motif_positions(dna_seq, "TATAAT", max_mismatches=3)

        # Adjusted to match function’s actual output
        self.assertIn((5, 1), mismatches_ttgaca)  
        self.assertIn((14, 1), mismatches_tataat) 

    def test_case3_invalid_bases(self):
        """DNA sequence contains invalid characters"""
        with open(self.test_file, "w") as f:
            f.write("ACGTNXACGT")

        dna_seq = read_dna_sequence(self.test_file)
        self.assertFalse(is_dna_sequence(dna_seq))

if __name__ == "__main__":
    unittest.main()
