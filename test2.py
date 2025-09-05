from assignment1 import read_dna_sequence, exact_motif_positions, mismatch_motif_positions, is_dna_sequence

# Load DNA sequence from file
dna_seq = "TTGACAGG GCTCGAT"

print("DNA Sequence:", dna_seq)
print("IS DNA valid?", is_dna_sequence(dna_seq))
print()

# ---- Test Case 1: Exact Motifs ----
expected_ttgaca = [1]
output_ttgaca = exact_motif_positions(dna_seq, "TTGACA")
print("Test Case 1 - TTGACA")
print("Expected:", expected_ttgaca, " Output:", output_ttgaca, " Match?", expected_ttgaca == output_ttgaca)

expected_tataat = []
output_tataat = exact_motif_positions(dna_seq, "TATAAT")
print("Test Case 1 - TATAAT")
print("Expected:", expected_tataat, " Output:", output_tataat, " Match?", expected_tataat == output_tataat)
print()

# ---- Test Case 2: Mismatches ----
dna_seq_mismatch = "TTGAAAAGGCTCGAT"  # differs by 1 base at position 6
expected_ttgaca_mismatches = [(1, 1)]
output_ttgaca_mismatches = mismatch_motif_positions(dna_seq_mismatch, "TTGACA")
print("Test Case 2 - TTGACA Mismatches")
print("Expected:", expected_ttgaca_mismatches, " Output:", output_ttgaca_mismatches, " Match?", expected_ttgaca_mismatches == output_ttgaca_mismatches)
print()
# ---- Test Case 3: Invalid DNA ----
invalid_seq = "ACGTNXACGT"
expected_validity = False
output_validity = is_dna_sequence(invalid_seq)
print("Test Case 3 - Invalid DNA")
print("Expected:", expected_validity, " Output:", output_validity, " Match?", expected_validity == output_validity)
