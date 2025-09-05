
"""
BioInformatics Assignment 1a
Ritika Lama and Anna Nguyen
-------------------------
This program identifies common prokaryotic promoter motifs in a DNA sequence,
specifically the -35 (TTGACA) and -10 (TATAAT) consensus sequences, which
are key regulatory elements in prokaryotic transcription.

Flow of the code:
1. Reads a DNA sequence from a text file, removing whitespace and newlines and converting
   all characters to uppercase.
2. Validates that the sequence contains only valid nucleotides (A, C, G, T).
3. Searches for exact matches of the promoter motifs in the sequence.
4. Searches for near matches allowing up to 2 nucleotide mismatches.
5. Calculates match scores for near matches based on the proportion of matching
   nucleotides(mismatches/6).
6. Highlights exact matches in green and mismatches in red in the output terminal in a readable format.
7. Reports motif positions and match scores to provide a clear visualization
   of potential promoter regions.

This program combines biological knowledge of promoter sequences with
computational string analysis to support motif detection and sequence
interpretation.
"""

"""
    This function reads the DNA sequence from a text file.

    Parameters:
        text (str): Path to the text file containing the DNA sequence.

    Returns:
        str: DNA sequence in uppercase with whitespace and newlines removed.
             Returns an empty string if the file is not found.
"""
def read_dna_sequence(text):
    try:
        with open(text, "r") as file:
            sequence = file.read()
        dna_seq = sequence.replace("\n", "").replace(" ", "").upper()
        print(f"DNA Sequence: {dna_seq}")
        return dna_seq
    except FileNotFoundError:
        print("File not found. Please make sure that the file is correct andexists in the current directory.")
        return ""

"""
    This function validates that the sequence contains only valid DNA bases(A,T,G,C).

    Parameters:
        dna_seq (str): DNA sequence to validate that was read from the text file.

    Returns:
        boolean: False if the sequence contains any nucleotide except A, C, G, T. False otherwise.
"""

def is_dna_sequence(dna_seq):
    base = "ACGT"
    for nucleotide in dna_seq:
        if nucleotide not in base:
            return False
    return True

"""
    This function calculates the number of mismatched nucleotides between two motifs.

    Parameters:
        seq1 (str): common motif.
        seq2 (str): Second motif extracted from DNA sequence.
    Returns:
        int: Number of positions where nucleotides differ.
"""
def calculate_mismatches(seq1, seq2):
    #return error if the lengths of the two motifs are not equal
    if len(seq1) != len(seq2):
        print("Strings must be the same length.")

    #set total mismatches to 0 and loop through the sequences to count mismatches
    total = 0
    #zip function is used to pair up the nucleotides from both sequencces
    for a, b in zip(seq1, seq2):
        #if the nucleotides are different, increment the mismatch count
        if a != b:
            total += 1
    return total

"""
    This function finds all exact positions of the common promoter motif in a DNA sequence.

    Parameters:
        dna_seq (str): DNA sequence to search.
        motif (str): Motif sequence to locate.

    Returns:
        list: Positions of all exact motif matches.
"""
def exact_motif_positions(dna_seq, motif):
    positions = []
    start = 0
    while True:
        # start shows the index of the first character of the motif in the dna sequence
        start = dna_seq.find(motif, start)
        # if the position is not found then exit the loop
        if start == -1:
            break
        # if found then append the position to the list
        positions.append(start + 1)  # +1 for 1-based indexing
        start += 1  # Move to the next character for overlapping matches
    return positions

"""
    This function finds positions of motif matches with mismatches allowed.

    Parameters:
        dna_seq (str): DNA sequence to search.
        motif (str): Motif sequence to locate.
        max_mismatches (int): Maximum number of mismatches allowed.

    Returns:
        list: Tuples of (position, mismatch_count) for motif matches
              that differ by up to max_mismatches nucleotides.
"""
def mismatch_motif_positions(dna_seq, motif, max_mismatches=2):
    mismatch_positions = []
    #check for parts in the dna sequence that have the matching first 2 nucleotides with the motif
    for i in range(len(dna_seq) - len(motif) + 1):
        #get that part of the dna sequence
        segment = dna_seq[i:i+len(motif)]
        #if the first 2 nucleotides match, check for mismatches
        if segment[:2] == motif[:2]:
            #count mismatches
            mismatches = calculate_mismatches(segment, motif)
            #if there are mismatches within the allowed limit, store the position and mismatch count
            if 0 < mismatches <= max_mismatches:
                mismatch_positions.append((i + 1, mismatches))
    return mismatch_positions
    
"""
    This function highlights detected motifs in the DNA sequence with brackets.

    Parameters:
        dna_seq (str): Original DNA sequence.
        motif (str): Motif sequence to highlight.
        positions (list): Positions (1-based indexing) where motifs occur.

    Returns:
        list: Modified DNA sequences with motifs enclosed in brackets.
"""
def highlighted_sequence(dna_seq, motif, positions):
    GREEN = "\033[92m"
    RESET = "\033[0m"
    highlighted_seqs = []
    for pos in positions:
        highlighted_seq = (dna_seq[:pos-1] + GREEN + motif + RESET + dna_seq[pos-1+len(motif):])
        highlighted_seqs.append(highlighted_seq)
    return highlighted_seqs

"""
    Search for exact and mismatch occurrences of motifs in a DNA sequence.

    Parameters:
        dna_seq (str): DNA sequence to search.
        motifs (list): List of motifs to check.

    Output:
        Prints positions of exact matches, mismatch matches, and highlights
        sequences with detected motifs for each motif in the list.
"""
def motif_checker(dna_seq, motifs):
    RESET = "\033[0m"
    RED = "\033[31m"
    for motif in motifs:
        # Find exact positions of the motif
        exact_positions = exact_motif_positions(dna_seq, motif)
        if exact_positions:
            print(f"Motif '{motif}' is found at {exact_positions} position.")
            highlighted_seqs = highlighted_sequence(dna_seq, motif, exact_positions)
            for seq in highlighted_seqs:
                print(f"Sequence with motif highlighted in green: {seq}")
        else:
            print(f"Motif '{motif}' not found in the sequence.")
        
        # Find potential mismatches of the motif
        mismatch_positions = mismatch_motif_positions(dna_seq, motif)
        if mismatch_positions:
            match_scores = {}
            for pos, mismatch_count in mismatch_positions:
                matches = len(motif) - mismatch_count
                score_percent = matches / len(motif)
                match_scores[pos] = round(score_percent, 2)

            print()
            print("There are some potential mismatches found in the DNA sequence.")
            print(f"Overall match scores for motif '{motif}'(given in position: match score format): {match_scores}")
            print()
            
            for pos, mismatch_count in mismatch_positions:
                segment = dna_seq[pos-1:pos-1+len(motif)]
                highlighted_seq = (dna_seq[:pos-1] + RED + segment + RESET + dna_seq[pos-1+len(motif):])
                print(f"Sequence with potential mismatch highlighted in green: {highlighted_seq}")
            print()
        else:
            print(f"No mismatches found")

    
"""
    This is the main program workflow.

    The steps followed in this function by calling other functions are as follows:
        1. Read DNA sequence from file (dna_sequence.txt).
        2. Validate that the sequence is a proper DNA sequence.
        3. Run motif search for exact and mismatch matches.
        4. Print results to the console.
"""
def main():
    motifs = ["TTGACA", "TATAAT"]
    text = "dna_sequence.txt"
    dna_sequence = read_dna_sequence(text)
    if dna_sequence and is_dna_sequence(dna_sequence):
        print("The sequence is a valid DNA sequence.\n")
        motif_checker(dna_sequence, motifs)
    else:
        print("The sequence is not a valid DNA sequence.")

if __name__ == "__main__":
    main()