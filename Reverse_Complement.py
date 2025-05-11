## Program to find the reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    # Define base pairing rules
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Create the complement strand
    comp_seq = ''.join([complement.get(base.upper(), base) for base in dna_sequence])

    # Reverse it
    rev_comp_seq = comp_seq[::-1]

    return rev_comp_seq

# Read sequence from FASTA
with open("sample_sequence.fasta", "r") as file:
    lines = file.readlines()
    dna = "".join(line.strip() for line in lines if not line.startswith(">"))
# Get reverse complement
rev_comp = reverse_complement(dna)

# Show result
print("Original sequence:     ", dna)
print("Reverse complement:    ", rev_comp)