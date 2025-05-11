# ORF Finder: Finds all ORFs in 3 reading frames

codon_table = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C', 'UGC': 'C',
    'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T',
    'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N',
    'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V',
    'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G',
    'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'UAA': '*',
    'UAG': '*', 'UGA': '*'  # Stop codons
}

stop_codons = ['UAA', 'UAG', 'UGA']

# Step 1: Transcribe DNA to RNA
def transcribe_dna(dna):
    return dna.upper().replace('T', 'U')

# Step 2: Find ORFs in a single reading frame
def find_orfs(rna, frame):
    orfs = []
    i = frame
    while i < len(rna) - 2:
        codon = rna[i:i+3]
        if codon == 'AUG':  # Found a start codon
            protein = ''
            start = i
            for j in range(i, len(rna)-2, 3):
                codon = rna[j:j+3]
                amino = codon_table.get(codon, 'X')
                if codon in stop_codons:
                    orfs.append((start, j+3, protein))  # Save (start, end, protein)
                    break
                protein += amino
            i = j  # Skip to end of this ORF
        else:
            i += 3
    return orfs

# Step 3: Read FASTA and run all frames
with open("sample_sequence.fasta", "r") as file:
    lines = file.readlines()
    dna_seq = "".join(line.strip() for line in lines if not line.startswith(">"))

rna_seq = transcribe_dna(dna_seq)

# Step 4: Run in all 3 reading frames
for frame in range(3):
    print(f"\nReading Frame {frame+1}:")
    found_orfs = find_orfs(rna_seq, frame)
    if not found_orfs:
        print("  No ORFs found.")
    for start, end, protein in found_orfs:
        print(f"  ORF from {start} to {end} → Protein: {protein}")