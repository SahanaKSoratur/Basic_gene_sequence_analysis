# DNA → RNA → Protein Converter (with AUG start detection)

# Step 1: Codon table (RNA codons to amino acids)
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


# Step 2: Transcription (DNA to RNA)
def transcribe_dna(dna):
    return dna.upper().replace('T', 'U')


# Step 3: Translation (RNA to Protein, starts at first AUG)
def translate_rna(rna):
    protein = ""
    start = rna.find("AUG")
    if start == -1:
        return "No start codon (AUG) found."

    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        amino_acid = codon_table.get(codon, 'X')  # 'X' = unknown
        if amino_acid == '*':
            break  # stop at first stop codon
        protein += amino_acid

    return protein


# Step 4: Read DNA from FASTA
with open("sample_sequence.fasta", "r") as file:
    lines = file.readlines()
    dna_seq = "".join(line.strip() for line in lines if not line.startswith(">"))

# Step 5: Run transcription and translation
rna_seq = transcribe_dna(dna_seq)
protein_seq = translate_rna(rna_seq)

# Step 6: Print results
print("DNA Sequence:     ", dna_seq)
print("RNA Sequence:     ", rna_seq)
print("Protein Sequence: ", protein_seq)