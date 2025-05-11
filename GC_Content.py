# Simple GC content calculator
def calculate_gc(dna_sequence):
    g = dna_sequence.upper().count('G')
    c = dna_sequence.upper().count('C')
    gc_percent = (g + c) / len(dna_sequence) * 100
    return gc_percent

with open("sample_sequence.fasta", "r") as file:
    lines = file.readlines()
    dna = "".join(line.strip() for line in lines if not line.startswith(">"))

gc = calculate_gc(dna)
print(f"GC Content: {gc:.2f}%")
