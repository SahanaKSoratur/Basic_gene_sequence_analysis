import csv

# GC content calculator
def calculate_gc(dna_sequence):
    g = dna_sequence.upper().count('G')
    c = dna_sequence.upper().count('C')
    gc_percent = (g + c) / len(dna_sequence) * 100
    return gc_percent

#  Reverse complement function
def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(complement.get(base.upper(), base) for base in reversed(seq))

# PAM site scanner function (forward or reverse strand)
def scan_pam_sites(dna, strand_label, results):
    for i in range(len(dna) - 23):  # 20bp guide + 3bp PAM
        guide = dna[i:i+20]
        pam = dna[i+20:i+23]

        if pam[1:] == "GG":
            gc = calculate_gc(guide)
            if 40 <= gc <= 60:
                results.append({
                    "Position": i + 1,
                    "Strand": strand_label,
                    "Guide": guide,
                    "PAM": pam,
                    "GC%": round(gc, 2)
                })

# Step 1: Read DNA from FASTA
with open("sample_sequence.fasta", "r") as file:
    lines = file.readlines()
    dna = "".join(line.strip() for line in lines if not line.startswith(">")).upper()

# Step 2: Scan both strands
all_hits = []
scan_pam_sites(dna, "Forward", all_hits)
scan_pam_sites(reverse_complement(dna), "Reverse", all_hits)

# Step 3: Print and save to CSV
if not all_hits:
    print("No valid CRISPR sites found.")
else:
    print("CRISPR Hits (GC 40–60%):\n")
    print("Pos\tStrand\t\tGuide RNA\t\tPAM\tGC%")
    print("-" * 60)
    for hit in all_hits:
        print(f"{hit['Position']}\t{hit['Strand']:7}\t{hit['Guide']}\t{hit['PAM']}\t{hit['GC%']}")

    # Save to CSV
    with open("crispr_hits.csv", "w", newline='') as csvfile:
        fieldnames = ["Position", "Strand", "Guide", "PAM", "GC%"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_hits)

    print("\n Results saved to 'crispr_hits.csv'")
