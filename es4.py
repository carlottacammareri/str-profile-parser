from Bio import SeqIO
import re
import csv

#CONFIGURATION
input_fasta = "concatamer.fasta" # Input FASTA file containing concatenated reads
output_csv = "str_profiles.csv" # Output CSV file for STR motif counts
str_motifs = ["AGAT", "AATG", "TCTA"]  ## STR motifs to search for

#FUNCTION TO COUNT STR MOTIFS
def count_str_motifs(seq, motifs):
    counts = {}
    for motif in motifs:
        # Regex pattern to find consecutive repeats of the motif
        pattern = f"(?:{motif})+"
        matches = re.findall(pattern, str(seq))
        # Count total repeats by dividing the match length by motif length
        counts[motif] = sum(len(m)//len(motif) for m in matches)
    return counts

#Process sequences
results = []
# Read sequences from the FASTA file
for record in SeqIO.parse(input_fasta, "fasta"):
    # Split each sequence into fragments using stretches of N (2 or more) as delimiters
    fragments = re.split("N{2,}", str(record.seq)) 
    for i, frag in enumerate(fragments):
        counts = count_str_motifs(frag, str_motifs)
        results.append({
            "fragment_id": f"{record.id}_frag{i+1}",
            **counts
        })

#Create CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["fragment_id"] + str_motifs)
    writer.writeheader()
    writer.writerows(results)

