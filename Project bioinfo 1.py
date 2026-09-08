

def read_fasta(filename):
    sequences = {}
    with open(filename, "r") as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq)
                seq_id = line[1:]  # remove ">"
                seq = []
            else:
                seq.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq)
    return sequences

from Bio import SeqIO

for record in SeqIO.parse("sequence.fasta","fasta-pearson"):
   print(record.id)  
   print(record.seq)  


for record in SeqIO.parse("sequence.fasta","fasta-pearson"):
    seq = record.seq
    gc_content = 100 * float(seq.count("G") + seq.count("C")) / len(seq)
    print(f"{record.id} : GC content = {gc_content:.2f}%")

from Bio import SeqIO

motif = "ATG" #start codon
for record in SeqIO.parse("sequence.fasta","fasta-pearson"):
    seq = str(record.seq)
    positions = []
    for i in range(len(seq)  - len(motif) + 1):
        if seq[i:i+len(motif)] == motif:
            positions.append(i+1)  #+1 for biological indexing 
    print(f"{record.id} : motif '{motif}' found at positions {positions}")

        
from Bio import SeqIO

for record in SeqIO.parse("sequence.fasta", "fasta-pearson"):
    dna_seq = record.seq
    print(f"{record.id}:")
    print(f"DNA: {dna_seq}")
    protein_seq = dna_seq.translate()
    print(f"protein:{protein_seq}")


from Bio import SeqIO

# Read sequences
records = list(SeqIO.parse("sequence.fasta", "fasta-pearson"))

# Filter: keep only sequences longer than 100 bases
filtered_records = [rec for rec in records if len(rec.seq) > 100]

# Write filtered sequences to new FASTA file
SeqIO.write(filtered_records, "filtered.fasta", "fasta")
print("Filtered sequences saved to filtered.fasta")

from Bio import SeqIO
# Calculate average length 
lengths = [len(rec.seq) for rec in records]

print("Average length:", sum(lengths)/len(lengths))

reference = records[0]

mutations = []
for i, seq in enumerate(records[1:], start=1):
    for j in range(len(reference.seq)):
        if seq.seq[j] != reference.seq[j]:
            mutations.append((i, j, reference.seq[j], seq.seq[j]))


# gc value 
gc_Values = []
for record in records:
    seq = record.seq
    gc_content = 100 * float(seq.count("G") + seq.count("C")) / len(seq)
    gc_Values.append(gc_content)

print("GC Values:", gc_Values)

print("Sample mutations:", mutations[:10])

import matplotlib.pyplot as plt

plt.hist(gc_Values)
plt.title("GC Content Distribution")
plt.xlabel("GC %")
plt.ylabel("Frequency")
plt.show()