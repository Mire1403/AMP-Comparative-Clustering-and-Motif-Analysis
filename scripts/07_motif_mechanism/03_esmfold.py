import torch
import esm
from Bio import SeqIO
import os

FASTA_FILE = "results/statistics/09_esmfold_input/esm_input.fasta"
OUTDIR = "results/statistics/09_esmfold_input/pdbs"

print("Loading ESMFold model...")
model = esm.pretrained.esmfold_v1()
model = model.eval().cpu()

# reduce memory usage (important for laptop)
model.set_chunk_size(64)

os.makedirs(OUTDIR, exist_ok=True)

for record in SeqIO.parse(FASTA_FILE, "fasta"):
    seq = str(record.seq)
    name = record.id

    print(f"Processing {name}...")

    with torch.no_grad():
        pdb = model.infer_pdb(seq)

    with open(f"{OUTDIR}/{name}.pdb", "w") as f:
        f.write(pdb)

print("✅ Done")