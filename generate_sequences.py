import random

# Fixed seed for reproducibility
random.seed(42)

def random_dna(length):
    return ''.join(random.choice("ATGC") for _ in range(length))

def mutate_sequence(seq):
    seq = list(seq)

    # --- Substitutions (2–7) ---
    num_sub = random.randint(2, 7)
    for _ in range(num_sub):
        pos = random.randint(0, len(seq)-1)
        original = seq[pos]
        choices = [b for b in "ATGC" if b != original]
        seq[pos] = random.choice(choices)

    # --- Deletions (1–4) ---
    num_del = random.randint(1, 4)
    for _ in range(num_del):
        if len(seq) > 1:
            pos = random.randint(0, len(seq)-1)
            seq.pop(pos)

    # --- Insertions (1–3) ---
    num_ins = random.randint(1, 3)
    for _ in range(num_ins):
        pos = random.randint(0, len(seq))
        seq.insert(pos, random.choice("ATGC"))

    return ''.join(seq)

# Generate base sequence
seq1 = random_dna(20)

# Generate modified sequence
seq2 = mutate_sequence(seq1)

# Write FASTA files
with open("seq1.fasta", "w") as f:
    f.write(">Sequence_1\n")
    f.write(seq1 + "\n")

with open("seq2.fasta", "w") as f:
    f.write(">Sequence_2\n")
    f.write(seq2 + "\n")

print("Sequence 1:", seq1)
print("Sequence 2:", seq2)
print("Files seq1.fasta and seq2.fasta created.")
