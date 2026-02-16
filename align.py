import sys

# ----------------------------
# FASTA Reader
# ----------------------------
def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence.upper()


# ----------------------------
# Parameter Reader
# ----------------------------
def read_parameters(filename):
    params = {}
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=")
            params[key.strip()] = float(value.strip())
    return params



# ----------------------------
# Smith-Waterman (Local Alignment)
# with Affine Gap Penalties
# ----------------------------
def local_alignment_affine(seq1, seq2, match, mismatch, gap_open, gap_extend):
    n = len(seq1)
    m = len(seq2)

    # DP matrices
    M = [[0]*(m+1) for _ in range(n+1)]
    Ix = [[0]*(m+1) for _ in range(n+1)]
    Iy = [[0]*(m+1) for _ in range(n+1)]

    # Traceback matrices
    trace_M = [[None]*(m+1) for _ in range(n+1)]
    trace_Ix = [[None]*(m+1) for _ in range(n+1)]
    trace_Iy = [[None]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_position = (0, 0, "M")

    # Fill matrices
    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            # ---- M matrix ----
            candidates = [
                (0, None),
                (M[i-1][j-1] + score, "M"),
                (Ix[i-1][j-1] + score, "Ix"),
                (Iy[i-1][j-1] + score, "Iy")
            ]
            best = max(candidates, key=lambda x: x[0])
            M[i][j] = best[0]
            trace_M[i][j] = best[1]

            # ---- Ix matrix (gap in seq2) ----
            candidates = [
                (0, None),
                (M[i-1][j] + gap_open, "M"),
                (Ix[i-1][j] + gap_extend, "Ix")
            ]
            best = max(candidates, key=lambda x: x[0])
            Ix[i][j] = best[0]
            trace_Ix[i][j] = best[1]

            # ---- Iy matrix (gap in seq1) ----
            candidates = [
                (0, None),
                (M[i][j-1] + gap_open, "M"),
                (Iy[i][j-1] + gap_extend, "Iy")
            ]
            best = max(candidates, key=lambda x: x[0])
            Iy[i][j] = best[0]
            trace_Iy[i][j] = best[1]

            # Track maximum score
            for matrix_name, value in [("M", M[i][j]), ("Ix", Ix[i][j]), ("Iy", Iy[i][j])]:
                if value > max_score:
                    max_score = value
                    max_position = (i, j, matrix_name)

    # ----------------------------
    # Traceback
    # ----------------------------
    aligned1 = ""
    aligned2 = ""

    i, j, matrix = max_position

    while i > 0 and j > 0:
        if matrix == "M":
            if M[i][j] == 0:
                break
            prev = trace_M[i][j]
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
            matrix = prev

        elif matrix == "Ix":
            if Ix[i][j] == 0:
                break
            prev = trace_Ix[i][j]
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
            matrix = prev

        elif matrix == "Iy":
            if Iy[i][j] == 0:
                break
            prev = trace_Iy[i][j]
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1
            matrix = prev

        else:
            break

    return aligned1, aligned2, max_score


# ----------------------------
# Main Program
# ----------------------------
if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python align.py seq1.fasta seq2.fasta parameters.txt")
        sys.exit(1)

    seq1_file = sys.argv[1]
    seq2_file = sys.argv[2]
    param_file = sys.argv[3]

    # Read inputs
    seq1 = read_fasta(seq1_file)
    seq2 = read_fasta(seq2_file)
    params = read_parameters(param_file)

    match = params["match"]
    mismatch = params["mismatch"]
    gap_open = params["gap_open"]
    gap_extend = params["gap_extend"]

    # Run alignment
    aligned1, aligned2, score = local_alignment_affine(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    # Write output
    with open("output.txt", "w") as f:
        f.write(f"Best Alignment Score: {score}\n")
        f.write(aligned1 + "\n")
        f.write(aligned2 + "\n")

    print("Alignment written to output.txt")

