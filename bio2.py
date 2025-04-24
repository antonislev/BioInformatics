def needleman_wunsch(seq1, seq2, alpha=1):
    match_score = 1
    mismatch_penalty = -alpha / 2
    gap_penalty = -alpha

    m, n = len(seq1), len(seq2)
    dp = [[0] * (n+1) for _ in range(m+1)]
    traceback = [[None] * (n+1) for _ in range(m+1)]

    for i in range(1, m+1):
        dp[i][0] = dp[i-1][0] + gap_penalty
        traceback[i][0] = 'up'
    for j in range(1, n+1):
        dp[0][j] = dp[0][j-1] + gap_penalty
        traceback[0][j] = 'left'

    for i in range(1, m+1):
        for j in range(1, n+1):
            score = match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty
            diag = dp[i-1][j-1] + score
            up = dp[i-1][j] + gap_penalty
            left = dp[i][j-1] + gap_penalty
            dp[i][j] = max(diag, up, left)

            if dp[i][j] == diag:
                traceback[i][j] = 'diag'
            elif dp[i][j] == up:
                traceback[i][j] = 'up'
            else:
                traceback[i][j] = 'left'

    aligned1, aligned2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        move = traceback[i][j]
        if move == 'diag':
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif move == 'up':
            aligned1 = seq1[i-1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
        elif move == 'left':
            aligned1 = '-' + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2

def multiple_alignment(sequences, alpha=1):
    aligned_seqs = [sequences[0]]

    for seq in sequences[1:]:
        # Στοίχιση του κάθε επόμενου στη μέχρι τώρα στοίχιση
        new_aligned = []
        ref = aligned_seqs[0]
        ref, new_seq = needleman_wunsch(ref, seq, alpha)

        # Ενημέρωση και των υπολοίπων παλιών σειρών ώστε να ταιριάζουν με ref
        for old_seq in aligned_seqs:
            aligned_old, aligned_ref = needleman_wunsch(old_seq, ref, alpha)
            new_aligned.append(aligned_old)

        new_aligned.append(new_seq)
        aligned_seqs = new_aligned

    return aligned_seqs
datasetA = [
    "CCTTTAGAACGCATTTAGGTCTAAAATCAGTCT",
    "GGATGAGAACGCATTTGGACTCAGATTTCAGTAA",
    "TATAGACGCATTAGGACTCAAATTTAGTT",
    "TGCATTAGAACGCATTTAGACTCAATTTCCGTCT",
    "AGTTTGAACGCATTTAGGACTCAAATTTCGTA",
    "TATTAGAACGCATTTAGGACTCAATTTCAGTA",
    "ATTAAAGGCATTTAGGACTCAAAGTTCAGTGA",
    "AAATTAGAACGCCTTAGGACTAAAATTTCAGTT",
    "AATTAGACGCATTTAGGACTAAAATTTCAGTTC",
    "GGATAGAACGCATTTAGGACTCAAATTCAGTTA"
]

aligned = multiple_alignment(datasetA, alpha=1)

# Εκτύπωση ευθυγραμμισμένων σειρών
for i, seq in enumerate(aligned):
    print(f"Seq {i+1}: {seq}")

