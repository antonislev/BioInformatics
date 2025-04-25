import numpy as np

# Ορισμός του datasetC
datasetC = [
    "GGTTAGAACGCATTAGGACTCAAATTTCAGTT",
    "CGCATTAGAACGCATTTAGGATACAAATTTCTGTTA",
    "TATTAGAACGCATTTAGGACTCAAATTTCAGTAT",
    "TCTTAGGACGCATAAGGACTCAAATTTGAGTA",
    "AGTAATAGAACGCATTTAGGACTCAAGTTTCATAC",
    "TCCATTAGAACGCATTTAGGACTAAATTTCAGTC",
    "CCATTAGAAGCATTTAGCCTCAAACTTCATA",
    "CTATAGCACGCATTTAGGACTCGAATTTCAGTGT",
    "GAGACTAGAACGCTCTAGGACTCAAATTCAGG",
    "GAATAAGAACGCATTTTGGACTCAAATTTCAGT",
    "CCTATTAGAACGAATTACGACTCAATCTCAGTT",
    "TGATTAGAACGCATTTAGGACTCAAATTTCAGTT",
    "AGATAGTACGCATTTAGGTCTGAAATTTCCGTG",
    "ATCTTAGAACGCATTTAGGACTCAACTTTCAGTCT",
    "CCATTAGAACGCATTTGGACGCAAATTTCAGTC",
    "CTTAGAACGCATTTAGACTCAAATTTAGTGA",
    "AGATAGAACGCATTTAGGACTCAATTCAGTA",
    "GTACTTGAACGCATTTAGGACTCAATTCAGTAA",
    "TATTAAACGCTTAGGAGTCACATTCATA",
    "TTATTAACGCATTTAGACTCAAATTCCAGAGC"
]
datasetCA = [
    "AGCTGATCGAGGTCGATCGTAAAGCTGAGTTCG",
    "GCTAGCTAGGCTAGCAGTCCATGATAGCGTCTA",
    "ATCGTAGCGTACGGTACATGCTAGGTCGATACG",
    "GTGCTAGCATGCATGCTAGAGGTCATGTAGTCA",
    "TGCATAGCTAGAGTGCATAGCGTAGGACATCGAT",
    "AGCTAGCTAGGCTAGCTAGGTCATAGCAGGTA",
    "CTAGATCGGATCGTACGATCGGATCTAGGTAGC",
    "GATCGTAGCGAGTCAGGTACGATGCTAAGTTC",
    "GTAGCGTACGAGCTAGCGATAGTCAGGACTGCT",
    "TACGTAGCTGAGTCAGGAGCATTAGCTACGATC",
    "CTAGGTCAGGCTAGTCGAGGATTCAGGATGACT",
    "CGTATCGATCGTACGGTAGAGGACGTAAGTACG",
    "GACGTAGCTAGATCGTATCAGTAGTGCATGCA",
    "AGCATTAGTCAGGACAGCTAGGTCATCGTACGA",
    "GCTAGTACGTAGAGTCATCGGATCGTAGGACT",
    "TACGATGCATAGCAGGTAGGCTATGACGTAGC",
    "GTACGTAGCTAGCTAGGACTCGTAGACTGACAT",
    "TAGGCTAGACGATCGTAGTCAGTAGCAGTGTG",
    "CTAGGATCGTACGAGTCATAGCGTACAGTACGT"
]


# scores
match_score = 1
mismatch_score = -1
gap_penalty = -2

# calc alignment for seq pair
def calculate_alignment(seq1, seq2):
    # DP
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    dp = np.zeros((len_seq1 + 1, len_seq2 + 1))

    # gap penalty
    for i in range(1, len_seq1 + 1):
        dp[i][0] = dp[i-1][0] + gap_penalty
    for j in range(1, len_seq2 + 1):
        dp[0][j] = dp[0][j-1] + gap_penalty

    # scores DP
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match, delete, insert)

    # remake of alignment path
    alignment1 = []
    alignment2 = []
    i, j = len_seq1, len_seq2
    while i > 0 or j > 0:
        current_score = dp[i][j]
        if i > 0 and dp[i-1][j] + gap_penalty == current_score:
            alignment1.append(seq1[i-1])
            alignment2.append("-")
            i -= 1
        elif j > 0 and dp[i][j-1] + gap_penalty == current_score:
            alignment1.append("-")
            alignment2.append(seq2[j-1])
            j -= 1
        else:
            alignment1.append(seq1[i-1])
            alignment2.append(seq2[j-1])
            i -= 1
            j -= 1

    # reverse results for correct order
    alignment1.reverse()
    alignment2.reverse()

    return ''.join(alignment1), ''.join(alignment2), dp[len_seq1][len_seq2]

# calc alignment 
for i in range(len(datasetC)):
    for j in range(i+1, len(datasetC)):  # check for evry pair
        print(f"Aligning Sequence {i+1} and Sequence {j+1}:")
        aligned_seq1, aligned_seq2, score = calculate_alignment(datasetC[i], datasetC[j])
        print("Alignment Score:", score)
        print("Aligned Sequences:")
        print(aligned_seq1)
        print(aligned_seq2)
        print("\n" + "-"*50 + "\n")



