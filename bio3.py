from collections import defaultdict, Counter
import pprint

alignment_A = [
    "-CCTTTAGAACGCATTTAGGTCT-A-AAATCAGTCT",
    "-GGATGAGAACGCATTT-GGACTCAGATTTCAGTAA",
    "---TATAG-ACGCA-TTAGGACTCA-AATTTAGT-T",
    "TGCATTAGAACGCATTTA-GACTCA-ATTTCCGTCT",
    "-AGTTT-GAACGCATTTAGGACTCA-AATTTCGT-A",
    "--TATTAGAACGCATTTAGGACTCA-ATTTCAGT-A",
    "---ATTA-AAGGCATTTAGGACTCAAAGTTCAGTGA",
    "-AAATTAGAACGC-CTTAGGACTAAAATTTCAGT-T",
    "--AATTAG-ACGCATTTAGGACTAAAATTTCAGTTC",
    "--GGATAGAACGCATTTAGGACTCA-AATTCAGTTA"
]

# Dataset B 
alignment_B = [
"TATTAGAACGCATTGGAGTCAAAATCAGTC",
"CATTAGAACGAGTTTGGGACTAAAGTTCAGTTT",
"AAATTAGAACGCATAGGACTCAATTTCAGTCG",
"GTATAGAACGCATTAGGACTCAAATTCAGCAT",
"CCTTTAAACGCATTAGGACTCAATTTCATG",
"CAATTAGAACGCATTAGGACTCAAATTTCAGTCC",
"TTATTGACCGCATTTAGGACTCAATTCATGG",
"TATCAGAACGCATTAGGACTCGATTTTCAGTCA",
"AATTAGAACGCATTTAGACTCAATTTCAGTA",
"CTAGAAGCATTGGACTCAAATTCTGTT",
"CGATAGAACGCATTTAGGACTCAAGTCTCAGTA",
"CCTATGAACGCATTCAGGACTCAAATTTCAGTTA",
"CATGAACGCATTTGGGACTCACATTCAGTG",
"GATATTGAACCATTTCGGTCTCAAACTTCGTT",
"CATTGAACGCATGAGGACTCAATTTAGTA",
"AATTAGAACGCATCAGGACTCAAATTCAGTT",
"GTAATTGAAGCATTAAGCTCAAATTTAGTT",
"TATTTGGACGCATTTAGGACACAATTCAGTTC",
"CACATTAGAACGCATTTAGGATCAAATCTCAGGGA",
"ACATTAGAACGTATTTAGGACTAAATTTCAGTG",
"CCTAATAGAACGCATTTAGACTCAAATTTCAGTGG",
"AAGTAGAACGCATTTTGGCTCAATTTCAGTG",
"CATTAGAACGCATTTAGGACTCAAATTCCATTC",
"GCGATTAGAACGCATTTAGGACTAAATTCGTG",
"TATTAGAACGTATTAGGACTCAAATTTCAGTGA",
"GATTGCAACGATTTAGGACCAAATTTAGTCG",
"TAGAGATCGCATTAGGATCAAATTCAGTCC",
"TTGTTGACGGATTTAGGACTCAAATTTCAGTGT",
"GTATTAAGACATTTAGGACTCAAATTTCAGTC",
"CATTGAACGCATTTAGGACTCAATTCAGTTA",
"CAATTAGACGCATTAAGGACTCAAATTTCATGC",
"CTATTTGCGCATTTAGGACCAGATTTCAGTC",
"AATTAACCGCATTTAGACTCAAATTTCATTC",
"TAATTAGACCGGATTTAGGACCAAATTTCCGTAA",
"AATTAGAACGCATTTAGGACTCAAATGTCAGTTT",
"CATTAGAACGCATTAGGACTCCAATTTCAGTCG",
"CACATTCGAACATTTAGGACTCAAAATTCAGT",
"CCGATTGGAACGCATTTAGGACACAATTTCAGTC",
"CGATTAGAACGCATTTAGTACTGAAATTTCGTG",
"CCATTAGACGCATTTGGGACCAAATTCAGTT",
"GTATTAGAACACATTAGGACTCAAATTTCATAG",
"AATTAACGCATTAAGACTCAAATTTCAGTC",
"CAAATAGAACGCATTAGGACTTAAATTTCGTCT",
"TATTAGAACCCTTTAGGACTCAATTTTCAGTT",
"CATTAGACGCATTGAGGACTCAAATTTCAGTGA",
"GAAACAGAACGCATTTAGGACTCAAATTTCAGTA",
"TGTTAGAACGCAATTAGGACTCAAATTTCAGTA",
"ATATTGAAGGCATTTAGGCTCAAATTTCAGCG",
"GAATTAGAACGCATTTAGGAACAATTTCGTA",
"AATAGAAGCATTTAGGACTAAATTTCGTT",
"CATTAGAACGCATTTAGGACTCAAATTTAAGTAA",
"ATTCGAACACATTTAGGACTCAATTGCAGTGG",
"AATTAGAAGCATTGAGGACTCAAATTCAGTGT",
"TTCATAGACGCATTTAGGACTCAAATTCCAGTC",
"GGTAGAACGCATTTATGACTCAATTTCAGCA",
"CATTAACGCATTTAGGAATCATATTCAGTT",
"GTATTAGAACGCATTTTGGACTCACACTTCAGTTA",
"TAAAGTAGAACGCATGTAGCACTCAAATTACAGTG",
"GGTATTGACGCATGTAAGACTCCAATTTCAGTT",
"CATTATAACGCCTTGGACTCAAATTTCAGTT",
"GATTATACGCATTTAGGACTCAAATTTCGGTA",
"GATTAGACGCATTAGGTCAAATTTCAGTC",
"GGATAAACGCATTAGACTCAAATTTCAGTC",
"AGATTAGAACGCCTCTAGGACTCAATTCAGTT",
"GATAGGCGTATTTAGGACTCAAATTTCAGTT",
"CAGTAGAACGCATTTAGGACTCAAATTTCAGTA",
"GCTTAGAACGCATTTAGACTCCAGTTCAGTTG",
"TATGGAAGCATTAGGACTCAATTTCAGTGT",
"ACCTTTGAACGATTTAGGACTCAATTTCAGTTT",
"TATCTGAACGCATTTACAACTCAAATTTCAGACT"
]

num_seqs_A = len(alignment_A)
align_len = len(alignment_A[0])
match_threshold = 0.6

#1: identify match positions (from A)
match_positions = []
for col in range(align_len):
    non_gaps = sum(seq[col] != '-' for seq in alignment_A)
    if non_gaps / num_seqs_A >= match_threshold:
        match_positions.append(col)

#2: label sequences (from datasetA structure)
def label_sequence(seq, match_positions):
    labeled = []
    m_idx = 1
    for i in range(len(seq)):
        if i >= len(match_positions):
            break
        pos = match_positions[i]
        if pos < len(seq):
            char = seq[pos]
            if char != '-':
                labeled.append((f"M{m_idx}", char))
            else:
                labeled.append((f"D{m_idx}", None))
            m_idx += 1
    return labeled

#3: for datasetA — build HMM from alignment
labeled_sequences = []
for seq in alignment_A:
    labeled = []
    m_idx = 1
    for i, char in enumerate(seq):
        if i in match_positions:
            state = f"M{m_idx}" if char != '-' else f"D{m_idx}"
            if char != '-':
                labeled.append((state, char))
            else:
                labeled.append((state, None))
            m_idx += 1
        else:
            if char != '-':
                state = f"I{m_idx-1}"
                labeled.append((state, char))
    labeled_sequences.append(labeled)

transitions = defaultdict(Counter)
emissions = defaultdict(Counter)

for labeled_seq in labeled_sequences:
    prev_state = 'S'
    for state, symbol in labeled_seq:
        transitions[prev_state][state] += 1
        if symbol:
            emissions[state][symbol] += 1
        prev_state = state
    transitions[prev_state]['E'] += 1

#4: train with B
for seq in alignment_B:
    labeled = label_sequence(seq, match_positions)
    prev_state = 'S'
    for state, symbol in labeled:
        transitions[prev_state][state] += 1
        if symbol:
            emissions[state][symbol] += 1
        prev_state = state
    transitions[prev_state]['E'] += 1

#5: normalize

def normalize(counter):
    total = sum(counter.values())
    return {k: round(v / total, 3) for k, v in counter.items()}

transitions_normalized = {k: normalize(v) for k, v in transitions.items()}
emissions_normalized = {k: normalize(v) for k, v in emissions.items()}

#print results
print("=== Match Positions (indices) ===")
print(match_positions)
print("\n=== Transitions (Probabilities) ===")
pprint.pprint(transitions_normalized)
print("\n=== Emissions (Probabilities) ===")
pprint.pprint(emissions_normalized)
