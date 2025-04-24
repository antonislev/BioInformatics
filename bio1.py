import random

# Βασικό αλφάβητο
alphabet = ['A', 'C', 'G', 'T']

# Τα 4 patterns
patterns = ['ATTAGA', 'ACGCATTT', 'AGGACTCAA', 'ATTTCAGT']

def mutate_pattern(pattern):
    pattern = list(pattern)
    num_mutations = random.randint(0, 2)
    positions = random.sample(range(len(pattern)), num_mutations)
    for pos in positions:
        if random.random() < 0.5:  # 50% να γίνει διαγραφή
            pattern[pos] = ''
        else:
            pattern[pos] = random.choice(alphabet)
    return ''.join(pattern)

def synthesize_sequence():
    # 1. Τυχαία αρχικά σύμβολα (1 έως 3)
    result = ''.join(random.choices(alphabet, k=random.randint(1, 3)))

    # 2. Για κάθε pattern, κάνουμε μετασχηματισμό και το προσθέτουμε
    for p in patterns:
        result += mutate_pattern(p)

    # 3. Τυχαία σύμβολα στο τέλος (1 ή 2)
    result += ''.join(random.choices(alphabet, k=random.randint(1, 2)))

    return result

# Δημιουργία 100 συμβολοσειρών
all_sequences = [synthesize_sequence() for _ in range(100)]

# Ανακάτεμα για τυχαία επιλογή
random.shuffle(all_sequences)

# Κατανομή στα σύνολα
datasetA = all_sequences[:10]
datasetB = all_sequences[10:80]
datasetC = all_sequences[80:]

# Αποθήκευση των datasets σε αρχεία
with open("datasetA.txt", "w") as f:
    f.write("\n".join(datasetA))

with open("datasetB.txt", "w") as f:
    f.write("\n".join(datasetB))

with open("datasetC.txt", "w") as f:
    f.write("\n".join(datasetC))

print("Τα δεδομένα αποθηκεύτηκαν σε αρχεία datasetA.txt, datasetB.txt, datasetC.txt")
