import random
import math

# Variable declarations
k = 10
max_iterations = 1000

# Methods
def print_mat(mat):
    for x in mat:
        print(x)

def read_file(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines()] # To get rid of the trailing \n
    return lines

def RandomlySelectKmers(dna, k):
    t = len(dna)
    kmer_list = []
    for i in range(t):
        start = random.randint(0, len(dna[i]) - k)
        kmer_list.append(dna[i][start:start + k])
    return kmer_list

def RandomlySelectElimMotif(dna):
    t = len(dna)
    return random.randint(0, t-1)


def Profile(motifs):
    t = len(motifs)
    k = len(motifs[0])
    profile = [[1 for i in range(k)] for j in range(4)] # Laplace's Rule of Succession (Init at 1)
    for i in range(t):
        for j in range(k):
            if (motifs[i][j] == 'A'):
                profile[0][j] += 1
            elif (motifs[i][j] == 'C'):
                profile[1][j] += 1
            elif (motifs[i][j] == 'G'):
                profile[2][j] += 1
            elif (motifs[i][j] == 'T'):
                profile[3][j] += 1
    for i in range(4):
        for j in range(k):
            profile[i][j] /= (t + 4)
    return profile

def MOTIFS(text, k, profile):
    max_prob = -1
    kmer = text[:k]
    for i in range(len(text) - k + 1):
        prob = 1
        for j in range(k):
            if (text[i + j] == 'A'):
                prob *= profile[0][j]
            elif (text[i + j] == 'C'):
                prob *= profile[1][j]
            elif (text[i + j] == 'G'):
                prob *= profile[2][j]
            elif (text[i + j] == 'T'):
                prob *= profile[3][j]
        if (prob > max_prob):
            max_prob = prob
            kmer = text[i:i + k]
    return kmer

# def Score(motifs):
#     t = len(motifs)
#     k = len(motifs[0])
#     score = 0
#     for i in range(k):
#         count = [0, 0, 0, 0]
#         for j in range(t):
#             if (motifs[j][i] == 'A'):
#                 count[0] += 1
#             elif (motifs[j][i] == 'C'):
#                 count[1] += 1
#             elif (motifs[j][i] == 'G'):
#                 count[2] += 1
#             elif (motifs[j][i] == 'T'):
#                 count[3] += 1
#         score += (t - max(count))
#     return score

def Score(motifs):
    t = len(motifs)
    k = len(motifs[0])
    score = 0
    for j in range(k):
        count = {'A': 1, 'C': 1, 'G': 1, 'T': 1}  
        # Laplace's Rule of Succession (Init at 1)
        for i in range(t):
            count[motifs[i][j]] += 1
        entropy = 0
        for nucleotide, nucleotide_count in count.items():
            probability = nucleotide_count / (t + 4)  
            # Pseudocount for Laplace's Rule of Succession
            entropy -= probability * math.log2(probability)
        score += 2 - entropy  # Maximized entropy is 2 for 4 nucleotides
    return score


def entropy(motifs):
    t = len(motifs)
    k = len(motifs[0])
    entropy_score = 0
    for j in range(k):
        count = {'A': 1, 'C': 1, 'G': 1, 'T': 1}  # Laplace's Rule of Succession with pseudocounts
        for i in range(t):
            count[motifs[i][j]] += 1
        entropy = 0
        for nucleotide, nucleotide_count in count.items():
            probability = nucleotide_count / (t + 4)  # Pseudocount for Laplace's Rule of Succession
            if probability > 0:  # Avoid log(0)
                entropy -= probability * math.log2(probability)
        entropy_score += entropy
    return entropy_score / k  # Normalization by dividing by the total number of positions


# Algorithm for Gibbs Sampling
def GibbsSampler(dna,k,N):
    
    motifs = RandomlySelectKmers(dna, k)
    best_motifs = motifs

    for i in range(N):
        motifs = RandomlySelectKmers(dna, k)
        motif_num = RandomlySelectElimMotif(dna)
        elemMotif = motifs[motif_num]
        profileMotif = motifs[:]
        del profileMotif[motif_num]
        
        profile = Profile(profileMotif)
        
        motifs[motif_num] = MOTIFS(dna[motif_num], k, profile) # P - most probable k-mer in the i-th string
        if (Score(motifs) < Score(best_motifs)):
            best_motifs = motifs

     # Generate report
    num_sequences = len(dna)
    print("Number of sequences available:", num_sequences)
    print("-" * 47)
    print("{:<12} {:<15} {:<20}".format("Sequence", "Motif", "Sequence length"))
    print("-" * 47)
    for i, (sequence, motif) in enumerate(zip(dna, best_motifs)):
        print("{:<12} {:<15} {:<20}".format(f"Sequence {i+1}", motif, len(sequence)))
    print("-" * 47)
    print("Entropy of the generated motifs:", entropy(best_motifs))

    return best_motifs
    
dna = read_file("hm03.txt")
best_motifs = GibbsSampler(dna,k,max_iterations)       
# print_mat(best_motifs)
# print(Score(best_motifs))
