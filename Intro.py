import matplotlib.pyplot as plt

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if not t[i + j] == p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mis = 0
        for j in range(len(p)):
            if not t[i + j] == p[j]:
                mis += 1
        if mis <= 2:
            occurrences.append(i)
    return occurrences


def naive_with_rc(p, t):
    r = reverse(p)
    if r == p:
        return naive(p, t)
    occurrences = naive(p, t)
    for i in range(len(t) - len(r) + 1):
        match = True
        for j in range(len(r)):
            if not t[i + j] == r[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences


def reverse(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


lambda_virus = readGenome('lambda_virus.fa')
occur = naive_with_rc('AGTCGA', lambda_virus)
# print('len of lambda_virus: ', len(lambda_virus))
# print('leftmost occurrences: %d' % min(occur))
# print('# occurrences: %d' % len(occur))

occur2 = naive_2mm('AGGAGGTT', lambda_virus)
# print('# occurrences: %d' % len(occur2))
# print('leftmost occurrences: %d' % min(occur))

def findGCByPos(reads):
    gc = [0] * 100
    totals = [0] * 100

    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])

    return gc


def readFastq(filename):
    sequences = []
    qualites = []
    with open(filename, 'r') as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break;
            sequences.append(seq)
            qualites.append(qual)
    return sequences, qualites
seqs, quals = readFastq('ERR037900_1.first1000.fastq')

gc = findGCByPos(seqs)
# plt.plot(range(len(gc)), gc)
# plt.show()
print(gc[60:70])