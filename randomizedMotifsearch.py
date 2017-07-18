# -*- coding: utf-8 -*-
def profileMostProbablekmer(text, k , matrix):
    '''Input: A string Text, an integer k, and a 4 × k matrix Profile.
     Output: A Profile-most probable k-mer in Text.'''
    maxp = None
    probablekmer = None
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        pt = 1
        for j in range(k):
            p = matrix[kmer[j]][j]
            pt *=p
        if maxp == None or pt > maxp:
            maxp = pt
            probablekmer = kmer
    return probablekmer
    
def motifsToProfile(motifs):
    '''Given a list of motifs, the functions converts it into dictionary containing frequency profile of each nucleotide.'''
    z = zip(*motifs)
    d = {}
    n = float(len(motifs))
    for i in range(len(z)):
        d.setdefault('A', []).append(z[i].count('A')/n)
        d.setdefault('C', []).append(z[i].count('C')/n)
        d.setdefault('G', []).append(z[i].count('G')/n)
        d.setdefault('T', []).append(z[i].count('T')/n)
    return d
    
def scoreMotifs(motifs):
    '''This function computes the score of list of motifs'''
    z = zip(*motifs)
    totalscore = 0
    for string in z:
        score = len(string)-max([string.count('A'),string.count('C'), string.count('G'), string.count('T')])
        totalscore += score
    return totalscore 
    
def motifsToProfileLaplace(motifs):
    '''Given a list of motifs, the functions converts it into dictionary containing frequency profile of each nucleotide.'''
    z = zip(*motifs)
    d = {}
    n = float(len(motifs))
    for i in range(len(z)):
        d.setdefault('A', []).append((z[i].count('A')+1)/(2*n))
        d.setdefault('C', []).append((z[i].count('C')+1)/(2*n))
        d.setdefault('G', []).append((z[i].count('G')+1)/(2*n))
        d.setdefault('T', []).append((z[i].count('T')+1)/(2*n))
    return d
    
def motifsFromProfile(dna, k , profile):
    '''Computes the profileMostProbablekmer for each text in dna and returns it as a list.'''
    return [profileMostProbablekmer(text, k , profile) for text in dna]
    
def random_kmers(dna, k):
    '''Returns a collection of randomly chosen k-mers in Dna. t is the number of items in dna list'''
    import random
    motifs = []
    for text in dna:
        i = random.randrange(0,len(dna[0])-k+1)
        motifs.append(text[i:i+k])
    return motifs
    
def randomizedMotifSearch(dna, k):
    motifs = random_kmers(dna, k)
    bestMotifs = motifs
    while True:
        profile = motifsToProfileLaplace(motifs)
        motifs = motifsFromProfile(dna, k , profile)
        if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
            bestMotifs = motifs
        else:
            return bestMotifs

def looprandomizedMotifSearch(dna, k, nloop):
    bestMotifs = None
    bestScore = None
    n = 0
    while n <= nloop:
        n += 1
        motifs = randomizedMotifSearch(dna, k)
        score = scoreMotifs(motifs)
        if bestScore is None or score < bestScore:
            bestScore = score
            bestMotifs = motifs
    return bestMotifs, bestScore


dna = [line.strip() for line in open('dataset_161_5.txt', 'r')]
bestMotifs,bestScore = looprandomizedMotifSearch(dna, 15, 2000)
x = ' '.join(str(i) for i in bestMotifs)
