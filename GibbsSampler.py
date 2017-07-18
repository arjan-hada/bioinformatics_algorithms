# -*- coding: utf-8 -*-
import random
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

def scoreMotifs(motifs):
    '''This function computes the score of list of motifs'''
    z = zip(*motifs)
    totalscore = 0
    for string in z:
        score = len(string)-max([string.count('A'),string.count('C'), string.count('G'), string.count('T')])
        totalscore += score
    return totalscore   

def random_kmers(dna, k):
    '''Returns a collection of randomly chosen k-mers in Dna. t is the number of items in dna list'''
    motifs = []
    for text in dna:
        i = random.randrange(0,len(dna[0])-k+1)
        motifs.append(text[i:i+k])
    return motifs      
                  
def Random(p):
    wheel = [0]*(len(p)+1)
    #Normalization of probability distribution
    s = sum(p)
    if s!= float(1):
        p = [float(i)/sum(p) for i in p]
    #Make a roulette wheel
    for index in range(len(p)):
        wheel[index+1] = wheel[index] + p[index]
    #print wheel
    #Rolling a biased dice
    r = random.random() #Randomly generate a number between 0 and 1
    #print r
    for i in range(len(wheel)-1):
        if wheel[i] < r < wheel[i+1]:
            result = i
    return result
    
def profileRandomlyGeneratedkmer(text, k , matrix):
    '''Input: A string Text, an integer k, and a 4 × k matrix Profile.
     Output: A Profile-most probable k-mer in Text.'''
    prob = []
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        pt = 1
        for j in range(k):
            p = matrix[kmer[j]][j]
            pt *=p
        prob.append(pt)
    #print prob
    kmer_index = Random(prob)
    randomkmer = text[kmer_index:kmer_index+k]
    return randomkmer
    
def gibbsSampler(dna, k, t, N):
    motifs = random_kmers(dna, k)
    bestMotifs = motifs
    for j in range(N):
        i = random.randrange(0,t)
        del motifs[i]
        profile = motifsToProfileLaplace(motifs)
        motifi = profileRandomlyGeneratedkmer(dna[i], k , profile)
        motifs.insert(i, motifi)
        if scoreMotifs(motifs) < scoreMotifs(bestMotifs):
            bestMotifs = motifs
    return bestMotifs
    
def loopGibbsSampler(dna, k, t, N):
    bestMotifs = None
    bestScore = None
    n = 0
    while n <= 20:
        n += 1
        motifs = gibbsSampler(dna, k, t, N)
        score = scoreMotifs(motifs)
        if bestScore is None or score < bestScore:
            bestScore = score
            bestMotifs = motifs
    return bestMotifs, bestScore
    
dna = [line.strip() for line in open('dataset_163_4.txt', 'r')]
bestMotifs,bestScore = loopGibbsSampler(dna, 15, 20, 2000)
x = ' '.join(str(i) for i in bestMotifs)