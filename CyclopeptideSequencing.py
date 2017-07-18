aminoacid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
aminoacidMass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103, 'L':113, 'N':114, 'D':115, 'K':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}
def linearSpectrum(peptide):
    """Input: An amino acid string Peptide.
     Output: The linear spectrum of Peptide."""
    prefixMass = [0]*((len(peptide)+1))
    for i in range(len(peptide)):
        prefixMass[i+1] = prefixMass[i] + aminoacidMass[peptide[i]]
    #print 'prefixMass', prefixMass
    linear_spectrum = [0]
    for i in range(len(prefixMass)-1):
        for j in range(i+1, len(prefixMass)):
            linear_spectrum.append(prefixMass[j] - prefixMass[i])
    return sorted(linear_spectrum)    

def cyclicSpectrum(peptide):
    """Input: An amino acid string Peptide.
     Output: The cyclic spectrum of Peptide."""
    prefixMass = [0]*((len(peptide)+1))
    for i in range(len(peptide)):
        prefixMass[i+1] = prefixMass[i] + aminoacidMass[peptide[i]]
    peptideMass = prefixMass[len(peptide)]
    cyclic_spectrum = [0]
    for i in range(len(prefixMass)-1):
        for j in range(i+1, len(prefixMass)):
            cyclic_spectrum.append(prefixMass[j] - prefixMass[i])
            if i > 0 and j < (len(prefixMass)-1):
                cyclic_spectrum.append(peptideMass - (prefixMass[j] - prefixMass[i]))
    return sorted(cyclic_spectrum) 

def expand(peptides):
    expanded = []
    for i in peptides:
        for j in aminoacid:
            expanded.append(i+j)
    return expanded

def mass(peptide):
    """Calculates the mass of peptide using the aminoacidMass dictionary"""
    massOfPeptide = 0
    for i in peptide:
        massOfPeptide += aminoacidMass[i]
    return massOfPeptide  

from collections import Counter
def isSubset(list1, list2):
    c1, c2 = Counter(list1), Counter(list2)
    for k, n in c1.items():
        if n > c2[k]:
            return False
    return True
    
def cycloPeptideSequencing(spectrum):
    peptides = [k for k,v in aminoacidMass.items() if v in spectrum]
    while len(peptides) != 0:
        peptides = expand(peptides)
        copy = peptides
        for peptide in copy:
            if mass(peptide) == max(spectrum):
                if cyclicSpectrum(peptide) == spectrum:
                    return peptide
                peptides.remove(peptide)
            else:
                if isSubset(linearSpectrum(peptide), spectrum) == False:
                    peptides.remove(peptide)    

#Read the data file and convert to list of numbers as integers
f = open('input/cycloseq_data.txt', 'r') 
for line in f: #line is a string
    numbers = line.split() # split the string on white-space and return a list of numbers as strings
    spectrum = map(int, numbers) #convert numbers to integers

#spectrum = [0, 113, 128, 186, 241, 299, 314, 427]   
#Calling the function
ans = cycloPeptideSequencing(spectrum)