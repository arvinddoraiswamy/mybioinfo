# Input:  Two strings, Pattern and Genome
# Output: A list containing all starting positions where Pattern appears as a substring of Genome

def find_oris(Pattern, Genome):
    len_pattern= len(Pattern)
    len_Genome=  len(Genome)
    diff= len_Genome - len_pattern + 1 

    all_kmers= [Genome[i:i+len_pattern] for i in range(0, diff)]
    return all_kmers

def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    # your code here
    all_kmers= find_oris(Pattern, Genome)
    positions= [offset for offset, pos in enumerate(all_kmers) if Pattern == pos]
        
    return positions


### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
lines = sys.stdin.read().splitlines()
print(' '.join([str(i) for i in PatternMatching(lines[0],lines[1])]))
