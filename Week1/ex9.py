# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):
    revComp = '' # output variable
    # your code here
    comp= complement(Pattern)
    revComp= ''.join(comp)[::-1].upper()
    return revComp


# Copy your reverse function from the previous step here.


# HINT:   Filling in the following function is optional, but it may come in handy when solving ReverseComplement
# Input:  A character Nucleotide
# Output: The complement of Nucleotide
def complement(Nucleotide):
    nucleotide_map= {'a':'t', 't':'a', 'c':'g', 'g':'c'}
    comp = '' # output variable
    # your code here
    split_Nucleotide= list(Nucleotide)
    comp= [nucleotide_map[char.lower()] for char in split_Nucleotide if char.lower() in nucleotide_map.keys()]
    return comp



### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
print(ReverseComplement(sys.stdin.read().strip()))
