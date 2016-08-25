# Copy your PatternMatching function below this line.

# Call PatternMatching with Pattern equal to "CTTGATCAT" and Genome equal to v_cholerae,
# and store the output as a variable called positions
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

# The following lines will automatically read in the Vibrio cholerae genome for you and store it in a variable named v_cholerae
import sys                              # needed to read the genome
with open('Vibrio_cholerae.txt','r') as f:
#with open('t_petrophila_oriC.txt','r') as f:
    text= f.read()

patterns= ['ATGATCAAG', 'CTTGATCAT']
#Pattern='ATGATCAAG'
#Pattern='CTTGATCAT'
total= 0
for Pattern in patterns:
    positions= PatternMatching(Pattern, text)
    total += len(positions)
    print(total)
