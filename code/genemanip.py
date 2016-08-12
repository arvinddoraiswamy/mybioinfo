'''
This is a simple utility function that takes a file as an input and returns a list of strings that can then be worked upon.
'''
def openfile(filename):
    file_content= []
    with open(filename) as f:
        tempvar_1= f.readlines()

    for string in tempvar_1:
        file_content.append(string)

    return file_content

'''
This function takes a genome as an input and finds the most frequent pattern of a specific length in it.
'''
def find_kmers_in_genome(genome, k):
    all_kmers= []
    unique_kmers= []
    kmer_map= {}

    len_kmer= k
    len_genome= len(genome)
    len_diff= len_genome - len_kmer + 1
    all_kmers= [string[i:i+len_kmer] for i in range(0, len_diff)]
    return all_kmers

'''
This takes a bunch of kmers and searches for the most frequent one.
'''
def get_most_freq_kmer(all_kmers):
    for kmer in all_kmers:
        if kmer in unique_kmers:
            kmer_map[kmer] += 1
        else:
            unique_kmers.append(kmer)
            kmer_map[kmer] = 1

    #Get the most frequent kmer
    kmer_count_by_value= sorted(kmer_map.values())
    kmer_most_freq_key= [key for key,value in kmer_map.items() if kmer_count_by_value[-1] == value]
    return kmer_most_freq_key

'''
This gets all the offsets for a specific pattern.
'''
def find_offsets_for_pattern(all_kmers, pattern):
    offsets= []
    all_kmers= find_oris(pattern, genome)
    positions= [offset for offset, kmer in enumerate(all_kmers) if pattern == kmer]
        
    return positions

'''
Gets a complement of a specific pattern
'''
def complement(pattern):
    nucleotide_map= {'a':'t', 't':'a', 'c':'g', 'g':'c'}
    comp = ''
    tempvar_split_pattern= list(pattern)
    comp= [nucleotide_map[char.lower()] for char in tempvar_split_pattern if char.lower() in nucleotide_map.keys()]
    return comp

'''
Reverses a pattern
'''
def reverse_complement(pattern):
    rev_comp = ''
    comp= complement(pattern)
    rev_comp= ''.join(comp)[::-1].upper()
    return rev_comp
