'''
Get minimum count from a pre-created skew array for a given genome and specific patterns
'''
def min_skew_array(genome, pattern1, pattern2):
    positions= []
    skew= get_skew_array(genome, pattern1, pattern2)
    min_value= min(skew.values())
    sorted_skew= sorted(skew.items(), key=lambda x: x[1])
    positions= [tempvar_1[0] for tempvar_1 in sorted_skew if tempvar_1[1] == min_value]
    return positions

'''
Get skew array for a given genome and a couple of specific patterns
'''
def get_skew_array(genome, pattern1, pattern2):
    temp_var1= list(genome)
    pattern1= pattern1.lower()
    pattern2= pattern2.lower()
    skew_count= 0
    skew= {}
    skew[skew_count]= 0
    for count, char in enumerate(temp_var1, start=1):
        char= char.lower()
        if char == pattern1:
            skew_count -= 1
            skew[count]= skew_count
        elif char == pattern2:
            skew_count += 1
            skew[count]= skew_count
        else:
            skew[count]= skew_count

    return skew

'''
Search for a single nucleotide
'''
def search_for_pattern(genome, pattern):
    pattern_map = {}
    genome= genome.lower()
    pattern= pattern.lower()
    original_len_genome= len(genome)
    count= 0

    #First do some work on the first half of the genome. Once done, add it to the end and do the rest of the stuff.
    first_half_genome= genome[0:(original_len_genome/2)]
    for char in list(first_half_genome):
        if char == pattern:
            count += 1
        else:
            continue
    pattern_map[0]= count

    #Start work by sliding windows bit by bit
    window_start= 1
    window_end= original_len_genome/2
    genome= genome + first_half_genome

    while window_start < original_len_genome: 
        if genome[window_start-1] == pattern and count != 0:
            count -= 1
        if genome[window_end] == pattern:
            count += 1

        pattern_map[window_start] = count
        window_start += 1
        window_end   += 1

    return pattern_map
        
'''
This is a simple utility function that takes a file as an input and returns a list of strings that can then be worked upon.
'''
def openfile(filename):
    with open(filename) as f:
        genome= f.read()
    return genome

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
    all_kmers= find_kmers_in_genome(genome, len(pattern))
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

