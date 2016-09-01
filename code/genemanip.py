from __future__ import division
import itertools
import random
import math

#Everything below this was all up until Week4
'''
This uses the 2 functions below to try and choose the best motif for each DNA string. It then uses this list as an input to the profile function and calculates its profile. Then it uses this profile to find the most probable kmers. It scores this new list against the best motif list we chose at the start. The moment the new list has a worse score than the previous list, we stop the search. If we don't the program will never end. 
'''
def randomized_motif_search(dna, kmer_length, no_of_motifs, no_of_iterations):
    count= 0

    #Randomly select a string of a certain length from random offsets in every string. Set this set of randomly chosen strings to be the Best Motif
    current_motif= get_random_kmer_from_dna(dna, kmer_length, no_of_motifs)
    best_motifs= current_motif

    while count < no_of_iterations:
        #Generate a profile matrix for the randomly selected kmers
        profile= generate_profile_matrix_with_pseudocounts(current_motif)

        #Use this profile to select a new set of most probable kmers from DNA
        current_motif= get_most_probable_kmers_from_profile(profile, dna)

        #The lower the score, the better the motif is, since score is a result of all the mismatches and you want more matches
        if score(current_motif) < score(best_motifs):
            best_motifs= current_motif
        else:
            return best_motifs
        count += 1

'''
Given a set of DNA strings and a kmer length, choose a random string of length kmer_length from each of the strings and store those
'''
def get_random_kmer_from_dna(dna, kmer_length, no_of_motifs):
    random_kmer= []
    for each_string in dna:
        offset= random.randint(1, len(each_string) - kmer_length)
        random_kmer.append(each_string[offset:offset+kmer_length])
    return random_kmer

'''
Given a preset profile and a bunch of Dna strings we start building towards implementing the RandomizedSearchAlgorithm.
'''
def get_most_probable_kmers_from_profile(profile, dna):
    most_probable_kmers= []
    for each_string in dna:
        most_probable_kmers.append(get_most_probable_motif(each_string, len(profile['A']), profile))
    return most_probable_kmers

'''
The basic idea of the greedy search is that you want to find a series of motifs (1 per string of DNA) that resembles all the others. Use the first motif as the fixed set. For each kmer in the first motif, go through all the other kmers in all the other DNA strings and find 1 string that matches the first one the most. Once you finish going through all the DNA strings, score the entire set. Now go to the next kmer in the first set and repeat this process. Compare the score of all the sets at the end. The "set" with the best score wins. Meaning, for each string of DNA that string was the most likely to be the motif.
'''
def greedysearch_with_pseudocounts(motifs, kmer_length, no_of_motifs):
    #Get the first kmer-length strings of all the motifs. These are the starting best motifs before we do any analysis. This is what will get updated as we move
    #through the string.
    best_motifs= []
    for motif_num in range(no_of_motifs):
        best_motifs.append(motifs[motif_num][0:kmer_length])

    #Get all the kmers of kmer length from the first motif
    all_kmers= find_kmers_in_genome(motifs[0], kmer_length)

    #Iterate through each kmer
    for kmer_num,kmer in enumerate(all_kmers):
        #Set each kmer from the first DNA string to be current. It is against this that you score all the other kmers from all the other DNA strings.
        current_motif= []
        current_motif.append(kmer)

        #Goes over all the DNA strings, finding the kmer_per_string that matches the kmer from the first string. This creates a complete set which is then scored.
        for j in range(1, no_of_motifs):
            profile = generate_profile_matrix_with_pseudocounts(current_motif[0:j])
            current_motif.append(get_most_probable_motif(motifs[j], kmer_length, profile))

        #If the score for the current motif is better, it becomes the new 'preferred set' and we move on to the next kmer from the first set after this.
        if score(current_motif) < score(best_motifs):
            best_motifs = current_motif

    return best_motifs

'''
The new profile matrix needs to be generated from the pseudocount matrix. Ensure that you calculate the correct per_column_totals before calculating profiles. The sum total of each column must always be 1, remember - hence the adjustment.
'''
def generate_profile_matrix_with_pseudocounts(motifs):
    len_each_motif= len(motifs[0])
    no_of_motifs= len(motifs)
    symbols=['A','C','G','T']
    profile= {}

    for symbol in symbols:
        profile[symbol]= []
        for symbol_cell in range(len_each_motif):
            profile[symbol].append(0)

    count_matrix= generate_count_matrix_with_pseudocounts(motifs)
    total= []
    total = [ sum(col) for col in itertools.izip_longest(*count_matrix.values(), fillvalue=0) ]

    for symbol,nucleotide_count_list in count_matrix.items():
        for symbol_count in range(len(nucleotide_count_list)):
            profile[symbol][symbol_count]= count_matrix[symbol][symbol_count] / total[symbol_count]

    return profile

'''
If any of the cells in the count matrix have a 0 in them, then any guessed motifs despite being really close could get rejected. So we increment each cell by 1 so nothing is ever 0.
'''
def generate_count_matrix_with_pseudocounts(motifs):
    symbols=['A','C','G','T']
    count  ={}
    final= {}
    len_each_motif= len(motifs[0])

    for symbol in symbols:
        count[symbol]= []
        for symbol_cell in range(len_each_motif):
            count[symbol].append(0)

    for motif in motifs:
        tempvar_2= list(motif)
        for offset in range(len(tempvar_2)):
            count[tempvar_2[offset]][offset] += 1

    for key,value_list in count.items():
        final[key]= []
        for value in value_list:
            final[key].append(value+1)

    return final

#Everything below this was all up until Week3
'''
Given a string of motifs extracted from a gene identify the number of A,G,C,T in each column and generate a matrix with those numbers as a dictionary of lists. Each entry in the list is a column.
'''
def generate_count_matrix(motifs):
    symbols=['A','C','G','T']
    count  ={}
    len_each_motif= len(motifs[0])

    for symbol in symbols:
        count[symbol]= []
        for symbol_cell in range(len_each_motif):
            count[symbol].append(0)

    for motif in motifs:
        tempvar_2= list(motif)
        for offset in range(len(tempvar_2)):
            count[tempvar_2[offset]][offset] += 1

    return count

'''
Given a count matrix, find out how much A,C,G and T occur in each column. Each column must add up to 1.
'''
def generate_profile_matrix(motifs):
    no_of_motifs= len(motifs)
    len_each_motif= len(motifs[0])
    symbols=['A','C','G','T']
    profile= {}

    for symbol in symbols:
        profile[symbol]= []
        for symbol_cell in range(len_each_motif):
            profile[symbol].append(0)

    count_matrix= generate_count_matrix(motifs)
    for symbol,nucleotide_count_list in count_matrix.items():
        for symbol_count in range(len(nucleotide_count_list)):
            profile[symbol][symbol_count]= count_matrix[symbol][symbol_count] / no_of_motifs

    return profile

'''
Which of A,C,G or T appeared the most in each column? Extract just that one nucleotide from each column. Join them all to get one of the
most likely consensuses. It is possible that there are multiple valid answers here if 2 nucleotides in a column have the same count.
'''
def consensus(motifs):
    count_matrix= generate_count_matrix(motifs)
    len_each_motif= len(motifs[0])
    symbols=['A','C','G','T']

    t1= 0
    predicted_motif= ''
    for position in range(len_each_motif):
        t1= 0
        t2= ''
        for symbol in symbols:
            if count_matrix[symbol][position] > t1:
                t1= count_matrix[symbol][position] 
                t2= symbol
        predicted_motif += t2

    return predicted_motif

'''
The consensus gives you the most_likely nucleotide per column. This function looks down the rest of the column and identifies the entries that do not match. This is repeated for_all_columns and for_all_genes, incrementing the score for each non match.
'''
def score(motifs):
    predicted_motif= consensus(motifs)
    score= 0
    for offset in range(0, len(motifs[0])):
        for motifnum in range(0, len(motifs)):
            if motifs[motifnum][offset] != predicted_motif[offset]:
                score += 1
            else:
                continue

    return score

'''
Given a profile and a target motif what is the probability that this motif is the one that actually exists in a genome.
'''
def get_probability_motif(target_motif, profile):
    probability= []
    for position in range(len(profile['A'])):
        if target_motif[position] in profile.keys():
            probability.append(profile[target_motif[position]][position])
        else:
            continue

    final_value= probability[0]
    for position in range(1,len(probability)):
        final_value *= probability[position]
    return final_value

'''
Given a list of genomes, kmer length and a profile, identify what the most likely motif is. In other words, you calculate a number of
probabilities and choose the one with the highest value.
'''
def get_most_probable_motif(genome, kmer_length, profile):
    all_kmers= find_kmers_in_genome(genome, kmer_length)
    all_probabilities= {}
    existing_values_list= []

    for kmer in all_kmers:
        value= get_probability_motif(kmer, profile)
        if value not in existing_values_list:
            existing_values_list.append(value)
            if kmer not in all_probabilities:
                all_probabilities[kmer]= value

    sorted_probs= sorted(all_probabilities.items(), key=lambda x: x[1])
    return sorted_probs[-1][0]

'''
The basic idea of the greedy search is that you want to find a series of motifs (1 per string of DNA) that resembles all the others. Use the first motif as the fixed set. For each kmer in the first motif, go through all the other kmers in all the other DNA strings and find 1 string that matches the first one the most. Once you finish going through all the DNA strings, score the entire set. Now go to the next kmer in the first set and repeat this process. Compare the score of all the sets at the end. The "set" with the best score wins. Meaning, for each string of DNA that string was the most likely to be the motif.
'''
def greedysearch(motifs, kmer_length, no_of_motifs):
    #Get the first kmer-length strings of all the motifs. These are the starting best motifs before we do any analysis. This is what will get updated as we move
    #through the string.
    best_motifs= []
    for motif_num in range(no_of_motifs):
        best_motifs.append(motifs[motif_num][0:kmer_length])

    #Get all the kmers of kmer length from the first motif
    all_kmers= find_kmers_in_genome(motifs[0], kmer_length)

    #Iterate through each kmer
    for kmer_num,kmer in enumerate(all_kmers):
        current_motif= []
        current_motif.append(kmer)

        #Goes over all the DNA strings and adds to the list, each time calculating profiles and then the probability of that profile in the next string.
        for j in range(1, no_of_motifs):
            profile = generate_profile_matrix(current_motif[0:j])
            current_motif.append(get_most_probable_motif(motifs[j], kmer_length, profile))
        if score(current_motif) < score(best_motifs):
            best_motifs = current_motif

    return best_motifs

def calculate_entropy(profile):
    len_each_motif= len(profile['A'])
    symbols=['A','C','G','T']
    logarithm_base= 2
    values= []
    for position in range(len_each_motif):
        for symbol in symbols:
            if profile[symbol][position] != 0.0:
                values.append(profile[symbol][position] * math.log(profile[symbol][position], logarithm_base))
            else:
                values.append(0.0)

    entropy= abs(sum(values))
    return entropy

#Everything below this was all up until Week2
'''
Get approximate matches for a specific pattern in a genome. The number of differences should not be greater than max_number_of_mismatches
'''
def approximate_match(genome, pattern, max_number_of_mismatches):
    positions= []
    all_kmers= find_kmers_in_genome(genome, len(pattern))
    for offset,kmer in enumerate(all_kmers):
        distance= hamming_distance(kmer, pattern)
        if distance <= max_number_of_mismatches:
            positions.append(offset)
        else:
            continue

    return positions

'''
Get hamming distance between 2 strings
'''
def hamming_distance(string1, string2):
    tempvar_1= list(string1)
    tempvar_2= list(string2)

    distance= 0
    for count, char in enumerate(tempvar_1):
        if char != tempvar_2[count]:
            distance += 1
        else:
            continue

    return distance

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
        
#Everything below this was all up until Week1
'''
This is a simple utility function that takes a file as an input and returns a list of strings that can then be worked upon.
'''
def openfile(filename):
    with open(filename) as f:
        genome= f.read()
    return genome

'''
This function takes a genome as an input and returns all patterns of a specific length in it.
'''
def find_kmers_in_genome(genome, k):
    all_kmers= []
    len_kmer= k
    len_genome= len(genome)
    len_diff= len_genome - len_kmer + 1

    all_kmers= [genome[i:i+len_kmer] for i in range(0, len_diff)]
    return all_kmers

'''
This takes a bunch of kmers of a certain length and searches for the most frequent one.
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
This gets all the offsets for a specific pattern in a genome.
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

