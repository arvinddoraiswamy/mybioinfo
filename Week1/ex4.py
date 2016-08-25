def find_oris(string):
    len_pattern= 2
    len_string = len(string)
    diff= len_string - len_pattern + 1

    all_2mers= [string[i:i+len_pattern] for i in range(0, diff)]
    return all_2mers
        
if __name__ == "__main__":
    string= 'GATCCAGATCCCCATAC'
    all_2mers= find_oris(string)

    unique_2mers= []
    kmer_map= {}

    for kmer in all_2mers:
        if kmer in unique_2mers:
            kmer_map[kmer] += 1
        else:
            unique_2mers.append(kmer)
            kmer_map[kmer] = 1

    kmer_count_by_value= sorted(kmer_map.values())
    kmer_most_freq_key= [key for key,value in kmer_map.items() if kmer_count_by_value[-1] == value]
    print(kmer_most_freq_key[0])
