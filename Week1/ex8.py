def find_oris(string):
    diff= len_string - len_pattern + 1
    all_2mers= [string[i:i+len_pattern] for i in range(0, diff)]
    return all_2mers

if __name__ == "__main__":
    string= 'ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'
    len_pattern= 10
    len_string = len(string)
    all_2mers= find_oris(string)

    unique_2mers= []
    kmer_map= {}
    kmer_map_offsets= {}
    offset= 0

    for kmer in all_2mers:
        if kmer in unique_2mers:
            kmer_map[kmer] += 1
        else:
            unique_2mers.append(kmer)
            kmer_map[kmer]   = 1

        kmer_map_offsets.setdefault(kmer, []).append(offset)
        offset += 1

    all_offsets= sorted(kmer_map_offsets.values(), key=len, reverse=True)
    len_all_offsets= [len(offset) for offset in all_offsets]
    
    import itertools
    final_offsets= [o1 for o1 in all_offsets if len(o1) == max(len_all_offsets)]
    l1= sorted(list(itertools.chain(*final_offsets)))

    words=[]
    for offset in l1:
        ori= string[offset:offset+len_pattern]
        if ori not in words:
            words.append(ori)
        else:
            continue
    print(words)
