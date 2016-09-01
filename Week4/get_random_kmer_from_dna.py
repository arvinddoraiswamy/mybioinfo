import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'motifs.txt'
    dna= genemanip.openfile(filename)
    kmer_length= 3
    print genemanip.get_random_kmer_from_dna(dna.splitlines(), kmer_length, len(dna[:-1]))
