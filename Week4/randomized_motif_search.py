import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'DosR.txt'
    dna= genemanip.openfile(filename).splitlines()
    kmer_length= 15
    no_of_iterations= 100
    print genemanip.randomized_motif_search(dna, kmer_length, len(dna), no_of_iterations)

