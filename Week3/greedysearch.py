import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'DosR.txt'
    motifs= genemanip.openfile(filename).splitlines()
    no_of_motifs= len(motifs)
    kmer_length= 15
    Motifs= genemanip.greedysearch(motifs, kmer_length, no_of_motifs)
    print Motifs
    print genemanip.score(Motifs)
