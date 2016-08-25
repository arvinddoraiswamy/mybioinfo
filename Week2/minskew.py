import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'E_coli.txt'
    genome= genemanip.openfile(filename)
    #genome= 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
    pattern1 = 'C'
    pattern2 = 'G'
    skew= {}
    positions= genemanip.min_skew_array(genome, pattern1, pattern2)

    for position in positions:
        print position
