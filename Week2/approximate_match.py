import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    pattern = 'CCGTCATCC'
    genome  = 'CCGTCATCCGTCATCCTCGCCACGTTGGCATGCATTCCGTCATCCCGTCAGGCATACTTCTGCATATAAGTACAAACATCCGTCATGTCAAAGGGAGCCCGCAGCGGTAAAACCGAGAACCATGATGAATGCACGGCGATTGC'
    max_number_of_mismatches = 3

    positions= genemanip.approximate_match(genome, pattern, max_number_of_mismatches)
    print positions
