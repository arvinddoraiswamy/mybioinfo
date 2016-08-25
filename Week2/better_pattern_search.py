import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'E_coli.txt'
    pattern = 'C'

    genome= genemanip.openfile(filename)
    pattern_map= genemanip.search_for_pattern(genome, pattern)
    print min(pattern_map.items(), key=lambda x: x[1])
    print max(pattern_map.items(), key=lambda x: x[1])
