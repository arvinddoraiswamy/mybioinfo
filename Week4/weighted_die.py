import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    probabilities= {'A': 0.2, 'C': 0.3, 'G': 0.4, 'T': 0.1}
    key= genemanip.weighted_die(probabilities)
    print key
