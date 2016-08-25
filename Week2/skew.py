import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    genome= 'GATACACTTCCCGAGTAGGTACTG'
    pattern1 = 'C'
    pattern2 = 'G'
    skew= {}
    skew= genemanip.get_skew_array(genome, pattern1, pattern2)
    print max(skew.items(), key=lambda x: x[1])

    #for key,value in skew.items():
    #    print key, value
