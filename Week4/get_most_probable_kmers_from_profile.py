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
    profile= {'A':[0.8,0.0,0.0,0.2],'C':[0.0,0.6,0.2,0.0], 'G':[0.2,0.2,0.8,0.0], 'T':[0.0,0.2,0.0,0.8]}
    genemanip.get_most_probable_kmers_from_profile(profile, dna.splitlines())

