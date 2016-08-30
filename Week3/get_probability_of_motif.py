import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    filename= 'motifs.txt'
    target_motif= 'AAGTTC'

    #motifs= genemanip.openfile(filename)
    #profile= genemanip.generate_profile_matrix(motifs.splitlines())
    profile= {'A':[0.4,0.3,0.0,0.1,0.0,0.9], 'C':[0.2,0.3,0.0,0.4,0.0,0.1], 'G':[0.1,0.3,1.0,0.1,0.5,0.0],  'T':[0.3,0.1,0.0,0.4,0.5,0.0]}
    print profile
    prob= genemanip.get_probability_motif(target_motif, profile)
    print prob
