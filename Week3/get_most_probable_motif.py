import sys
import os

#Adding directory to the path where Python searches for modules
module_folder = os.path.dirname('/opt/Courses/TechCourses/Bioinformatics/code/')
sys.path.insert(0, module_folder)

#Importing genemanipulating module. This has a lot of common functions.
import genemanip

if __name__ == "__main__":
    import sys
    lines = sys.stdin.read().splitlines()
    Text = lines[0]
    k = int(lines[1])
    A = [float(c) for c in lines[2].split()]
    C = [float(c) for c in lines[3].split()]
    G = [float(c) for c in lines[4].split()]
    T = [float(c) for c in lines[5].split()]
    profile = {'A':A, 'C':C, 'G':G, 'T':T}
    #profile= genemanip.generate_profile_matrix(motifs.splitlines())
    print genemanip.get_most_probable_motif(Text, k, profile)
