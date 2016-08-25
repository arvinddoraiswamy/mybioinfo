'''
def find_oris(string):
    l1= len(string)
    l2= len(pattern)
    diff= l1-l2+1

    offsets= []

    for i in range(0, diff):
        t1= string[i:i+l2]
        if t1.lower()== pattern.lower():
            offsets.append(i)
    
    return offsets

if __name__ == "__main__":
    string = 'AAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAATAATTACAGAGTACACAACATCCAT'
    pattern= 'AAA'
    offsets= find_oris(string)
    print(len(offsets))
'''

# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
lines = sys.stdin.read().splitlines()
print(PatternCount(lines[1],lines[0]))
