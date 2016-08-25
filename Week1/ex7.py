# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    FrequentPatterns = [] # output variable
    # your code here
    d1= CountDict(Text,k)
    all_offsets= sorted(d1.values(), reverse=True)
    m1= max(all_offsets)

    for k1,v1 in d1.items():
        t1= Text[k1:k1+k]
        if v1 == m1 and t1 not in FrequentPatterns: 
            FrequentPatterns.append(t1)
        else:
            continue
    
    return FrequentPatterns

# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
# HINT:   This code should be identical to when you last implemented CountDict
def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Pattern, Text): 
    count = 0
    for i in range(len(Text)-len(Pattern)+1): 
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys
lines = sys.stdin.read().splitlines()
print(' '.join(FrequentWords(lines[0],int(lines[1]))))
