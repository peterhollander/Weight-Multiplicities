import itertools
from updated_alt_sets import currentAltSets2

#getting coefficients
bigToSmall = {
    "A" : {"a", "e", "t"}, "B" : {"b", "e", "t"}, "C" : {"a", "f", "t"},
    "D" : {"a", "e", "l"}, "E" : {"c", "f", "t"}, "F" : {"b", "g", "t"},
    "G" : {"a", "h", "l"}, "H" : {"b", "e", "l"}, "I" : {"a", "f", "r"},
    "J" : {"c", "g", "t"}, "K" : {"d", "h", "l"}, "L" : {"b", "i", "l"},
    "M" : {"a", "j", "r"}, "N" : {"b", "g", "s"}, "O" : {"c", "f", "r"},
    "P" : {"a", "h", "o"}, "Q" : {"c", "g", "s"}, "R" : {"a", "j", "o"},
}

#contradictions of type II
zeroToOne = {
    "b":{"a"}, "c":{"a", "b"}, "d":{"a", "b", "c"}, "f":{"e"}, "g":{"e", "f"},
    "h":{"e","f"}, "i":{"e","f","g","h"}, "j":{"e","f","h"}, "l":{"t"},
    "r":{"t","l"}, "s":{"t","l","r"}, "o":{"t","l","r"}, "a":set(),
    "e": set(), "t": set(), "p": set(),
}

#coefficients,, just in case (idk if it's used yet)
A = {"a", "e", "t"} 
D = {"a", "e", "l"} 
C = {"a", "f", "t"} 
I = {"a", "f", "r"}
G = {"a", "h", "l"}
P = {"a", "h", "o"}
M = {"a", "j", "r"}
R = {"a", "j", "o"}
B = {"b", "e", "t"}
H = {"b", "e", "l"}
F = {"b", "g", "t"}
N = {"b", "g", "s"}
L = {"b", "i", "l"}
E = {"c", "f", "t"}
O = {"c", "f", "r"}
J = {"c", "g", "t"}
Q = {"c", "g", "s"}
K = {"d", "h", "l"}

allElementList = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                  "M", "N", "O", "P", "Q", "R"]
allElementSet = set(allElementList)


#turns any imputed set of big letters into small ones
def takeBigToSmall(combination):
    allSmall = set()
    for element in combination:
        #print("element", element)
        #print("bigToSmall[element]", bigToSmall[element])
        allSmall = allSmall.union(bigToSmall[element])
        #print("allSmall.union(bigToSmall[element])", allSmall.union(bigToSmall[element]))
    #print("combination", combination)
    #print("allSmall",allSmall)
    return allSmall



#gives small letters that need to be >=0 or there's a type I contradiction
def giveContradictionsI(combination):
    return takeBigToSmall(combination)


#gives small letters that need to be >= 0 or there's a type II contradiction
def giveContradictionsII(combination):
    allSmall = takeBigToSmall(combination)

    #we are assuming this combination is included
    smallContradictions = set() #these also need to be >=0
    for small in allSmall:
        smallContradictions = smallContradictions.union(zeroToOne[small])
        #print("smallContradictions.union(zeroToOne[small])", smallContradictions.union(zeroToOne[small]))
    #print("smallContradictions", smallContradictions)

    return smallContradictions


def giveNeeded(combination): #combination = an alternation set
    #getting togher all small letters that need to be >= 0
    all_contradictions = giveContradictionsI(combination).union(giveContradictionsII(combination))
    
    #gettin all 3-element subsets of all_contradictions
    contraSubsets = list(itertools.combinations(all_contradictions,3))
    #contraSubsets = list(itertools.combinations(smallContradictions,3))

    needed = []
    #turning tuples into sets
    contraChangedSubsets = list(set(subset) for subset in contraSubsets)
    #print("contrasubsets", contraChangedSubsets)

    #check what big letters need to be in the alt set
    for subset in contraChangedSubsets:
        for element in allElementList:
            if subset == bigToSmall[element] and element not in combination:
            #if subset == bigToSmall[element]:
                needed.append(element)
    return needed


def removeSubsets(altSets):
    newAltSet = []
    for subset in altSets:
        if len(giveNeeded(subset)) == 0:
            newAltSet.append(subset)
    print(len(newAltSet))
    return newAltSet

print(removeSubsets(currentAltSets2))

