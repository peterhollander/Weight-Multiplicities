import itertools
from updated_alt_sets import currentAltSets1

bigToSmall = {
    "A" : {"a", "e", "t"}, "B" : {"b", "e", "t"}, "C" : {"a", "f", "t"},
    "D" : {"a", "e", "l"}, "E" : {"c", "f", "t"}, "F" : {"b", "g", "t"},
    "G" : {"a", "h", "l"}, "H" : {"b", "e", "l"}, "I" : {"a", "f", "r"},
    "J" : {"c", "g", "t"}, "L" : {"b", "i", "l"},
    "M" : {"a", "j", "r"}, "N" : {"b", "g", "s"}, "O" : {"c", "f", "r"},
    "P" : {"a", "h", "o"}, "Q" : {"c", "g", "s"}, "R" : {"a", "j", "o"},
}
#contradictions
zeroToOne = {
    "b":{"a"}, "c":{"a", "b","g"}, "d":{"a", "b", "c"}, "f":{"e"}, "g":{"e", "f"},
    "h":{"e","f"}, "i":{"e","f","g","h"}, "j":{"e","f","h"}, "l":{"t"},
    "r":{"t","l"}, "s":{"t","l","r"}, "o":{"t","l","r","j"}, "a":set(),
    "e": set(), "t": set(), "p": set(),
    }
allElementList = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "L", "M", "N", "O",    "P", "Q", "R"]

#turns any imputed set of big letters (summation term) into small ones (conditions)
def takeBigToSmall(combination):
    allSmall = set()
    for element in combination:
        allSmall = allSmall.union(bigToSmall[element])
    return allSmall

#gives small letters (conditions) that need to be >= 0 or there's a contradiction
def giveContradictions(combination):
    allSmall = takeBigToSmall(combination)
    #we are assuming this combination is included in the alt set
    smallContradictions = set() #these also need to be >=0
    for small in allSmall:
        smallContradictions = smallContradictions.union(zeroToOne[small])
    if {"a","g"}.issubset(allSmall):
        smallContradictions = smallContradictions.union("t")
    return smallContradictions

def giveNeeded(combination):
    #gettin all 3-element subsets of the 'contradicting' small leters (conditions)
    contraSubsets = list(itertools.combinations(giveContradictions(combination),3))
    needed = []
    #turning tuples into sets
    contraChangedSubsets = list(set(subset) for subset in contraSubsets)
    #check what big letters (summation terms) are needed
    for subset in contraChangedSubsets:
        for element in allElementList:
            if subset == bigToSmall[element] and element not in combination:
                needed.append(element)
    #account for s_0 and o_0
    if all(x in takeBigToSmall(combination) for x in ['o','s']):
        needed.append('K') #because K is a no-longer used letter
    return needed

#put this all together to discard contradicting alt sets
def removeSubsets(altSets):
    newAltSet = []
    for subset in altSets:
        if len(giveNeeded(subset)) == 0:
            if not(all(x in subset for x in ['A','N']) and 'J' not in subset):
                newAltSet.append(subset)
    print(len(newAltSet))
    return newAltSet

print(removeSubsets(currentAltSets1))
