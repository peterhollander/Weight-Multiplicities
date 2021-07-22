import itertools
from updated_alt_sets import currentAltSets2

#getting coefficients
bigToSmall = {
    "A" : {"a", "d", "j"}, "B" : {"b", "d", "j"}, "C" : {"a", "e", "j"},
    "D" : {"a", "d", "l"}, "E" : {"c", "e", "j"}, "F" : {"b", "f", "j"},
    "G" : {"a", "g", "l"}, "H" : {"b", "d", "l"}, "I" : {"a", "e", "o"},
    "J" : {"c", "f", "j"}, "K" : {"b", "h", "l"},
    "L" : {"a", "i", "o"}, "M" : {"b", "f", "p"}, "N" : {"c", "e", "o"},
    "O" : {"a", "g", "q"}, "P" : {"c", "f", "p"}, "Q" : {"a", "i", "q"},
}

#contradictions of type II
zeroToOne = {
    "b":{"a"}, "c":{"a", "b", "f"}, "e":{"d"}, "f":{"d", "e"},
    "g":{"d","e"}, "h":{"d","e","f","g"}, "i":{"d","e","g"}, "l":{"j"},
    "o":{"j","l"}, "p":{"j","l","o"}, "q":{"j","l","o","i"}, "a":set(),
    "d": set(), "j": set(), "p": set(),
}

allElementList = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q"]
allElementSet = set(allElementList)

#turns any imputed set of big letters into small ones
def takeBigToSmall(combination):
    allSmall = set()
    for element in combination:
        allSmall = allSmall.union(bigToSmall[element])
        allSmall.union(bigToSmall[element]))
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
    if {"a","g"}.issubset(allSmall):
        smallContradictions = smallContradictions.union("t")
        smallContradictions.union(zeroToOne[small]))

    return smallContradictions

def giveNeeded(combination): #combination = an alternation set
    #getting togher all small letters that need to be >= 0
    all_contradictions = giveContradictionsI(combination).union(giveContradictionsII(combination))
    
    #gettin all 3-element subsets of all_contradictions
    contraSubsets = list(itertools.combinations(all_contradictions,3))

    needed = []
    #turning tuples into sets
    contraChangedSubsets = list(set(subset) for subset in contraSubsets)

    #check what big letters need to be in the alt set
    for subset in contraChangedSubsets:
        for element in allElementList:
            if subset == bigToSmall[element] and element not in combination:
                needed.append(element)
    #account for s_0 and o_0
    if all(x in takeBigToSmall(combination) for x in ['o','s']):
        needed.append('K') #because K is a no-longer used letter
    return needed

def removeSubsets(altSets):
    newAltSet = []
    for subset in altSets:
        if len(giveNeeded(subset)) == 0:
            if not(all(x in subset for x in ['A','N']) and 'J' not in subset):
                newAltSet.append(subset)
    print(len(newAltSet))
    return newAltSet

print(removeSubsets(currentAltSets2))
