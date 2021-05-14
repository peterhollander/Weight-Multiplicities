import itertools

bigToSmall = {
    "A" : {"a", "e", "t"}, "B" : {"b", "e", "t"}, "C" : {"a", "f", "t"},
    "D" : {"a", "e", "l"}, "E" : {"c", "f", "t"}, "F" : {"b", "g", "t"},
    "G" : {"a", "h", "l"}, "H" : {"b", "e", "l"}, "I" : {"a", "f", "r"},
    "J" : {"c", "g", "t"}, "K" : {"d", "h", "l"}, "L" : {"b", "i", "l"},
    "M" : {"a", "j", "r"}, "N" : {"b", "g", "s"}, "O" : {"c", "f", "r"},
    "P" : {"a", "h", "o"}, "Q" : {"c", "g", "s"}, "R" : {"a", "j", "o"},
}
allElementList = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                  "M", "N", "O", "P", "Q", "R"]
allElementSet = set(allElementList)

def giveAltSet():
    altSet = []
    setCombinations = [] #holds all possible alt sets, before taking things out
    #get all possible alternation sets (subsets of allElementList)
    for length in range(len(allElementSet)+1):
        setCombinations.extend(list(itertools.combinations(allElementSet,length)))
        #turn subsets into lists
    for i in range(len(setCombinations)):
        setCombinations[i] = list(setCombinations[i])

    #workin on a specific subset
    for combination in setCombinations:
        allElements = set() #holds all small elements (ie. conditions) from combination
        neededElements = [] # holds elements that are needed
        allElementsSubsets = [] #holds all 3-element substets of allElements

        #put all of the small elements (conditions) together
        for element in combination:
            allElements = allElements.union(set(bigToSmall[element]))

        #now we make all possible 3-element subsets of the conditions
        allElementsSubsets = list(itertools.combinations(allElements,3))
        #and turn subsets from tuples into sets
        for i in range(len(allElementsSubsets)):
            allElementsSubsets[i] = set(allElementsSubsets[i])

        #now we use the 3-element subsets to see what elements are needed
        for subset in allElementsSubsets:
            for allElement in allElementList:
                #check if a big letter equals it and if it is
                #one of the letters already in alt set
                    if subset == bigToSmall[allElement] and allElement not in combination:
                        #add the big letter to a list
                        neededElements.append(allElement)

        #if there are no needed elements, ie this is a valid altset, we add it
        if len(neededElements) == 0:
            altSet.append(combination)
    print("number of sets:", len(altSet))
    return altSet

print(giveAltSet())
