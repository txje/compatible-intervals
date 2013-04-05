# Jeremy Wang
# Oct 12, 2009
#
# Compute interval covers using heterozygous genotype data
# leftScan: non-overlapping left-to-right (simply reverse data to get R->L)
# uberScan: overlapping maximal-size intervals
# To compute Max-k intervals from Uber intervals, see compatinv.py


# returns l->r intervals for snpList (hap or genotype)
def leftScan(snpList, pessimistic=False):
    # 0 - incompatible, 1 - out of phase, 2 - compatible, 3 - in phase, 4 - either
    start = 0
    intervals = []
    inv_groupings = []
    recent_percent = 0.01
    while start < len(snpList)-1:
        if float(start)/len(snpList) > recent_percent:
            print str(recent_percent*100) + "%"
            recent_percent += 0.01
        genGroups = [[] for i in xrange(len(snpList[0]))]
        blueGroups = [[] for i in xrange(len(snpList[0]))]
        for end in range(start + 1, len(snpList)):
            incomp = False
            problem = ""
            for c in range(start, end):
                compat = compatible(snpList[c], snpList[end], pessimistic)
                if compat == 0:
                    problem = "red between " + str(c) + " and " + str(end)
                    incomp = True
                elif compat == 2:
                    continue
                elif compat == 1:
                    o = orangeComp(c, end, genGroups, snpList[c], snpList[end])
                    genGroups = o[1]
                    if not o[0]:
                        incomp = True
                        problem = "orange between " + str(c) + " and " + str(end)
                elif compat == 3:
                    g = greenComp(c, end, genGroups, snpList[c], snpList[end])
                    genGroups = g[1]
                    if not g[0]:
                        incomp = True
                        problem = "green between " + str(c) + " and " + str(end)
                elif compat == 4:
                    l = blueComp(c, end, blueGroups, snpList[c], snpList[end])
                    blueGroups = l[1]
                    if not l[0]:
                        incomp = True
                        problem = "blue between " + str(c) + " and " + str(end)
                b = checkBlue(genGroups, blueGroups, snpList)
                genGroups = b[1]
                blueGroups = b[2]
                if not b[0]:
                    incomp = True
                if incomp:
                    break
            if incomp:
                intervals.append((start,end-1))
                #print problem
                #removeSNP(genGroups, end) #SHIT!
                inv_groupings.append(removeSNP(genGroups, end))
                start = end
                break
        if not incomp or start == len(snpList)-1:
            intervals.append((start,end))
            #print problem
            if not incomp:
                inv_groupings.append(genGroups)
            else:
                inv_groupings.append([[] for i in xrange(len(snpList[0]))])
            start = end
    return intervals, inv_groupings

# checks and extends consistency of blue edges (imperfect and probably really slow)
def checkBlue(genGroups, blueGroups, snpList):
    for index in xrange(len(blueGroups)):
        for pair in blueGroups[index]:
            relationship = None
            oneloc, twoloc = locations(pair[0], pair[1], genGroups[index])
            if not oneloc == None and not twoloc == None and oneloc[0] == twoloc[0]:
                if oneloc[1] == twoloc[1]:
                    relationship = 0
                else:
                    relationship = 1
            if not relationship == None:
                for blueindex in xrange(len(blueGroups)):
                    if not pair in blueGroups[blueindex] or blueindex == index:
                        continue
                    compare_oneloc, compare_twoloc = locations(pair[0], pair[1], genGroups[blueindex])
                    if not compare_oneloc == None and not compare_twoloc == None and compare_oneloc[0] == compare_twoloc[0]: # it has already established a relationship, just check that it is consistent
                        if (relationship == 0 and not compare_oneloc[1] == compare_twoloc[1]) or (relationship == 1 and compare_oneloc[1] == compare_twoloc[1]):
                            return (False, genGroups, blueGroups)
                    else:
                        if compare_oneloc == None and compare_twoloc == None:
                            # add this relationship
                            if relationship == 0:
                                genGroups[blueindex].append([[pair[0], pair[1]],[]])
                            else: # relationship == 1
                                genGroups[blueindex].append([[pair[0]], [pair[1]]])
                        elif compare_oneloc == None and not compare_twoloc == None:
                            genGroups[blueindex][compare_twoloc[0]][compare_twoloc[1]-relationship].append(pair[0]) # funny logic, assuming index -1 == index 1 in a len(2) array
                        elif not compare_oneloc == None and compare_twoloc == None:
                            genGroups[blueindex][compare_oneloc[0]][compare_oneloc[1]-relationship].append(pair[1]) # funny logic, assuming index -1 == index 1 in a len(2) array
                        else: # not compare_oneloc[0] == compare_twoloc[0]:
                            # merge them
                            genGroups[blueindex][compare_oneloc[0]][compare_oneloc[1]].extend(genGroups[blueindex][compare_twoloc[0]][compare_twoloc[1]-relationship])
                            genGroups[blueindex][compare_oneloc[0]][1-compare_oneloc[1]].extend(genGroups[blueindex][compare_twoloc[0]][1-compare_twoloc[1]-relationship])
                            genGroups[blueindex].pop(compare_twoloc[0])
    return (True, genGroups, blueGroups)

# alters genGroups, adding a blue edge between c and end
def blueComp(c, end, blueGroups, csnp, endsnp):
    for index in range(len(blueGroups)):
        if csnp[index] == '2' and endsnp[index] == '2':
            blueGroups[index].append((c,end))
    return (True, blueGroups)

# alters genGroups, adding a green edge between c and end
def greenComp(c, end, genGroups, csnp, endsnp):
    for index in range(len(genGroups)):
        if csnp[index] == '2' and endsnp[index] == '2':
            cloc, endloc = locations(c, end, genGroups[index])
            if cloc == None and endloc == None:
                #create a new subgraph
                genGroups[index].append([[c, end],[]])
            elif cloc == None and not endloc == None:
                #c does not exist
                genGroups[index][endloc[0]][endloc[1]].append(c)
            elif not cloc == None and endloc == None:
                #end does not exist
                    genGroups[index][cloc[0]][cloc[1]].append(end)
            elif not cloc[0] == endloc[0]:
                #they are in separate existing subgraphs, merge them
                if cloc[1] == endloc[1]:
                    genGroups[index][cloc[0]][0].extend(genGroups[index][endloc[0]][0])
                    genGroups[index][cloc[0]][1].extend(genGroups[index][endloc[0]][1])
                else:
                    genGroups[index][cloc[0]][0].extend(genGroups[index][endloc[0]][1])
                    genGroups[index][cloc[0]][1].extend(genGroups[index][endloc[0]][0])
                genGroups[index].pop(endloc[0])
            elif cloc[0] == endloc[0]:
                #they are in the same subgraph
                if cloc[1] == endloc[1]:
                    pass
                else:
                    return (False, genGroups)
    return (True, genGroups)

# alters genGroups, adding an orange edge between c and end
def orangeComp(c, end, genGroups, csnp, endsnp):
    for index in range(len(genGroups)):
        if csnp[index] == '2' and endsnp[index] == '2':
            cloc, endloc = locations(c, end, genGroups[index])
            if cloc == None and endloc == None:
                #create a new subgraph
                genGroups[index].append([[c],[end]])
            elif cloc == None and not endloc == None:
                #c does not exist
                genGroups[index][endloc[0]][1-endloc[1]].append(c)
            elif not cloc == None and endloc == None:
                #end does not exist
                genGroups[index][cloc[0]][1-cloc[1]].append(end)
            elif not cloc[0] == endloc[0]:
                #they are in separate existing subgraphs, merge them (opposite)
                if cloc[1] == endloc[1]:
                    genGroups[index][cloc[0]][0].extend(genGroups[index][endloc[0]][1])
                    genGroups[index][cloc[0]][1].extend(genGroups[index][endloc[0]][0])
                else:
                    genGroups[index][cloc[0]][0].extend(genGroups[index][endloc[0]][0])
                    genGroups[index][cloc[0]][1].extend(genGroups[index][endloc[0]][1])
                genGroups[index].pop(endloc[0])
            elif cloc[0] == endloc[0]:
                #they are in the same subgraph
                if not cloc[1] == endloc[1]:
                    pass
                else:
                    return (False, genGroups)
    return (True, genGroups)

# adds c, end to genGroups and blueGroups, returns if they are compatible
def genCompatible(c, end, snpList, genGroups, blueGroups, pessimistic=False):
    compat = compatible(snpList[c], snpList[end], pessimistic)
    if compat == 0:
        return (False, genGroups, blueGroups)
    elif compat == 2:
        pass
    elif compat == 1:
        o = orangeComp(c, end, genGroups, snpList[c], snpList[end])
        genGroups = o[1]
        if not o[0]:
            return (False, genGroups, blueGroups)
    elif compat == 3:
        g = greenComp(c, end, genGroups, snpList[c], snpList[end])
        genGroups = g[1]
        if not g[0]:
            return (False, genGroups, blueGroups)
    elif compat == 4:
        l = blueComp(c, end, blueGroups, snpList[c], snpList[end])
        blueGroups = l[1]
        if not l[0]:
            return (False, genGroups, blueGroups)
    b = checkBlue(genGroups, blueGroups, snpList)
    genGroups = b[1]
    blueGroups = b[2]
    if not b[0]:
        return (False, genGroups, blueGroups)
    return (True, genGroups, blueGroups)

# returns uber intervals of snpList (haplotype or genotype)
def uberScan(snpList, pessimistic=False):
    # 0 - incompatible, 1 - out of phase, 2 - compatible, 3 - in phase, 4 - either
    start = 0
    intervals = []
    inv_groupings = []
    while start < len(snpList):
        genGroups = [[] for i in range(len(snpList[0]))]
        blueGroups = [[] for i in xrange(len(snpList[0]))]
        compatible = True
        for end in range(start + 1, len(snpList)):
            for c in range(start, end):
                genCompat = genCompatible(c, end, snpList, genGroups, blueGroups, pessimistic)
                compatible = genCompat[0]
                genGroups = genCompat[1]
                blueGroups = genCompat[2]
                if not compatible:
                    break
            if not compatible:
                break
        if compatible:
            break
        
        intervals.append((start,end-1))
        #removeSNP(genGroups, end) #SHIT!!!!!!
        inv_groupings.append(removeSNP(genGroups, end))
        
        # backtrack
        genGroups = [[] for i in range(len(snpList[0]))]
        blueGroups = [[] for i in xrange(len(snpList[0]))]
        new_start = end-1 # this is very important in case the last interval had only a single SNP
        for new_start in xrange(end-1, start-1, -1):
            for c in xrange(end, new_start, -1):
                genCompat = genCompatible(c, new_start, snpList, genGroups, blueGroups, pessimistic)
                compatible = genCompat[0]
                genGroups = genCompat[1]
                blueGroups = genCompat[2]
                if not compatible:
                    break
            if not compatible:
                break
        start = new_start+1
        
    intervals.append((start, end))
    inv_groupings.append(genGroups)
    return intervals, inv_groupings


# -----------------------------------------
# Utilities
# -----------------------------------------

# returns compatibility of sdp1 vs sdp2
def compatible(sdp1, sdp2, pessimistic=False):    
    # Returns 0 if they are incompatible
    # 1 if they may be (must be 01, 10 to be compatible)
    # 2 if they are compatible
    # 3 if they may be (must be 00, 11 to be compatible)
    # 4 if they may be (could be either)
    gametes = {}
    ambiguous = 0
    sdp1 = str(sdp1)
    sdp2 = str(sdp2)
    r = None
    for i in range(len(sdp1)):
        if sdp1[i] == '0' and sdp2[i] == '0':
            gametes['00'] = True
        elif sdp1[i] == '1' and sdp2[i] == '0':
            gametes['10'] = True
        elif sdp1[i] == '1' and sdp2[i] == '1':
            gametes['11'] = True
        elif sdp1[i] == '0' and sdp2[i] == '1':
            gametes['01'] = True
        elif sdp1[i] == '2' and sdp2[i] == '0':
            gametes['10'] = True
            gametes['00'] = True
        elif sdp1[i] == '0' and sdp2[i] == '2':
            gametes['01'] = True
            gametes['00'] = True
        elif sdp1[i] == '2' and sdp2[i] == '1':
            gametes['11'] = True
            gametes['01'] = True
        elif sdp1[i] == '1' and sdp2[i] == '2':
            gametes['11'] = True
            gametes['10'] = True
        elif sdp1[i] == '2' and sdp2[i] == '2':
            ambiguous += 1
    if len(gametes) == 4:
        r = 0
    elif ambiguous > 0:
        if len(gametes) == 3:
            if not gametes.has_key('01') or not gametes.has_key('10'):
                r = 3
            elif not gametes.has_key('00') or not gametes.has_key('11'):
                r = 1
        elif gametes.has_key('00') and gametes.has_key('11'):
            r = 3
        elif gametes.has_key('01') and gametes.has_key('10'):
            r = 1
        elif ambiguous > 1:#else:
            r = 4
    if r == None:
        r = 2
    if pessimistic:
        if r == 1 or r == 3 or r == 0:
            return 0
        if r == 4 and ambiguous > 1:
            return 0
        else:
            return 2
    else:
        return r


# returns the locations of one and two in genGroupList
def locations(one, two, genGroupList):
    oneloc = None
    twoloc = None
    for groupIndex in range(len(genGroupList)):
        group = genGroupList[groupIndex]
        if one in group[0]:
            oneloc = (groupIndex, 0)
        elif one in group[1]:
            oneloc = (groupIndex, 1)
        if two in group[0]:
            twoloc = (groupIndex, 0)
        elif two in group[1]:
            twoloc = (groupIndex, 1)
        if not oneloc == None and not twoloc == None:
            break
    return oneloc, twoloc

def removeSNP(group, snp):
    for snpgroup in group:
        for subgroup in snpgroup:
            for side in subgroup:
                if snp in side:
                    side.remove(snp)
    return group
