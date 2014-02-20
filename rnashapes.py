"""
This is a homemade version of the Shape extractor function of RNAshapes
http://bibiserv2.cebitec.uni-bielefeld.de/rnashapes
"""

import argparse
import math

def checkStructure(structure):
    valid_chars = True
    for char in structure:
        if char not in ["(", ")", "."]:
            valid_chars = False
            break

    result = (structure.count("(") == structure.count(")")) and valid_chars
    if valid_chars == False:
        print structure
    return result

def find_shapes(structure):
    #assert given structure is adequate
    assert (checkStructure(structure))
    list_opener = []
    list_stems = []

    i = 0
    list_stem_end = []
    # separate in stems
    while i < len(structure):
        if structure[i] == "(":
            list_opener.append(i)
        elif structure[i] == ")":
            current_stem = dict(opener=[], closer=[], open_dict=dict())
            while i < len(structure):

                if structure[i] == ")":
                    opener = list_opener.pop(-1)
                    current_stem["opener"].append(opener)
                    current_stem["closer"].append(i)
                    current_stem["open_dict"][str(opener)] = i

                    if list_opener and list_opener[-1] in list_stem_end:
                        list_stem_end.append(list_opener[-1])
                        break

                elif structure[i] == "(":

                    if list_opener:
                        list_stem_end.append(list_opener[-1])
                    i -= 1
                    break
                i += 1
            list_stems.append(current_stem)
        i += 1

    # build the level1 for each stems
    range_occupied = []
    dict_lvl1 = dict()
    for stem in list_stems:
        range_open = range(min(stem['opener']), max(stem['opener'])+1)
        range_close = range(min(stem['closer']), max(stem['closer'])+1)
        range_occupied.extend(range_open+range_close)

        temp_lvl1_open = ""
        temp_lvl1_close = ""

        last_opener = None
        last_closer = None
        for opener in sorted(stem['opener']):
            if last_opener == None:
                temp_lvl1_open += "["
                temp_lvl1_close = "]" + temp_lvl1_close
            else:
                if math.fabs(opener - last_opener) != 1:
                    temp_lvl1_open += "_"
                if math.fabs(stem["open_dict"][str(opener)] - last_closer) != 1:
                    temp_lvl1_close = "_" + temp_lvl1_close
                if temp_lvl1_open.endswith("_") or temp_lvl1_close.startswith("_"):
                    temp_lvl1_open += "["
                    temp_lvl1_close = "]" + temp_lvl1_close
            last_opener = opener
            last_closer = stem["open_dict"][str(opener)]

        dict_lvl1[str(min(stem['opener']))] = dict(elem=temp_lvl1_open, lvl5="[")
        dict_lvl1[str(min(stem['closer']))] = dict(elem=temp_lvl1_close, lvl5="]")

    # assemble level1
    level1 = ""
    level5 = ""
    for i, element in enumerate(structure):
        if str(i) in dict_lvl1:
            level1 += dict_lvl1[str(i)]["elem"].strip()
            level5 += dict_lvl1[str(i)]["lvl5"]
        if element == "." and not level1.endswith("_") and not i in range_occupied:
            level1 += "_"

    level1 = level1.strip().replace("[_]", "[]").replace(" ", "")
    level3 = level1.replace("_", "")

    return level5, level3, level1

def benchmark(FILE):
    #file is supposed to be organized in groupe of 4 rows
    #dot bracket
    #shape level 5
    #shape level 3
    #shape level 1
    ## next group
    file = open(FILE, "r")
    data = file.readlines()
    file.close()
    #clean the input
    data2 =[]
    for i in data:
        data2.append(i.strip())
    #print data2
    
    failed5 = []
    failed3 = []
    failed1 = []
    numFailed5 =0
    numFailed3 =0
    numFailed1 =0
    success = 0
    for i in range(0, len(data2)):
        if i%4 == 0:
            x = find_shapes(data2[i])
            if ((x[0] != data2[i+1])):
                conflict = []
                conflict.append(data2[i]) #dotbracket
                conflict.append(x[0]) #Python result
                conflict.append(data2[i+1]) #RNAshapes result
                failed5.append(conflict)
                numFailed5+=1
            elif(x[1] != data2[i+2]):
                conflict = []
                conflict.append(data2[i])
                conflict.append(x[1])
                conflict.append(data2[i+2])
                failed3.append(conflict)
                numFailed3+=1
            elif(x[2] != data2[i+3]):
                conflict = []
                conflict.append(data2[i])
                conflict.append(x[2])
                conflict.append(data2[i+3])
                failed1.append(conflict)
                numFailed1+=1
            else:
                success += 1
                
    mistakesCases = numFailed1 + numFailed3 + numFailed5
    numCases = len(data2)/4*3
    successCases = float(numCases-mistakesCases)/numCases
    return [failed5, failed3, failed1, numFailed5, numFailed3, numFailed1, successCases]
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--structure',action="store", required=True, dest="structure")

    ns = parser.parse_args()

    structure = ns.structure

    level5, level3, level1 = find_shapes(structure)
    print level1
    print level3
    print level5