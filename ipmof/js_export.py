def getEnergyMapAtomList(MOFs, atomList):
    eMapAtomList = []
    for mofIndex in range(len(MOFs)):
        eMapAtomList.append({'name': [], 'epsilon': [], 'sigma': [], 'eMapAtomIndex': []})
        for atomIndex in range(len(atomList['name'])):
            for MOFatomName in MOFs[mofIndex].uniqueAtomNames:
                if atomList['name'][atomIndex] == MOFatomName:
                    eMapAtomList[mofIndex]['name'].append(atomList['name'][atomIndex])
                    eMapAtomList[mofIndex]['sigma'].append(atomList['sigma'][atomIndex])
                    eMapAtomList[mofIndex]['epsilon'].append(atomList['epsilon'][atomIndex])
                    eMapAtomList[mofIndex]['eMapAtomIndex'].append(atomIndex)

    return eMapAtomList


def exportGLOBALEnergyMapjs(eMap, avgEnergyLimit, eMapAtomList, MOFindex, exportDir):
    eMapFile = open(exportDir, 'w')

    eMapFile.write("GLOBAL.eMapAtomNames = [];\n")
    for atomIndex, atom in enumerate(eMapAtomList[MOFindex]['name']):
        eMapFile.write("eMapAtomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    eMapFile.write("GLOBAL.avgEnergyLimit = " + str(avgEnergyLimit) + ";\n")

    eMapFile.write("GLOBAL.eMapAtomIndex = [];\n")
    for atomIndex in range(len(eMapAtomList[MOFindex]['name'])):
        eMapFile.write("eMapAtomIndex[" + str(atomIndex) + "] = " + str(atomIndex+3) + ";\n")

    eMapFile.write("GLOBAL.eMap = [];\n")
    for eMapIndex, line in enumerate(eMap):
        eMapFile.write("eMap[" + str(eMapIndex) + "] = ")
        for lineIndex, eMapAtomIndex in enumerate(eMapAtomList[MOFindex]['eMapAtomIndex']):
            if lineIndex == 0:
                eMapFile.write("[" + str(line[lineIndex]) + ", ")
                eMapFile.write(str(line[lineIndex+1]) + ", ")
                eMapFile.write(str(line[lineIndex+2]) + ", ")
                eMapFile.write(str(line[eMapAtomIndex+3]) + ", ")
            elif lineIndex == len(eMapAtomList[MOFindex]['eMapAtomIndex'])-1:
                eMapFile.write(str(line[eMapAtomIndex+3]) + "];\n")
            else:
                eMapFile.write(str(line[eMapAtomIndex+3]) + ", ")
    eMapFile.close()


def exportUniqueEnergyMapjs(eMap, avgEnergyLimit, eMapAtomList, MOFindex, exportDir):
    eMapFile = open(exportDir, 'w')

    eMapFile.write("var eMapAtomNames = [];\n")
    for atomIndex, atom in enumerate(eMapAtomList[MOFindex]['name']):
        eMapFile.write("eMapAtomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    eMapFile.write("var avgEnergyLimit = " + str(avgEnergyLimit) + ";\n")

    eMapFile.write("var eMapAtomIndex = [];\n")
    for atomIndex in range(len(eMapAtomList[MOFindex]['name'])):
        eMapFile.write("eMapAtomIndex[" + str(atomIndex) + "] = " + str(atomIndex+3) + ";\n")

    eMapFile.write("var eMap = [];\n")
    for eMapIndex, line in enumerate(eMap):
        eMapFile.write("eMap[" + str(eMapIndex) + "] = ")
        for lineIndex, eMapAtomIndex in enumerate(eMapAtomList[MOFindex]['eMapAtomIndex']):
            if lineIndex == 0:
                eMapFile.write("[" + str(line[lineIndex]) + ", ")
                eMapFile.write(str(line[lineIndex+1]) + ", ")
                eMapFile.write(str(line[lineIndex+2]) + ", ")
                eMapFile.write(str(line[eMapAtomIndex+3]) + ", ")
            elif lineIndex == len(eMapAtomList[MOFindex]['eMapAtomIndex'])-1:
                eMapFile.write(str(line[eMapAtomIndex+3]) + "];\n")
            else:
                eMapFile.write(str(line[eMapAtomIndex+3]) + ", ")
    eMapFile.close()


def exportMOFjs(MOF, exportDir):
    MOFfile = open(exportDir, 'w')

    MOFfile.write("var MOFatomNames = [];\n")
    for atomIndex, atom in enumerate(MOF.uniqueAtomNames):
        MOFfile.write("MOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    MOFfile.write("var MOF_UCsize = [];\n")
    for ucIndex, uc in enumerate(MOF.UCsize):
        MOFfile.write("MOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("var MOF_UCangle = [];\n")
    for ucIndex, uc in enumerate(MOF.UCangle):
        MOFfile.write("MOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("var MOF_edgePoints = [];\n")
    for edgeIndex, edge in enumerate(MOF.edgePoints):
        MOFfile.write("MOF_edgePoints[" + str(edgeIndex) + "] = " + str(edge) + ";\n")

    MOFfile.write("var MOF = [];\n")
    for MOFindex in range(len(MOF.atomName)):
        MOFfile.write("MOF[" + str(MOFindex) + "] = [")
        MOFfile.write(str(MOF.atomCoor[MOFindex][0]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][1]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][2]) + ", ")
        MOFfile.write("'" + str(MOF.atomName[MOFindex]) + "'" + "];\n")

    MOFfile.close()


def exportGLOBALMOFjs(MOF, exportDir):
    MOFfile = open(exportDir, 'w')

    MOFfile.write("GLOBAL.MOFatomNames = [];\n")
    for atomIndex, atom in enumerate(MOF.uniqueAtomNames):
        MOFfile.write("MOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    MOFfile.write("GLOBAL.MOF_UCsize = [];\n")
    for ucIndex, uc in enumerate(MOF.UCsize):
        MOFfile.write("MOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("GLOBAL.MOF_UCangle = [];\n")
    for ucIndex, uc in enumerate(MOF.UCangle):
        MOFfile.write("MOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("GLOBAL.MOF = [];\n")
    for MOFindex in range(len(MOF.atomName)):
        MOFfile.write("MOF[" + str(MOFindex) + "] = [")
        MOFfile.write(str(MOF.atomCoor[MOFindex][0]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][1]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][2]) + ", ")
        MOFfile.write("'" + str(MOF.atomName[MOFindex]) + "'" + "];\n")

    MOFfile.close()


def exportBaseMOFjs(MOF, exportDir):
    MOFfile = open(exportDir, 'w')

    MOFfile.write("var baseMOFatomNames = [];\n")
    for atomIndex, atom in enumerate(MOF.uniqueAtomNames):
        MOFfile.write("baseMOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    MOFfile.write("var baseMOF_UCsize = [];\n")
    for ucIndex, uc in enumerate(MOF.UCsize):
        MOFfile.write("baseMOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("var baseMOF_UCangle = [];\n")
    for ucIndex, uc in enumerate(MOF.UCangle):
        MOFfile.write("baseMOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("var baseMOF_edgePoints = [];\n")
    for edgeIndex, edge in enumerate(MOF.edgePoints):
        MOFfile.write("MOF_edgePoints[" + str(edgeIndex) + "] = " + str(edge) + ";\n")

    MOFfile.write("var baseMOF = [];\n")
    for MOFindex in range(len(MOF.atomName)):
        MOFfile.write("baseMOF[" + str(MOFindex) + "] = [")
        MOFfile.write(str(MOF.atomCoor[MOFindex][0]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][1]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][2]) + ", ")
        MOFfile.write("'" + str(MOF.atomName[MOFindex]) + "'" + "];\n")

    MOFfile.close()


def exportGLOBALBaseMOFjs(MOF, exportDir):
    MOFfile = open(exportDir, 'w')

    MOFfile.write("GLOBAL.baseMOFatomNames = [];\n")
    for atomIndex, atom in enumerate(MOF.uniqueAtomNames):
        MOFfile.write("baseMOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")

    MOFfile.write("GLOBAL.baseMOF_UCsize = [];\n")
    for ucIndex, uc in enumerate(MOF.UCsize):
        MOFfile.write("baseMOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("GLOBAL.baseMOF_UCangle = [];\n")
    for ucIndex, uc in enumerate(MOF.UCangle):
        MOFfile.write("baseMOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")

    MOFfile.write("GLOBAL.baseMOF = [];\n")
    for MOFindex in range(len(MOF.atomName)):
        MOFfile.write("baseMOF[" + str(MOFindex) + "] = [")
        MOFfile.write(str(MOF.atomCoor[MOFindex][0]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][1]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][2]) + ", ")
        MOFfile.write("'" + str(MOF.atomName[MOFindex]) + "'" + "];\n")

    MOFfile.close()
