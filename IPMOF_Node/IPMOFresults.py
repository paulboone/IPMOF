# Functions used for reading results file of the IPMOF simulations
import os

def readStructures(resultsDir):
    resultsFileName = 'results.txt'
    resDir = os.path.join(resultsDir, resultsFileName)
    resFile = open(resDir, 'r')
    IPstructures = {'numAtoms': [], 'structureIndex': [], 'IPtrial': [], 'rotation': [], 'structureEnergy': [], 'xyz': []}
    newStructure = True
    lineCount = 0
    xyz = []
    for line in resFile:

        if 'Initial' in line:
            break

        if '---' in line:
            IPstructures['numAtoms'].append(numAtoms)
            IPstructures['structureIndex'].append(structureIndex)
            IPstructures['IPtrial'].append(IPtrial)
            IPstructures['rotation'].append(rotation)
            IPstructures['structureEnergy'].append(structureEnergy)
            IPstructures['xyz'].append(xyz)
            newStructure = True
            lineCount = 0
            xyz = []
            continue

        if newStructure:
            if lineCount == 0:
                numAtoms = float(line)
            if lineCount == 1:
                structureIndex = float(line.split()[1])
                IPtrial = float(line.split()[4])
            if lineCount == 2:
                rotation = [float(line.split()[2]), float(line.split()[4]), float(line.split()[6])]
            if lineCount == 3:
                structureEnergy = float(line.split()[1])
            if lineCount > 4:
                xyz.append(line[0:len(line)-2])
            lineCount += 1

    return IPstructures


def getMinEnergyStructures(IPstructures, numStructures):
    minEnergyStructures = []
    for strIdx in range(numStructures):
        minEnergy = sorted(IPstructures['structureEnergy'])[0]
        minIndex = IPstructures['structureEnergy'].index(minEnergy)
        minEnergyStructures.append(minIndex)
    return minEnergyStructures


def readSummary(resultsDir):
    resultsFileName = 'results.txt'
    resDir = os.path.join(resultsDir, resultsFileName)
    resFile = open(resDir, 'r')
    percent = []
    structureCount = []
    trialCount = []
    icCount = []
    for line in resFile:

        if 'Rotational Freedom' in line:
            rotFreedom = float(line.split()[2])

        if 'Rotation Limit' in line:
            rotLimit = float(line.split()[2])

        if 'Energy Scale' in line:
            energyScale = float(line.split()[2])

        if 'Percent' in line:
            percent.append(float(line.split()[1]))
            structureCount.append(float(line.split()[4]))
            trialCount.append(float(line.split()[7]))
            icCount.append(float(line.split()[10]))

    summary = {}
    summary['rotFreedom'] = rotFreedom
    summary['rotLimit'] = rotLimit
    summary['energyScale'] = energyScale
    summary['percent'] = percent
    summary['structureCount'] = structureCount
    summary['trialCount'] = trialCount
    summary['icCount'] = icCount

    return summary
