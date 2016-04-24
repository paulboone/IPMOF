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
    rotFreedom = 0
    rotLimit = 0
    energyScale = 0
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


def readJobFile(jobFileDir):
    # Read simulation time data from the job file
    jobFile = open(jobFileDir, 'r')
    for line in jobFile:
        if 'start_time' in line:
            startTime = line.split()[4]
        if 'end_time' in line:
            endTime = line.split()[4]
    jobFile.close()
    startTime = startTime.split(':')
    endTime = endTime.split(':')
    simulationTime = [0, 0, 0]
    i = 0
    for start, end in zip(startTime, endTime):
        simulationTime[i] = int(end) - int(start)
        i += 1
    textTime = str(simulationTime[0]) + ':'
    textTime += str(simulationTime[1]) + ':' + str(simulationTime[2])
    return textTime
