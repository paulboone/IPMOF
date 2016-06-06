# IPMOF nodeJS python library
import os


# --------------------------------------------------------------------------------------------------
# Functions for reading nodeJS output files
# --------------------------------------------------------------------------------------------------
def readStructures(resultsDir):
    resultsFileName = 'results.txt'
    resDir = os.path.join(resultsDir, resultsFileName)
    resFile = open(resDir, 'r')
    resLines = resFile.readlines()
    resFile.close()

    # Count the total number of structures found by looking at dashed lines
    totalStructureCount = 0
    for line in resLines:
        if 'Base' in line:
            break

        if '---' in line:
            totalStructureCount += 1
            lineIndex = resLines.index(line) + 1
            nextLine = resLines[lineIndex]
            if '---' in nextLine:
                totalStructureCount -= 1
            if 'Base' in nextLine:
                totalStructureCount -= 1

    # Create a dictionary for storing information about the structures discovered in simulations
    IPstructures = {'numAtoms': [], 'structureIndex': [], 'IPtrial': [], 'rotation': [], 'structureEnergy': [], 'xyz': []}
    newStructure = True
    structureCount = 0
    lineCount = 0
    xyz = []
    for line in resLines:

        if structureCount == totalStructureCount:
            break

        if '---' in line:
            IPstructures['numAtoms'].append(numAtoms)
            IPstructures['structureIndex'].append(structureIndex)
            IPstructures['IPtrial'].append(IPtrial)
            IPstructures['rotation'].append(rotation)
            IPstructures['structureEnergy'].append(structureEnergy)
            IPstructures['xyz'].append(xyz)
            newStructure = True
            structureCount += 1
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
                xyz.append(line[0:len(line) - 2])
            lineCount += 1

    return IPstructures


def readUnitCell(IPstructures, resultsDir):
    resultsFileName = 'results.txt'
    resDir = os.path.join(resultsDir, resultsFileName)
    resFile = open(resDir, 'r')
    resLines = resFile.readlines()
    resFile.close()

    for line in resLines:
        if 'Base MOF Unit Cell a:' in line:
            bMOF_a = float(line.split()[5])
            bMOF_b = float(line.split()[7])
            bMOF_c = float(line.split()[9])
        if 'Base MOF Unit Cell alpha:' in line:
            bMOF_alpha = float(line.split()[5])
            bMOF_beta = float(line.split()[7])
            bMOF_gamma = float(line.split()[9])
        if 'Mobile MOF Unit Cell a:' in line:
            mMOF_a = float(line.split()[5])
            mMOF_b = float(line.split()[7])
            mMOF_c = float(line.split()[9])
        if 'Mobile MOF Unit Cell alpha:' in line:
            mMOF_alpha = float(line.split()[5])
            mMOF_beta = float(line.split()[7])
            mMOF_gamma = float(line.split()[9])

    IPstructures['bMOF_UCsize'] = [bMOF_a, bMOF_b, bMOF_c]
    IPstructures['bMOF_UCangle'] = [bMOF_alpha, bMOF_beta, bMOF_gamma]
    IPstructures['mMOF_UCsize'] = [mMOF_a, mMOF_b, mMOF_c]
    IPstructures['mMOF_UCangle'] = [mMOF_alpha, mMOF_beta, mMOF_gamma]

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


# --------------------------------------------------------------------------------------------------
# Functions for exporting javascript files
# --------------------------------------------------------------------------------------------------
def exportIPMOFjs(fileNames, exportParameters, source_dir, export_dir):
    IPMOFdir = os.path.join(source_dir, 'IPMOF.js')
    newIPMOFdir = os.path.join(export_dir, 'IPMOF.js')

    energyScale = str(exportParameters[0])
    rotationFreedom = str(exportParameters[1])
    rotationAmount = str(exportParameters[2])
    outputStructures = exportParameters[3]

    IPMOF = open(IPMOFdir, 'r')
    newIPMOF = open(newIPMOFdir, 'w')

    for file in fileNames:
        newIPMOF.write('require("./' + file + '");\n')

    for line in IPMOF:
        if 'var energyScale' in line:
            line = 'var energyScale = ' + energyScale + ';\n'
        if 'var rotationFreedom' in line:
            line = 'var rotationFreedom = ' + rotationFreedom + ';\n'
        if 'var rotationLimit' in line:
            line = 'var rotationLimit = ' + rotationAmount + ';\n'
        if 'var outputStructures' in line:
            line = 'var outputStructures = ' + outputStructures + ';\n'
        newIPMOF.write(line)
    IPMOF.close()
    newIPMOF.close()


def generateExportList(energyScale, rotationFreedom, rotationLimit, outputStructure):
    totalParam = len(energyScale) * len(rotationFreedom) * len(rotationLimit) * len(outputStructure)
    exportParamList = [[] for i in range(totalParam)]
    exportFolderList = []
    parIdx = 0
    for ES in energyScale:
        for RF in rotationFreedom:
            for RL in rotationLimit:
                for XYZ in outputStructure:
                    es = "{:.0E}".format(ES)
                    if '+' in es:
                        es = '+' + es.split('+')[-1]
                    if '-' in es:
                        es = '-' + es.split('-')[-1]
                    if '0' in es:
                        es = es.replace('0', '')
                    if XYZ == 'true':
                        exportFolderList.append('E' + es + '_F' + str(RF) + '_L' + str(RL) + '_S')
                    else:
                        exportFolderList.append('E' + es + '_F' + str(RF) + '_L' + str(RL))
                    exportParamList[parIdx].append(ES)
                    exportParamList[parIdx].append(RF)
                    exportParamList[parIdx].append(RL)
                    exportParamList[parIdx].append(XYZ)
                    parIdx += 1

    return exportParamList, exportFolderList


def exportJobSH(jobDir, jobName, queue, wallTime):
    jobsh = []
    jobsh.append('#!/bin/bash')
    jobsh.append(' ')
    jobsh.append('#PBS -j oe')
    jobsh.append('#PBS -N ' + jobName)
    jobsh.append('#PBS -q ' + queue)
    jobsh.append('#PBS -l nodes=1:ppn=1')
    jobsh.append('#PBS -l walltime=' + wallTime)
    jobsh.append('#PBS -S /bin/bash')
    jobsh.append(' ')
    jobsh.append('echo JOB_ID: $PBS_JOBID JOB_NAME: $PBS_JOBNAME HOSTNAME: $PBS_O_HOST')
    jobsh.append('echo "start_time: `date`"')
    jobsh.append(' ')
    jobsh.append('module use --append /ihome/cwilmer/pab135/modules/modulefiles')
    jobsh.append('module load nodejs/5.7.0')
    jobsh.append(' ')
    jobsh.append('cd $PBS_O_WORKDIR')
    jobsh.append('node IPMOF.js > results.txt')
    jobsh.append('echo end_time: `date`')
    jobsh.append(' ')
    jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.OU $PBS_O_WORKDIR/$PBS_JOBID.out')
    jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.ER $PBS_O_WORKDIR/$PBS_JOBID.err')
    jobsh.append(' ')
    jobsh.append('exit')
    jobSHdir = os.path.join(jobDir, 'job.sh')
    jobSubmission = open(jobSHdir, 'w')
    for line in jobsh:
        jobSubmission.write(line + '\n')
    jobSubmission.close()
