# Functions used for reading results file of the IPMOF simulations
import os
import shutil
import matplotlib.pyplot as plt

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
                xyz.append(line[0:len(line)-2])
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

def getMinEnergyStructures(IPstructures, numStructures):
    minEnergyStructures = []
    if numStructures > len(IPstructures['structureEnergy']):
        numStructures = len(IPstructures['structureEnergy'])
        print('Omitting excess structures...')
    for strIdx in range(numStructures):
        minEnergy = sorted(IPstructures['structureEnergy'])[strIdx]
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
    totStartTime = int(startTime[0])*3600 + int(startTime[1])*60 + int(startTime[2])
    totEndTime = int(endTime[0])*3600 + int(endTime[1])*60 + int(endTime[2])
    seconds = totEndTime - totStartTime
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    textTime = "%02d:%02d:%02d" % (h, m, s)

    #return textTime
    return seconds

def exportIPxyz(IPstructures, xyzFileName, exportDir):
    xyzDir = os.path.join(exportDir, xyzFileName + '.xyz')
    xyzFile = open(xyzDir, 'w')

    # Number of atoms in the structure
    numAtoms = str(len(IPstructures['xyz']))
    xyzFile.write(numAtoms + '\n')

    # Name of the structure
    xyzFile.write(xyzFileName + '\n')

    # Atom coordinates
    for coor in IPstructures['xyz']:
        xyzFile.write(coor + '\n')

    xyzFile.close()


def exportIPpdb(IPstructures, pdbFileName, exportDir):
    pdbDir = os.path.join(exportDir, pdbFileName + '.pdb')
    pdbFile = open(pdbDir, 'w')

    pdbFile.write('HEADER    CSD ENTRY ' + pdbFileName + '\n')

    UCsize = str(IPstructures['bMOF_UCsize'][0]) + '    '
    UCsize += str(IPstructures['bMOF_UCsize'][1]) + '    '
    UCsize += str(IPstructures['bMOF_UCsize'][2])

    UCangle = str(IPstructures['bMOF_UCangle'][0]) + '  '
    UCangle += str(IPstructures['bMOF_UCangle'][1]) + '  '
    UCangle += str(IPstructures['bMOF_UCangle'][2])

    pdbFile.write('CRYST1    ' + UCsize + '  ' + UCangle + ' P1\n')

    # Atom coordinates
    for coorIndex, coor in enumerate(IPstructures['xyz']):
        coorLine = 'HETATM\t\t' + str(coorIndex) + '\t' + coor.split()[0] + '\tUNK\t1\t\t'
        coorLine += coor.split()[1] + '\t' + coor.split()[2] + '\t' + coor.split()[3]
        coorLine += '\t1.00\t0.00\t\t\t' + coor.split()[0]
        pdbFile.write(coorLine + '\n')

    pdbFile.close()


class PB:
    def source(self, sourceDir):
        self.sourceDir = sourceDir

    def input(MOFname, UCsize, UCangle, exportDir):
        PBinputDir = os.path.join(exportDir, 'input.dat')
        PBinputFile = open(PBinputDir, 'w')

        PBinputFile.write(MOFname + '\n')

        UCsizeLine = str(UCsize[0]) + ' ' + str(UCsize[1]) + ' ' + str(UCsize[2]) + '\n'
        PBinputFile.write(UCsizeLine)

        UCangleLine = str(UCangle[0]) + ' ' + str(UCangle[1]) + ' ' + str(UCangle[2]) + '\n'
        PBinputFile.write(UCangleLine)

        PBinputFile.close()

    def initialize(self, exportDir):
        sourceDir = self.sourceDir
        for inpFile in os.listdir(sourceDir):
            inpDir = os.path.join(sourceDir, inpFile)
            shutil.copy(inpDir, exportDir)

    def JobSH(jobDir, jobName, queue, wallTime):
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
        jobsh.append('cd $PBS_O_WORKDIR')
        jobsh.append('./ihome/cwilmer/kbs37/PoreBlazer/poreblazer.exe < input.dat > results.txt')
        jobsh.append('echo end_time: `date`')
        jobsh.append(' ')
        jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.OU $PBS_O_WORKDIR/$PBS_JOBID.out')
        jobsh.append('cp /var/spool/torque/spool/$PBS_JOBID.ER $PBS_O_WORKDIR/$PBS_JOBID.err')
        jobsh.append(' ')
        jobsh.append('exit')
        jobSHdir = os.path.join(jobDir, 'job.sh')
        jobSubmission = open(jobSHdir,'w')
        for line in jobsh:
            jobSubmission.write(line + '\n')
        jobSubmission.close()

    def readResults(PBdir):
        PBresultsDir = os.path.join(PBdir, 'results.txt')
        PBresultsFile = open(PBresultsDir, 'r')
        PBresultsLines = PBresultsFile.readlines()
        for line in PBresultsLines:
            if 'surface area per mass' in line:
                SA = float(line.split()[-1])
            if 'System density' in line:
                RO = float(line.split()[-1])
            if '(point accessible) volume in cm^3/g' in line:
                PV = float(line.split()[-1])
        PBres = {}
        PBres['SA'] = SA
        PBres['RO'] = RO
        PBres['PV'] = PV

        return PBres

    def PSD(resultsDir):
        psdDir = os.path.join(resultsDir, 'psd.txt')
        psdFile = open(psdDir, 'r')
        psdLines = psdFile.readlines()
        psdFile.close()

        poreSize = []
        frequency = []
        for i, line in enumerate(psdLines):
            if i > 0:
                poreSize.append(float(line.split()[0]))
                frequency.append(float(line.split()[1]))

        PSD = {}
        PSD['poreSize'] = poreSize
        PSD['frequency'] = frequency
        DPDindex = PSD['frequency'].index(max(PSD['frequency']))
        PSD['DPD'] = PSD['poreSize'][DPDindex]

        return PSD

    def plotPSD(PSD):
        import matplotlib.pyplot as plt
        plt.subplots(figsize=(8, 5))
        plt.plot(PSD['poreSize'], PSD['frequency'])
        plt.title('Pore Size Distribution', y=1.04, fontsize=18)
        plt.ylabel('Frequency', fontsize = 18)
        plt.xlabel('Pore Size (Angstrom)', fontsize = 18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.text(PSD['DPD']*1.1, max(PSD['frequency'])*0.97, 'DPD: ' + str(PSD['DPD']), fontsize=14)
        plt.show()
