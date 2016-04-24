import os
from IPMOFresults import *

# Folder containing results file
outputDir = "C:\\Kutay\\IPMOF_Output"

baseMOF = os.listdir(outputDir)
table = []
totalStructureCount = 0
completed = '---'
simCount = 0

for bMOF in baseMOF:
    mobileMOF = os.listdir(os.path.join(outputDir,bMOF))

    for mMOF in mobileMOF:
        simPar = os.listdir(os.path.join(outputDir,bMOF, mMOF))

        for sim in simPar:
            simDir = os.path.join(outputDir, bMOF, mMOF, sim)
            resultsDir = os.path.join(outputDir,bMOF, mMOF, sim, 'results.txt')

            # Check if job .out file exists
            for simFile in os.listdir(simDir):
                if '.out' in simFile:
                    jobFileDir = os.path.join(outputDir,bMOF, mMOF, sim, simFile)
                    runTime = readJobFile(jobFileDir)
                    #runTime = '01:00'
                    break
                else:
                    runTime = '---'

            # Check if results file exists
            if os.path.isfile(resultsDir):
                resultsSize = os.path.getsize(resultsDir)

                if resultsSize > 0:
                    summary = readSummary(simDir)
                    IPstructures = readStructures(simDir)
                    #minEstructures = getMinEnergyStructures(IPstructures, 5)
                    if len(summary['percent']) > 0:
                        completed = summary['percent'][-1]
                        totalStructureCount = summary['structureCount'][-1]
                    else:
                        completed = 'Interrupted'
                        totalStructureCount = IPstructures['structureIndex'][-1]
                else:
                    completed = 'Error'
                    totalStructureCount = 'Err'

            else:
                completed = 'NoResults'
                totalStructureCount = 'NoResults'

        table.append([])
        table[simCount].append(bMOF)
        table[simCount].append(mMOF)
        table[simCount].append(sim)
        table[simCount].append(completed)
        table[simCount].append(runTime)
        table[simCount].append(totalStructureCount)
        simCount += 1

from tabulate import tabulate

headers = ["Base MOF", "Mobile MOF", "SimPar", "Completed", "RunTime", "sCount"]
print(tabulate(table, headers))
