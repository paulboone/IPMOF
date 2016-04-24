import os
import shutil

# IPMOF.js source file directory
source_dir = "C:\\Kutay\\IPMOF\\IPMOF_Node"
# MOF structure files and eMap files directory
inputDir = "C:\\Kutay\\IPMOF_Input"
# Simulation files export directory
export_dir = "C:\\Kutay\\IPMOF_Frank"

os.chdir(source_dir)

from exportFrank import *
# --------------------------------------------------------------------------------------------------
# Energy limit scaling, rotationFreedom, rotationAmount, outputStructures
ES = [1E-2]
RF = [30]
RL = [50]
XYZ = ['true']
queue = 'shared'
wallTime = '01:00:00'
# --------------------------------------------------------------------------------------------------
# Generate lists for simulation parameters and file names to create for the simulations
exportParamList, exportFolderList = generateExportList(ES, RF, RL, XYZ)

# Directory for IPMOF_Functions source file to copy in each directory
IPfunDir = os.path.join(source_dir, 'IPMOF_Functions.js')

for baseMOF in os.listdir(inputDir):

    mobileMOFdir = os.path.join(inputDir, baseMOF)
    baseMOF_exportDir = os.path.join(export_dir, baseMOF)
    os.mkdir(baseMOF_exportDir)

    for mobileMOF in os.listdir(mobileMOFdir):

        inputFolderDir = os.path.join(mobileMOFdir, mobileMOF)
        mobileMOF_exportDir = os.path.join(baseMOF_exportDir, mobileMOF)
        os.mkdir(mobileMOF_exportDir)

        inputFiles = os.listdir(inputFolderDir)

        # Read input files created by initialize_IPMOFjs (MOF, baseMOF, eMap)
        IPinpFiles = []
        for inp in inputFiles:
            IPinpFiles.append(inp)
        # Add IPMOF_Functions to the file list (this list is used in IPMOF to require these files)
        IPinpFiles.append('IPMOF_Functions.js')

        # For each simulation parameter set and each MOF combination
        # Copy all the required files (MOF, baseMOF, eMap, IPMOF_Functions)
        # Export IPMOF with given parameters and export job.sh with given ob submission parameters
        jobNum = 0
        for expPar, expFold in zip(exportParamList, exportFolderList):
            expDir = os.path.join(mobileMOF_exportDir, expFold)
            os.mkdir(expDir)

            for inpFile in inputFiles:
                dstDir = os.path.join(expDir, inpFile)
                inpDir = os.path.join(inputFolderDir, inpFile)
                shutil.copy(inpDir, dstDir)

            IPfun_dstDir = os.path.join(expDir, 'IPMOF_Functions.js')

            shutil.copy(IPfunDir, IPfun_dstDir)

            exportIPMOFjs(IPinpFiles, expPar, source_dir, expDir)

            jobName = baseMOF + '-' + mobileMOF + str(jobNum)
            exportJobSH(expDir, jobName, queue, wallTime)
            jobNum += 1
