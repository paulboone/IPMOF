import os

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
    jobSubmission = open(jobSHdir,'w')
    for line in jobsh:
        jobSubmission.write(line + '\n')
    jobSubmission.close()
