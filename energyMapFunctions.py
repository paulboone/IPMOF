# Functions used in Interpenetrating MOFs project
# for the calculation energy maps of MOF unit cells
# as well as preparation of input files for IPMOF javascript code
# Date: April 2016
# Author: Kutay B. Sezginel

import os
import xlrd
import math
import numpy as np


def energyMap(MOF1, atomList, cutOffArray, gridSize):
    """
    Calculate energy map for a given MOF class.
    MOF -> base (map) | atomFFparameters -> sigma and epsilon values for given atoms
    cutOff -> cut-off value for LJ potential | gridSize -> grid size array for each dimension
    Packed coordinates for MOF1 must be defined before running the function.
    """
    cutOff = min(cutOffArray)

    sortedX = sorted(MOF1.edgePoints, key=lambda x: x[0], reverse=True)
    sortedY = sorted(MOF1.edgePoints, key=lambda y: y[1], reverse=True)
    sortedZ = sorted(MOF1.edgePoints, key=lambda z: z[2], reverse=True)
    eMapMax = [math.ceil(sortedX[0][0]), math.ceil(sortedY[0][1]), math.ceil(sortedZ[0][2])]
    eMapMin = [math.floor(sortedX[-1][0]), math.floor(sortedY[-1][1]), math.floor(sortedZ[-1][2])]

    xGrid = np.linspace(eMapMin[0], eMapMax[0], (eMapMax[0]-eMapMin[0])/gridSize+1)
    yGrid = np.linspace(eMapMin[1], eMapMax[1], (eMapMax[1]-eMapMin[1])/gridSize+1)
    zGrid = np.linspace(eMapMin[2], eMapMax[2], (eMapMax[2]-eMapMin[2])/gridSize+1)

    numAtomsMOF = len(MOF1.uniqueAtomNames)
    numAtomsEnergy = len(atomList['sigma'])

    # Initialize energy map according to grid size and coordinates plus number of unique atoms
    energyMap = np.zeros([len(xGrid)*len(yGrid)*len(zGrid), numAtomsEnergy+3])

    sig, eps = LBmix(MOF1.sigma, atomList['sigma'], MOF1.epsilon, atomList['epsilon'])

    mapIndex = 0
    V = np.zeros([numAtomsEnergy])

    for x in xGrid:
        for y in yGrid:
            for z in zGrid:
                energyMap[mapIndex][0:3] = [x, y, z]
                Vtotal = np.zeros([numAtomsEnergy])
                for unitCell in MOF1.packedCoor:
                    MOF1index = 0
                    for atomCoor in unitCell:
                        atomIndex1 = MOF1.uniqueAtomNames.index(MOF1.atomName[MOF1index])
                        MOF1index += 1
                        r = coorDist(atomCoor, [x, y, z])
                        if r > cutOff:
                            continue
                        if r == 0:
                            energyMap[mapIndex][3:(numAtomsEnergy+3)] = np.ones([1, numAtomsEnergy])*math.inf
                        else:
                            for atomIndex2 in range(numAtomsEnergy):
                                V[atomIndex2] = calculateLJ(r, sig[atomIndex1][atomIndex2], eps[atomIndex1][atomIndex2])
                                Vtotal[atomIndex2] = Vtotal[atomIndex2] + V[atomIndex2]
                energyMap[mapIndex][3:(numAtomsEnergy+3)] = Vtotal
                mapIndex += 1
    return energyMap


def readMol2(MOFfile):
    """
    Reads mol2 file after opening the file in python and returns unit cell size [a, b, c],
    unit cell angles [alpha, beta, gamma], atom names and atom coordinates [x, y, z]
    """
    lineCount = 0
    readCoordinates = False

    for line in MOFfile:
        if '@<TRIPOS>ATOM' in line:
            readCoordinates = True
            atomName = []
            atomX = []
            atomY = []
            atomZ = []
            atomCoor = []
        if '@<TRIPOS>BOND' in line:
            readCoordinates = False
        if readCoordinates and '@<TRIPOS>ATOM' not in line:
            atomName.append(line.split()[1])
            atomX.append(float(line.split()[2]))
            atomY.append(float(line.split()[3]))
            atomZ.append(float(line.split()[4]))
            atomCoor.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
        if '@<TRIPOS>CRYSIN' in line:
            lineCount += 1
        if lineCount == 1 and '@<TRIPOS>CRYSIN' not in line:
            a = float(line.split()[0])
            b = float(line.split()[1])
            c = float(line.split()[2])
            alpha = float(line.split()[3])
            beta = float(line.split()[4])
            gamma = float(line.split()[5])
            lineCount += 1

    UCsize = [a, b, c]
    UCangle = [alpha, beta, gamma]
    MOFfile.close()

    return UCsize, UCangle, atomName, atomCoor


def calculateEnergyLimits(baseMOFIndex, MOFs, MOFlist, eMapAtomList, eMap):
    from tabulate import tabulate
    headers = ["MOF Name", "Atoms", "Average Energy Limit"]
    table = [[] for mof in MOFs]
    averageEnergyLimits = []

    for MOFindex2 in range(len(MOFs)):
        energyLimit = []
        atomName = []
        for atomIndex in range(len(eMapAtomList[MOFindex2]['name'])):
            sortedMap = sorted(eMap, key=lambda x: x[atomIndex+3], reverse=True)
            energyLimit.append(sortedMap[len(MOFs[baseMOFIndex].atomCoor)][atomIndex+3])
            atomName.append(eMapAtomList[MOFindex2]['name'][atomIndex])
        avgEnergyLimit = sum(energyLimit)/len(energyLimit)
        averageEnergyLimits.append(avgEnergyLimit)
        MOFname = MOFlist[MOFindex2].split('.')[0]
        table[MOFindex2].append(MOFname)
        table[MOFindex2].append(atomName)
        # table[MOFindex2].append(energyLimit)
        table[MOFindex2].append(avgEnergyLimit)

    print(tabulate(table, headers))

    return averageEnergyLimits


def exportEnergyMapjs(eMap, atomList, exportDir):
    eMapFile = open(exportDir, 'w')

    eMapFile.write("var eMapAtomNames = [];\n")
    atomIndex = 0
    for atom in atomList['name']:
        eMapFile.write("eMapAtomNames[" + str(atomIndex) + "] = '" + atom + "';\n")
        atomIndex += 1

    eMapFile.write("var eMap = [];\n")
    eMapIndex = 0
    for line in eMap:
        eMapFile.write("eMap[" + str(eMapIndex) + "] = ")
        for i in range(len(line)):
            if i == 0:
                eMapFile.write("[" + str(line[i]) + ", ")
            elif i == len(line)-1:
                eMapFile.write(str(line[i]) + "];\n")
            else:
                eMapFile.write(str(line[i]) + ", ")
        eMapIndex += 1
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


def exportMOFjs(MOF, exportDir):
    MOFfile = open(exportDir, 'w')

    MOFfile.write("var MOFatomNames = [];\n")
    atomIndex = 0
    for atom in MOF.uniqueAtomNames:
        MOFfile.write("MOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")
        atomIndex += 1

    MOFfile.write("var MOF_UCsize = [];\n")
    ucIndex = 0
    for uc in MOF.UCsize:
        MOFfile.write("MOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")
        ucIndex += 1

    MOFfile.write("var MOF_UCangle = [];\n")
    ucIndex = 0
    for uc in MOF.UCangle:
        MOFfile.write("MOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")
        ucIndex += 1

    MOFfile.write("var MOF = [];\n")
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
    atomIndex = 0
    for atom in MOF.uniqueAtomNames:
        MOFfile.write("baseMOFatomNames[" + str(atomIndex) + "] = '" + atom + "';\n")
        atomIndex += 1

    MOFfile.write("var baseMOF_UCsize = [];\n")
    ucIndex = 0
    for uc in MOF.UCsize:
        MOFfile.write("baseMOF_UCsize[" + str(ucIndex) + "] = " + str(uc) + ";\n")
        ucIndex += 1

    MOFfile.write("var baseMOF_UCangle = [];\n")
    ucIndex = 0
    for uc in MOF.UCangle:
        MOFfile.write("baseMOF_UCangle[" + str(ucIndex) + "] = " + str(uc) + ";\n")
        ucIndex += 1

    MOFfile.write("var baseMOF = [];\n")
    for MOFindex in range(len(MOF.atomName)):
        MOFfile.write("baseMOF[" + str(MOFindex) + "] = [")
        MOFfile.write(str(MOF.atomCoor[MOFindex][0]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][1]) + ", ")
        MOFfile.write(str(MOF.atomCoor[MOFindex][2]) + ", ")
        MOFfile.write("'" + str(MOF.atomName[MOFindex]) + "'" + "];\n")

    MOFfile.close()


# Export coordinates to xyz file
def exportXYZ(packedCoor, atomNames, exportDir):
    xyzFile = open(exportDir, 'w')
    xyzFile.write(str(len(packedCoor)*len(packedCoor[0])))
    xyzFile.write(exportDir.split('/')[-1])
    for UC in packedCoor:
        for atomIndex in range(len(UC)):
            line = atomNames[atomIndex] + ' '
            line += str(UC[atomIndex][0]) + ' '
            line += str(UC[atomIndex][1]) + ' '
            line += str(UC[atomIndex][2]) + ' \n'
            xyzFile.write(line)
    xyzFile.close()


# Generates a list of MOF files in given directory with given file format (ex: "".mol2")
def generateMOFlist(fileDir, fileFormat):
    fileList = os.listdir(fileDir)
    MOFlist = []
    for fileName in fileList:
        if fileFormat in fileName:
            MOFlist.append(fileName)
    return MOFlist


def getUniqueAtomList(MOFs):
    """
    Gets atom name, epsilon, and sigma values for non-repeating (unique) atoms in a list of
    MOF classes.
    """
    atomList = {'name': [], 'sigma': [], 'epsilon': []}
    for MOFindex in range(len(MOFs)):
        for atomIndex in range(len(MOFs[MOFindex].uniqueAtomNames)):
            atomList['name'].append(MOFs[MOFindex].uniqueAtomNames[atomIndex])
            atomList['sigma'].append(MOFs[MOFindex].sigma[atomIndex])
            atomList['epsilon'].append(MOFs[MOFindex].epsilon[atomIndex])

    newAtomList = {'name': [], 'sigma': [], 'epsilon': []}
    newAtomList['name'] = list(set(atomList['name']))
    newAtomList['epsilon'] = [0]*len(newAtomList['name'])
    newAtomList['sigma'] = [0]*len(newAtomList['name'])

    for atomName, sig, eps in zip(atomList['name'], atomList['sigma'], atomList['epsilon']):
        if atomName in newAtomList['name']:
            uniqueListIndex = newAtomList['name'].index(atomName)
            newAtomList['epsilon'][uniqueListIndex] = eps
            newAtomList['sigma'][uniqueListIndex] = sig

    return newAtomList


# Separates different types of atoms using the atom coordinates and atom names
# Returns unique atom coordinates and unique atoms
def separateAtoms(atomCoor, atomNames):
    uniqueAtoms = []
    for atom in atomNames:
        if atom not in uniqueAtoms:
            uniqueAtoms.append(atom)

    uniqueAtomCoor = [[] for i in range(len(uniqueAtoms))]

    for atomIndex in range(len(atomCoor)):
        uniqueAtomIndex = 0
        for uniqueAtom in uniqueAtoms:
            if atomNames[atomIndex] == uniqueAtom:
                uniqueAtomCoor[uniqueAtomIndex].append(atomCoor[atomIndex])
            uniqueAtomIndex += 1

    return uniqueAtoms, uniqueAtomCoor


# Read force field parameters from an excel file according to force field selection
def readFFparameters(excelFileDir, FFselection):
    # Read Excel File
    forceField_data = xlrd.open_workbook(excelFileDir)
    # Read columns to acquire force field parameters
    atomNames = forceField_data.sheets()[0].col_values(0)[2:]
    UFFsigma = forceField_data.sheets()[0].col_values(1)[2:]
    UFFepsilon = forceField_data.sheets()[0].col_values(2)[2:]
    DREsigma = forceField_data.sheets()[0].col_values(3)[2:]
    DREepsilon = forceField_data.sheets()[0].col_values(4)[2:]
    UFF = [atomNames, UFFsigma, UFFepsilon]
    DRE = [atomNames, DREsigma, DREepsilon]
    if FFselection == 'UFF':
        return UFF
    if FFselection == 'DRE':
        return DRE
    else:
        print('No such force field')


def getFFparameters(atomNames, FFparameters):
    """
    Get force field parameters of the atom list and force field parameters you provide
    """
    atomFFparameters = []
    for atom in atomNames:
        for i in range(len(FFparameters[0])):
            if atom == FFparameters[0][i]:
                atomName = atom
                sigma = FFparameters[1][i]
                epsilon = FFparameters[2][i]
                atomFFparameters.append([atomName, sigma, epsilon])
    return atomFFparameters


def readAtomRadius():
    import xlrd
    import numpy as np

    radiusFileDir = 'C:\\Users\\kutay\\Dropbox\\Academic\\WilmerLab\\Research\\Web-IPMOF\\IPMOFvisualization\\atomicRadius.xlsx'
    radius_data = xlrd.open_workbook(radiusFileDir)

    atomNames = radius_data.sheets()[0].col_values(0)[:]
    emprical = radius_data.sheets()[0].col_values(1)[:]
    calculated = radius_data.sheets()[0].col_values(2)[:]
    vdW = radius_data.sheets()[0].col_values(3)[:]

    # Combine radius data for calculated and emprical radius values
    combined = []
    for radius, emp in zip(calculated, emprical):
        if radius == 'no data':
            if emp != 'no data':
                radius = emp
            else:
                radius = -1
        combined.append(radius)

    # Create atom list to write name and radius values for atoms
    atomList = {'name': [], 'radius': []}
    for name, radius in zip(atomNames, combined):
        if radius != -1:
            atomList['name'].append(name)
            atomList['radius'].append(radius)

    return atomList


def gridBox(p, grid):
    """
    Calculates surrounding grid points for a given point in energy map
    Input x,y,z coordinates of point in space and grid size to get floor and ceiling values
    """
    floor = []
    ceil = []
    index = 0
    for point in p:
        floor.append(math.floor(round(point / grid[index], 1)) * grid[index])
        ceil.append(math.floor(point / grid[index]) * grid[index] + grid[index])
        index += 1
    return floor, ceil


# MOF class that holds coordinate, atom name, and unit cell information
class MOF:
    def initialize(self):
        self.UCsize, self.UCangle, self.atomName, self.atomCoor = readMol2(self.file)
        self.uniqueAtomNames, self.uniqueAtomCoor = separateAtoms(self.atomCoor, self.atomName)

    def initializeFF(self, FFtype):
        FFparam = getFFparameters(self.uniqueAtomNames, FFtype)
        self.sigma = []
        self.epsilon = []
        for i in range(len(FFparam)):
            self.uniqueAtomNames[i] = FFparam[i][0]
            self.sigma.append(FFparam[i][1])
            self.epsilon.append(FFparam[i][2])


# Calculates distance between two given coordinates in list form
def coorDist(coor1, coor2):
    return math.sqrt((coor1[0] - coor2[0])**2 + (coor1[1] - coor2[1])**2 + (coor1[2] - coor2[2])**2)


# Lorentz-Berthelot mixing rules for given lists of sigma1, sigma2, epsilon1, and epsilon2
def LBmix(sigmaList1, sigmaList2, epsilonList1, epsilonList2):
    sig = np.zeros([len(sigmaList1), len(sigmaList2)])
    eps = np.zeros([len(epsilonList1), len(epsilonList2)])

    for index1 in range(len(sigmaList1)):
        for index2 in range(len(sigmaList2)):
            sig[index1][index2] = (sigmaList1[index1] + sigmaList2[index2]) / 2
            eps[index1][index2] = math.sqrt(epsilonList1[index1] * epsilonList2[index2])

    return sig, eps


# Calculate Lennard Jones potential for give distance, sigma, and epsilon values
def calculateLJ(r, sig, eps):
    return 4 * eps * ((sig/r)**12 - (sig/r)**6)


# Packing class containing functions used for packing unit cells
class Packing:

    def factor(UCsize, cutOff):
        packingFactor = []
        for cut_off, UCdimension in zip(cutOff, UCsize):
            packingFactor.append(math.ceil(cut_off/UCdimension)*2+1)
        return packingFactor

    def vectors(packingFactor, UCsize, UCangle):

        a = UCsize[0]
        b = UCsize[1]
        c = UCsize[2]
        alpha = math.radians(UCangle[0])
        beta = math.radians(UCangle[1])
        gamma = math.radians(UCangle[2])

        xV = [a, 0, 0]
        yV = [b*math.cos(gamma), b*math.sin(gamma), 0]
        zV = [0.0]*3
        zV[0] = c*math.cos(beta)
        zV[1] = (c*b*math.cos(alpha) - yV[0]*zV[0])/yV[1]
        zV[2] = math.sqrt(c*c - zV[0]*zV[0] - zV[1]*zV[1])

        UCvectors = [xV, yV, zV]

        packingVectors = []
        for x in range(packingFactor[0]):
            for y in range(packingFactor[1]):
                for z in range(packingFactor[2]):
                    packingVectors.append([x, y, z])

        translationVectors = []
        for vector in packingVectors:
            xTranslation = xV[0] * vector[0] + yV[0] * vector[1] + zV[0] * vector[2]
            yTranslation = xV[1] * vector[0] + yV[1] * vector[1] + zV[1] * vector[2]
            zTranslation = xV[2] * vector[0] + yV[2] * vector[1] + zV[2] * vector[2]
            translationVectors.append([xTranslation, yTranslation, zTranslation])

        return translationVectors, UCvectors

    def UC(translationVectors, packingFactors, UCvectors, atomCoor):

        xV = UCvectors[0]
        yV = UCvectors[1]
        zV = UCvectors[2]

        packedCoor = [[] for i in range(len(translationVectors))]

        translationFactor = []
        originTranslationVector = []
        dim = 0
        for factor in packingFactors:
            translationFactor.append((factor-1)/2)
            originTranslationVector.append((factor-1)/2 * xV[dim] + (factor-1)/2 * yV[dim] + (factor-1)/2 * zV[dim])
            dim += 1

        packingIndex = 0
        for translation in translationVectors:
            for coor in atomCoor:
                x = coor[0] + translation[0] - originTranslationVector[0]
                y = coor[1] + translation[1] - originTranslationVector[1]
                z = coor[2] + translation[2] - originTranslationVector[2]
                packedCoor[packingIndex].append([x, y, z])
            packingIndex += 1

        return packedCoor

    def edgePoints(UCvectors):
        UCedges = []
        UCedges.append([0, 0, 0])

        for vec in UCvectors:
            UCedges.append(vec)

        # (a, b, 0)
        UCedges.append([UCvectors[0][0]+UCvectors[1][0], UCvectors[0][1]+UCvectors[1][1],
                        UCvectors[0][2]+UCvectors[1][2]])
        # (0, b, c)
        UCedges.append([UCvectors[1][0]+UCvectors[2][0], UCvectors[1][1]+UCvectors[2][1],
                        UCvectors[1][2]+UCvectors[2][2]])
        # (a, 0, c)
        UCedges.append([UCvectors[0][0]+UCvectors[2][0], UCvectors[0][1]+UCvectors[2][1],
                        UCvectors[0][2]+UCvectors[2][2]])
        # (a, b, c)
        UCedges.append([UCvectors[0][0]+UCvectors[1][0]+UCvectors[2][0],
                        UCvectors[0][1]+UCvectors[1][1]+UCvectors[2][1],
                        UCvectors[0][2]+UCvectors[1][2]+UCvectors[2][2]])
        return UCedges


def ucv(MOF):

    a = MOF.UCsize[0]
    b = MOF.UCsize[1]
    c = MOF.UCsize[2]
    alp = math.radians(MOF.UCangle[0])
    bet = math.radians(MOF.UCangle[1])
    gam = math.radians(MOF.UCangle[2])
    V = a*b*c*math.sqrt(1-math.cos(alp)**2-math.cos(bet)**2-math.cos(gam)**2+2*math.cos(alp)*math.cos(bet)*math.cos(gam))
    return V


def car2frac(coor, UCsize, UCangle, UCV):
    v = UCV

    alp = UCangle[0] / 180 * math.pi
    bet = UCangle[1] / 180 * math.pi
    gam = UCangle[2] / 180 * math.pi

    a = UCsize[0]
    b = UCsize[1]
    c = UCsize[2]

    x = coor[0]
    y = coor[1]
    z = coor[2]

    xfrac = 1/a*x
    xfrac += - math.cos(gam)/(a*math.sin(gam))*y
    xfrac += (math.cos(alp)*math.cos(gam)-math.cos(bet))/(a*v*math.sin(gam))*z

    yfrac = 0
    yfrac += 1/(b*math.sin(gam))*y
    yfrac += (math.cos(bet)*math.cos(gam)-math.cos(alp))/(b*v*sin(gam))*z

    zfrac = 0
    zfrac += 0
    zfrac += math.sin(gam)/(c*v)*z

    return [xfrac, yfrac, zfrac]


def orthorhombic(MOF):
    ortho = False
    if MOF.UCangle[0] == 90 and MOF.UCangle[1] == 90 and MOF.UCangle[2] == 90:
        ortho = True
    return ortho


def atomicVolume(MOF, radiusList):

    uniqueAtomRadius = []
    for atomName, atomRadius in zip(radiusList['name'], radiusList['radius']):
        for uniqAtom in MOF.uniqueAtomNames:
            if atomName == uniqAtom:
                uniqueAtomRadius.append(atomRadius)

    atomVolumes = 0
    for atomIndex, atom in enumerate(MOF.atomName):
        for atomName, atomRadius in zip(MOF.uniqueAtomNames, uniqueAtomRadius):
            if atomName == atom:
                atomVolumes += 4/3 * math.pi * (atomRadius/100)**3

    return atomVolumes


# Plots atom coordinates for packed unit cells
def plotPackedCell(packedCoor, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(packedCoor)):
        for coor in packedCoor[i]:
            r = i / len(packedCoor)
            ax.scatter(coor[0], coor[1], coor[2], '-o', c=[r, 0.1, 0.1], s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


# Plots atom coordinates for packed unit cells
def plotXYZ(xyzCoor, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for coor in xyzCoor:
        ax.scatter(coor[0], coor[1], coor[2], '-o', s=10)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()


def plotEnergyMap(eMap, azim, elev):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import proj3d

    def orthogonal_proj(zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, a, b],
                        [0, 0, 0, zback]])
    proj3d.persp_transformation = orthogonal_proj

    xCoor = []
    yCoor = []
    zCoor = []
    energy = []
    for i in range(len(eMap)):
        xCoor.append(eMap[i][0])
        yCoor.append(eMap[i][1])
        zCoor.append(eMap[i][2])
        energy.append(eMap[i][3])
    scale = 1E20
    colors = []
    rgb = np.zeros([len(energy), 3])
    for i in range(len(energy)):
        if energy[i] > max(energy)*0.5/scale:
            energy[i] = 1
        elif energy[i] > max(energy)*0.25/scale:
            energy[i] = 0.9
        elif energy[i] > max(energy)*0.1/scale:
            energy[i] = 0.8
        elif energy[i] > max(energy)*0.075/scale:
            energy[i] = 0.7
        elif energy[i] > max(energy)*0.05/scale:
            energy[i] = 0.6
        elif energy[i] > max(energy)*0.025/scale:
            energy[i] = 0.5
        elif energy[i] > max(energy)*0.01/scale:
            energy[i] = 0.4
        else:
            energy[i] = 0.3
        colors.append(energy[i])
        rgb[i][0] = colors[i]

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(xCoor, yCoor, zCoor, 'o', c=rgb, edgecolor='')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.azim = azim
    ax.elev = elev

    plt.show()
