# Functions used in Interpenetrating MOFs project
# for the calculation energy maps of MOF unit cells
# as well as preparation of input files for IPMOF javascript code
# Date: April 2016
# Author: Kutay B. Sezginel

import os
import xlrd
import math
import numpy as np



def calculateEnergyLimits(baseMOFIndex, MOFs, eMapAtomList, eMap):
    """
    Calculates average energy limit and tabulates results.
    The average energy limit is calculated by sorting all the energy values in the eMap and
    taking the Nth highest energy values from the top where N is the total number of atoms in the
    baseMOF. Then the energy value is calculated as the average of all the energy values belonging
    to different types of atoms in the energy map.
    """
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
        MOFname = MOFs[MOFindex2].name
        table[MOFindex2].append(MOFname)
        table[MOFindex2].append(atomName)
        # table[MOFindex2].append(energyLimit)
        table[MOFindex2].append(avgEnergyLimit)

    print(tabulate(table, headers))

    return averageEnergyLimits


def readAtomRadius(radiusFileDir):
    import xlrd
    import numpy as np

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
                atomVolumes += 4 / 3 * math.pi * (atomRadius / 100)**3

    return atomVolumes


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
