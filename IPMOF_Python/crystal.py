# IPMOF Crystal Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import math
from forcefield import *


class MOF:
    """
    MOF class that holds coordinate, atom name, and unit cell information
    """
    def initialize(self):
        self.uc_size, self.uc_angle, self.atom_names, self.atom_coors = readMol2(self.file)
        self.uniq_atom_names_names
        self.uniq_atom_names_coors = separate_atoms(self.atom_coors, self.atom_names)

    def initializeFF(self, FFtype):
        FFparam = getFFparameters(self.uniq_atom_names_names, FFtype)
        self.sigma = []
        self.epsilon = []
        for i in range(len(FFparam)):
            self.uniq_atom_names_names[i] = FFparam[i][0]
            self.sigma.append(FFparam[i][1])
            self.epsilon.append(FFparam[i][2])


def separate_atoms(atom_coors, atom_names):
    """
    Identifies different types of atoms in given atom coordinates and atom names
    Returns unique atom coordinates and unique atom names
    """"
    uniq_atom_names = []
    for atom in atom_names:
        if atom not in uniq_atom_names:
            uniq_atom_names.append(atom)

    uniq_atom_coors = [[] for i in range(len(uniq_atom_names))]
    for atom_index in range(len(atom_coors)):
        for uniq_atom_index, uniq_atom in enumerate(uniq_atom_names):
            if atom_names[atom_index] == uniq_atom:
                uniq_atom_coors[uniq_atom_index].append(atom_coors[atom_index])

    return uniq_atom_names, uniq_atom_coors


class Packing:
    """
    Packing class containing functions used for packing unit cells
    """
    def factor(uc_size, cutOff):
        packingFactor = []
        for UCdimension in uc_size:
            packingFactor.append(math.ceil(cutOff/UCdimension)*2+1)
        return packingFactor

    def vectors(packingFactor, uc_size, uc_angle):

        a = uc_size[0]
        b = uc_size[1]
        c = uc_size[2]
        alpha = math.radians(uc_angle[0])
        beta = math.radians(uc_angle[1])
        gamma = math.radians(uc_angle[2])

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

    def UC(translationVectors, packingFactors, UCvectors, atom_coors):

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
            for coor in atom_coors:
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
    a = MOF.uc_size[0]
    b = MOF.uc_size[1]
    c = MOF.uc_size[2]
    alp = math.radians(MOF.uc_angle[0])
    bet = math.radians(MOF.uc_angle[1])
    gam = math.radians(MOF.uc_angle[2])
    V = a*b*c*math.sqrt(1-math.cos(alp)**2-math.cos(bet)**2-math.cos(gam)**2+2*math.cos(alp)*math.cos(bet)*math.cos(gam))
    return V


def car2frac(coor, uc_size, uc_angle, UCV):
    v = UCV

    alp = uc_angle[0] / 180 * math.pi
    bet = uc_angle[1] / 180 * math.pi
    gam = uc_angle[2] / 180 * math.pi

    a = uc_size[0]
    b = uc_size[1]
    c = uc_size[2]

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

def frac_pbc(frac_coor):
    pbc_coor = []
    for coor in frac_coor:
        pbc_coor.append(coor - math.floor(coor))
    return pbc_coor
