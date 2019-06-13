#!/usr/bin/python

from Bio.PDB import *
from Bio.PDB.Polypeptide import *
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import *
from Bio.Seq import Seq
import numpy as np
import math
import os
import sys
import math

def getDihedral(p):
    ''' Reference: http://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python'''
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def getChains(fileName):
    ''' Get all chains in a protein

    Args:
        name (string), mmcif file name

    Return:
        chains (generator), all chains in the protein
    '''
    # Create a mmcif parser
    parser = MMCIFParser();
    # Create a structure object
    name = os.path.splitext(os.path.basename(fileName))[0];
    structure = parser.get_structure(name, fileName);
    # Get all chains
    for chain in structure.get_chains():
        yield chain;

def getResidues(chain):
    return list(chain.get_residues());

def getSCDist(atomDict, resName):
    ''' Get the distance between CA and the mass center of the side-chain

        Args:
            residue (residue)
    '''
    massCenter = getSCMassCenter(atomDict, resName);
    if len(massCenter) == 0: # no atoms
        return -1;
    return getDist(atomDict['CA'], massCenter);

def getBlockDist(atomDict, resName):
    ''' Get the distance between CA and the mass center of the side-chain

        Args:
            residue (residue)
    '''
    massCenter = getBlockMassCenter(atomDict, resName);
    if len(massCenter) == 0: # no atoms
        return -1;
    return getDist(atomDict['CA'], massCenter);

def getAtomMass(atomName):
    t = atomName[0];
    mN = 14.0067;
    mC = 12.0107;
    mO = 15.9994;
    mS = 32.065;
    if t == 'N':
        return mN;
    if t == 'C':
        return mC;
    if t == 'O':
        return mO;
    if t == 'S':
        return mS;

def getBlockMassCenter(atomDict, resName):
    ''' Get the block center

        Ref: http://formulas.tutorvista.com/physics/center-of-mass-formula.html
    '''
    l = [];
    if resName == 'ALA':
        atomList = ['CB'];
    if resName == 'ARG':
        atomList = ['NE', 'CZ', 'NH1', 'NH2'];
    if resName == 'ASN':
        atomList = ['CG', 'OD1', 'ND2'];
    if resName == 'ASP':
        atomList = ['CG', 'OD1', 'OD2'];
    if resName == 'CYS':
        atomList = ['CB', 'SG'];
    if resName == 'GLU':
        atomList = ['CD', 'OE1', 'OE2'];
    if resName == 'GLN':
        atomList = ['CD', 'OE1', 'NE2'];
    if resName == 'GLY':
        atomList = [];
    if resName == 'HIS':
        atomList = ['CG', 'ND1', 'CD2', 'CE1', 'NE2'];
    if resName == 'ILE':
        atomList = ['CD1'];
    if resName == 'LEU':
        atomList = ['CG', 'CD1', 'CD2'];
    if resName == 'LYS':
        atomList = ['CE', 'NZ'];
    if resName == 'MET':
        atomList = ['SD', 'CE'];
    if resName == 'PHE':
        atomList = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'];
    if resName == 'PRO':
        atomList = ['CB', 'CG', 'CD'];
    if resName == 'SER':
        atomList = ['OG'];
    if resName == 'THR':
        atomList = ['CB', 'OG1', 'CG2'];
    if resName == 'TRP':
        atomList = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'];
    if resName == 'TYR':
        atomList = ['OH'];
    if resName == 'VAL':
        atomList = ['CB', 'CG1', 'CG2'];
    if len(atomList) == 0:
        return l;
    massTemp = 0;
    xTemp = 0;
    yTemp = 0;
    zTemp = 0;
    for a in atomList:
        coord = atomDict[a];
        m = getAtomMass(a);
        massTemp += m;
        xTemp += coord[0]*m;
        yTemp += coord[1]*m;
        zTemp += coord[2]*m;
    return [xTemp/massTemp, yTemp/massTemp, zTemp/massTemp];

def getSCMassCenter(atomDict, resName):
    ''' Get the residue mass center

        Ref: http://formulas.tutorvista.com/physics/center-of-mass-formula.html
    '''
    l = [];
    if resName == 'ALA':
        atomList = ['CB'];
    if resName == 'ARG':
        atomList = ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'];
    if resName == 'ASN':
        atomList = ['CB', 'CG', 'OD1', 'ND2'];
    if resName == 'ASP':
        atomList = ['CB', 'CG', 'OD1', 'OD2'];
    if resName == 'CYS':
        atomList = ['CB', 'SG'];
    if resName == 'GLU':
        atomList = ['CB', 'CG', 'CD', 'OE1', 'OE2'];
    if resName == 'GLN':
        atomList = ['CB', 'CG', 'CD', 'OE1', 'NE2'];
    if resName == 'GLY':
        atomList = [];
    if resName == 'HIS':
        atomList = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'];
    if resName == 'ILE':
        atomList = ['CB', 'CG1', 'CG2', 'CD1'];
    if resName == 'LEU':
        atomList = ['CB', 'CG', 'CD1', 'CD2'];
    if resName == 'LYS':
        atomList = ['CB', 'CG', 'CD', 'CE', 'NZ'];
    if resName == 'MET':
        atomList = ['CB', 'CG', 'SD', 'CE'];
    if resName == 'PHE':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'];
    if resName == 'PRO':
        atomList = ['CB', 'CG', 'CD'];
    if resName == 'SER':
        atomList = ['CB', 'OG'];
    if resName == 'THR':
        atomList = ['CB', 'OG1', 'CG2'];
    if resName == 'TRP':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'];
    if resName == 'TYR':
        atomList = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'];
    if resName == 'VAL':
        atomList = ['CB', 'CG1', 'CG2'];
    if len(atomList) == 0:
        return l;
    massTemp = 0;
    xTemp = 0;
    yTemp = 0;
    zTemp = 0;
    for a in atomList:
        coord = atomDict[a];
        m = getAtomMass(a);
        massTemp += m;
        xTemp += coord[0]*m;
        yTemp += coord[1]*m;
        zTemp += coord[2]*m;
    return [xTemp/massTemp, yTemp/massTemp, zTemp/massTemp];

def getDist(coord_1, coord_2):
    dist = 0;
    for i in xrange(3):
        dist += (coord_1[i] - coord_2[i])*(coord_1[i] - coord_2[i]);
    return math.sqrt(dist);

def getRotamer(residues, residueList, dataSetXR, dataSetComp):
    ''' Collect the data for each specific residue

        Args:

            residues (list), residue list
            residueList (string list), list of the targeted residues
            dataSetXR(numpy.array list), XR_Block, XR_SC, XR_Phi, XR_Psi, XR_Chi1
    '''
    # open output files
    anomalyFile = open('anomaly.csv', 'a');
    #with open(name+'.csv', 'a') as outputFile:
    residueSize = len(residues);
    for index, residue in enumerate(residues):
        # Ignore the first residue and the last residue
        if index == 0 or index == residueSize-1:
            continue;
        # Calculate phi, psi and chi for the specific residues
        if residue.get_resname() in residueList:
            angles = getPhiPsiChi(residues, index);
            if len(angles) == 0:
                continue;
            # output protein id, chain id, residue index, residue name
            line = residue.get_parent().get_parent().get_parent().get_id()+','+str(residue.get_parent().get_id())+','+str(residue.get_id()[1])+','+residue.get_resname();
            # calculate dist between CA and side-chain center, dist between CA and block center
            atomDict = getResidueDict(residue);
            resName = residue.get_resname();
            SCDist = getSCDist(atomDict, resName);
            BlockDist = getBlockDist(atomDict, resName);
            line = line+','+str(SCDist);
            line = line+','+str(BlockDist);
            # output phi, psi, chi_1
            angles = angles[0:3];
            for angle in angles:
                line = line+','+str(angle);
            data = [BlockDist, SCDist] + angles;
            HBOS, COMP = getScore(data, residue.get_resname(), residueList, dataSetXR, dataSetComp);

            # get the score 
            score = ',';
            for i in HBOS:
                score += '{:10.4f},'.format(i);
            for i in COMP:
                score += '{:10.4f},'.format(i);
            HBOS = np.asarray(HBOS);
            COMP = np.asarray(COMP);
            featureList = ['Block', 'SC', 'Phi', 'Psi', 'Chi1'];
            score += featureList[np.argmax(HBOS)]+',';
            score += '{:10.4f},'.format(max(HBOS));
            score += featureList[np.argmax(COMP)]+',';
            score += '{:10.4f}'.format(max(COMP));
            line += score;
            line = line+'\n';
            anomalyFile.write(line);
            #outputFile.write(line);
    # close the output files
    anomalyFile.close();

def getScore(data, residueName, residueList, dataSetXR, dataSetComp):
    ''' get the score for a residue structure
    Args:
        data, blockDist, SCDist, phi, psi, chi1
        dataSetXR, blockXR, SCXR, phiXR, psiXR, chi1XR
        dataSetComp, blockComp, SCComp, phiComp, psiComp, chi1Comp
    Returns:
        score
    '''
    index = residueList.index(residueName)+1;
    block1 = readDistDataSet(index, data[0], dataSetXR[0]);#X-Ray, block
    SC1 = readDistDataSet(index, data[1], dataSetXR[1]);#X-Ray, SC
    phi1 = readAngleDataSet(index, data[2], dataSetXR[2]);#X-Ray, Phi
    psi1 = readAngleDataSet(index, data[3], dataSetXR[3]);#X-Ray, Psi
    chi1 = readAngleDataSet(index, data[4], dataSetXR[4]);#X-Ray, Chi1
    # get HBOS
    HBOS = [];
    # block, SC, phi, psi, chi1
    HBOS.append(processHBOSValue(block1));
    HBOS.append(processHBOSValue(SC1));
    HBOS.append(processHBOSValue(phi1));
    HBOS.append(processHBOSValue(psi1));
    HBOS.append(processHBOSValue(chi1));

    block2 = readDistDataSet(index, data[0], dataSetComp[0]);#EM-X-Ray, block
    SC2 = readDistDataSet(index, data[1], dataSetComp[1]);#EM-X-Ray, SC
    phi2 = readAngleDataSet(index, data[2], dataSetComp[2]);#EM-X-Ray, phi
    psi2 = readAngleDataSet(index, data[3], dataSetComp[3]);#EM-X-Ray, psi
    chi2 = readAngleDataSet(index, data[4], dataSetComp[4]);#EM-X-Ray, Chi2
    # get COMP
    COMP = [];
    COMP.append(processCOMPValue(block2));
    COMP.append(processCOMPValue(SC2));
    COMP.append(processCOMPValue(phi2));
    COMP.append(processCOMPValue(psi2));
    COMP.append(processCOMPValue(chi2));

    return HBOS, COMP;

def processCOMPValue(v):
    if v > 0:
        return v;
    else:
        return 0;

def processHBOSValue(v):
    if v < 0.001:
        return 5;
    else:
        return math.log(1/v, 10);

def readDistDataSet(columnIndex, value, dataSet):
    if value < 0 or value > 10:
        return 0;
    rowIndex = int(value/0.05);
    return dataSet[rowIndex, columnIndex];

def readAngleDataSet(columnIndex, value, dataSet):
    if value < 0 or value > 360:
        return 0;
    rowIndex = int(value/5);
    return dataSet[rowIndex, columnIndex];

def getResidueDict(residue):
    ''' Convert residue to a dict, key is atom name and value is coordinates '''
    rDict = {};
    atoms = list(residue.get_atom());
    for atom in atoms:
        rDict[atom.get_id()] = np.array(atom.get_coord());
    return rDict;

def getPhiPsiChi(residues, index):
    angles = [];
    if not checkResidue(residues[index-1], residues[index-1].get_resname()):
        return angles;
    prevAtoms = getResidueDict(residues[index-1]); # get last residue
    if not checkResidue(residues[index], residues[index].get_resname()):
        return angles;
    currAtoms = getResidueDict(residues[index]); # get the current residue
    if not checkResidue(residues[index+1], residues[index+1].get_resname()):
        return angles;
    nextAtoms = getResidueDict(residues[index+1]); # get the next residue
    atom_Cp = prevAtoms['C']; # C_i-1
    atom_N = currAtoms['N']; # N_i
    atom_CA = currAtoms['CA']; # CA_i
    atom_Cc = currAtoms['C']; # C_i
    atom_Nn = nextAtoms['N']; # N_i+1
    phi = getDihedral([atom_Cp, atom_N, atom_CA, atom_Cc]); # phi
    psi = getDihedral([atom_N, atom_CA, atom_Cc, atom_Nn]); # psi
    angles.append(getAngle(phi));
    angles.append(getAngle(psi));
    resName = residues[index].get_resname();
    atomList = [];
    if resName == 'ARG': # X1, X2, X3, X4
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
        atomList.append(currAtoms['NE']);
        atomList.append(currAtoms['CZ']);
    if resName == 'ASN': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'ASP': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'CYS': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['SG']);
    if resName == 'GLU': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
    if resName == 'GLN': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
    if resName == 'HIS': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'ILE': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG1']);
    if resName == 'LEU': # X1, X2
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD1']);
    if resName == 'LYS': # X1, X2, X3, X4
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['CD']);
        atomList.append(currAtoms['CE']);
        atomList.append(currAtoms['NZ']);
    if resName == 'MET': # X1, X2, X3
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
        atomList.append(currAtoms['SD']);
        atomList.append(currAtoms['CE']);
    if resName == 'PHE': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'PRO': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'SER': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['OG']);
    if resName == 'THR': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['OG1']);
    if resName == 'TRP': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'TYR': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG']);
    if resName == 'VAL': # X1
        atomList.append(currAtoms['N']);
        atomList.append(currAtoms['CA']);
        atomList.append(currAtoms['CB']);
        atomList.append(currAtoms['CG1']);
    angles = angles + getChi(atomList);
    return angles;

def getCoord(atom):
    ''' Convert atom coordinates to numpy array '''
    return np.array(atom.get_coord());

def getAngle(angle):
    ''' Convert -pi ~ pi to 0 to 2*pi '''
    if angle < 0:
        angle = 360 + angle;
    return angle;

def getChi(atomList):
    ''' Get all chi angles of the side-chain of a residue '''
    chi = [];
    if len(atomList) < 4:
        return chi;
    for i in xrange(len(atomList)-3):
        chi.append(getAngle(getDihedral([atomList[i], atomList[i+1], atomList[i+2], atomList[i+3]])));
    return chi;

def getAtom(atoms, name):
    for atom in atoms:
        if atom.get_id() == name:
            return atom;

def checkResidue(residue, name):
    ''' Validate residue
    Args:
        residue (Bio.PDB.Residue), residue
        name (string), residue name

    Return:
        True/False
    '''
    if not residue.has_id('N'):
        return False;
    if not residue.has_id('CA'):
        return False;
    if not residue.has_id('C'):
        return False;
    if not residue.has_id('O'):
        return False;
    if name == 'ALA':
        if not residue.has_id('CB'):
            return False;
    if name == 'ARG':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('NE'):
            return False;
        if not residue.has_id('CZ'):
            return False;
        if not residue.has_id('NH1'):
            return False;
        if not residue.has_id('NH2'):
            return False;
    if name == 'ASN':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('OD1'):
            return False;
        if not residue.has_id('ND2'):
            return False;
    if name == 'ASP':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('OD1'):
            return False;
        if not residue.has_id('OD2'):
            return False;
    if name == 'CYS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('SG'):
            return False;
    if name == 'GLU':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('OE1'):
            return False;
        if not residue.has_id('OE2'):
            return False;
    if name == 'GLN':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('OE1'):
            return False;
        if not residue.has_id('NE2'):
            return False;
    if name == 'HIS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('ND1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('NE2'):
            return False;
    if name == 'ILE':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
        if not residue.has_id('CD1'):
            return False;
    if name == 'LEU':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
    if name == 'LYS':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
        if not residue.has_id('CE'):
            return False;
        if not residue.has_id('NZ'):
            return False;
    if name == 'MET':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('SD'):
            return False;
        if not residue.has_id('CE'):
            return False;
    if name == 'PHE':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CZ'):
            return False;
    if name == 'PRO':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD'):
            return False;
    if name == 'SER':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('OG'):
            return False;
    if name == 'THR':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('OG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
    if name == 'TRP':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('NE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CE3'):
            return False;
        if not residue.has_id('CZ2'):
            return False;
        if not residue.has_id('CZ3'):
            return False;
        if not residue.has_id('CH2'):
            return False;
    if name == 'TYR':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG'):
            return False;
        if not residue.has_id('CD1'):
            return False;
        if not residue.has_id('CD2'):
            return False;
        if not residue.has_id('CE1'):
            return False;
        if not residue.has_id('CE2'):
            return False;
        if not residue.has_id('CZ'):
            return False;
        if not residue.has_id('OH'):
            return False;
    if name == 'VAL':
        if not residue.has_id('CB'):
            return False;
        if not residue.has_id('CG1'):
            return False;
        if not residue.has_id('CG2'):
            return False;
    return True;

def getCifFiles(directory):
    ''' Get all cif files in a specific directory
    
        Args:
            directory (string), a folder name
            
        Return:
            cif file list (generator)
        '''
    for path, subdirs, files in os.walk(directory):
            for name in files:
                fileTemp, extension = os.path.splitext(name);
                if extension == '.cif':
                    yield os.path.join(path, name)

def getResidueList():
    ''' Set a list of target residues

        Return:
            string list, residue list
    '''
    #l = ['THR'];
    l = ['ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'];
    return l;

def clearOutputs(l):
    for f in l:
        os.remove(f+'.csv');

def chain2seq(chain):
    seq = [];
    for residue in chain:
	if is_aa(residue.get_resname(), standard=True):
	    seq.append(three_to_one(residue.get_resname()))
	else:
	    seq.append("X");
    return Seq(str(''.join(seq)), IUPAC.protein)

def processFile(fileName, residueList, dataSetXR, dataSetComp):
    ''' Process a specific residue in a cif file

        Args:
            fileName (string), cif file name
            residueName (string), specific residue name
    '''
    # Get the chains in a protein
    chains = getChains(fileName);
    chains = list(chains);

    # Get the sequences for chains
    seqs = []
    for chain in chains:
        seqs.append(chain2seq(chain));

    print 'Seq: ', len(seqs);

    # Calculate identity
    seqLen = len(seqs);
    similar = set();
    for i in range(seqLen-1):
        for j in range(i+1, seqLen):
            alignments = pairwise2.align.globalxx(seqs[i], seqs[j]);
            #print alignments[0][2], alignments[0][4], alignments[0][2]/alignments[0][4];
            if alignments[0][2]/alignments[0][4] > 0.95:
                similar.add(j);

    similar = list(similar);
    #print 'Similar list: ', similar

    # Choose the chains without above 95% similarity
    c = [];
    for i, chain in enumerate(chains):
        if i not in similar:
            c.append(chain);

    chains = c;

    print '# of chains: ', len(chains);

    # Process each chain in the protein
    for chain in chains:
        residues = getResidues(chain);
        getRotamer(residues, residueList, dataSetXR, dataSetComp);
        #break; # process the first chain in a protein model

def readAnomalyData(folderName):
    print "Read X-Ray reference data ...";
    try:
        XR_Block = np.genfromtxt(os.path.join(os.path.join(folderName, 'Anomaly'), 'Block.csv'), delimiter=',');
	XR_SC = np.genfromtxt(os.path.join(os.path.join(folderName, 'Anomaly'), 'SC.csv'), delimiter=',');
	XR_Chi1 = np.genfromtxt(os.path.join(os.path.join(folderName, 'Anomaly'), 'Chi1.csv'), delimiter=',');
	XR_Phi = np.genfromtxt(os.path.join(os.path.join(folderName, 'Anomaly'), 'Phi.csv'), delimiter=',');
	XR_Psi = np.genfromtxt(os.path.join(os.path.join(folderName, 'Anomaly'), 'Psi.csv'), delimiter=',');
    except Exception:
	print 'Not able to read X-Ray reference data ...'
	sys.exit()
    return (XR_Block, XR_SC, XR_Phi, XR_Psi, XR_Chi1);

def readOverEstimate(folderName):
    print "Read the histogram comparison data ...";
    try:
        Compare_Block = np.genfromtxt(os.path.join(os.path.join(folderName, 'OverEstimate'), 'Block.csv'), delimiter=',');
	Compare_SC = np.genfromtxt(os.path.join(os.path.join(folderName, 'OverEstimate'), 'SC.csv'), delimiter=',');
	Compare_Chi1 = np.genfromtxt(os.path.join(os.path.join(folderName, 'OverEstimate'), 'Chi1.csv'), delimiter=',');
	Compare_Phi = np.genfromtxt(os.path.join(os.path.join(folderName, 'OverEstimate'), 'Phi.csv'), delimiter=',');
	Compare_Psi = np.genfromtxt(os.path.join(os.path.join(folderName, 'OverEstimate'), 'Psi.csv'), delimiter=',');
    except Exception:
	print 'Not able to read the histogram comparison data ...'
	sys.exit()
    return (Compare_Block, Compare_SC, Compare_Phi, Compare_Psi, Compare_Chi1);

def main(l):
    ''' Generate dihedral angles for a specific residue 
    
        Args:
            residueName (string), target residue
    '''
    # read the anomaly data
    XR_Block, XR_SC, XR_Phi, XR_Psi, XR_Chi1 = readAnomalyData(sys.argv[2]);
    dataSetXR = [XR_Block, XR_SC, XR_Phi, XR_Psi, XR_Chi1];
    print XR_Block.shape, XR_SC.shape, XR_Phi.shape, XR_Psi.shape, XR_Chi1.shape
    compare_Block, compare_SC, compare_Phi, compare_Psi, compare_Chi1 = readOverEstimate(sys.argv[2]);
    dataSetComp = [compare_Block, compare_SC, compare_Phi, compare_Psi, compare_Chi1];
    print compare_Block.shape, compare_SC.shape, compare_Phi.shape, compare_Psi.shape, compare_Chi1.shape

    # Get all cif files in a directory
    cifFiles = getCifFiles(sys.argv[1]);
    count = 0;

    for cifFile in cifFiles:
        count += 1;
        print str(count)+'th: '+cifFile;
        try:
            processFile(cifFile, l, dataSetXR, dataSetComp);
        except Exception, e:
            print str(e);

if __name__ == '__main__':
    ''' Generate the score for each residue
    	python -W ignore anomaly.py cifFolder Analysis
    '''
    l = getResidueList();# get the target residue list
    main(l); # check the targeted residues only

