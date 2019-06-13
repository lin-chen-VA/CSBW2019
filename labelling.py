#!/usr/bin/python3

import shutil
import urllib.request as request
from contextlib import closing

def xml_download(name):
    """Download xml file from PDB repository
    Args:
        name (str): residue name

    Returns:
        bool: True, download successfully; False, not able to download successfully
    """
    url = 'ftp://ftp.rcsb.org/pub/pdb/validation_reports/'+name[1:3]+'/'+name+'/'+name+'_validation.xml.gz'

    filename = name+'_validation.xml.gz'

    try:
        with closing(request.urlopen(url)) as r:
            with open(filename, 'wb') as f:
                shutil.copyfileobj(r, f)
    except Exception as err:
        print(err)
        return False

    return True

import gzip

def unpack_gz(gzfile, xmlfile):
    """Unpack xml.gz file to xml file

    Args:
        gzfile (str): xml.gz file name
        xmlfile (str): xml file name
    """
    with gzip.open(gzfile, 'rb') as f_in:
        with open(xmlfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            os.remove(gzfile)

import os
def file_exist(filename):
    """Check if a file exist

    Args:
        filename (str): file name

    Returns:
        bool: True: if the file exists; False, otherwise
    """
    exists = os.path.isfile(filename)
    if exists:
        return True
    else:
        return False

from xml.dom.minidom import parse, parseString
def get_resolution(xmlname):
    """Read the resolution from the xml file

    Args:
        xmlname (str): xml file name

    Returns:
        float: resolution
    """
    dom = parse(xmlname)
    entry = dom.getElementsByTagName("Entry")[0]
    resol = entry.getAttribute("PDB-resolution")
    return float(resol)

from os import listdir
from os.path import isfile, join
def get_files(directory):
    files = [f for f in listdir(directory) if isfile(join(directory, f))]
    files = [f for f in files if f.endswith('cif')]
    return files

import sys
import csv
def csv_dict():
    """Build a dict from csv file

    Returns:
        dict: key is the PDB ID, value is a list of residue ID, each residue ID contains chain ID, index, and residue name
    """
    csv_file = sys.argv[1]
    d = {}
    with open(csv_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if row[0] in d.keys():
                d[row[0]].append([row[1], row[2], row[3]])
            else:
                d[row[0]] = [[row[1], row[2], row[3]]]
    return d

def xml_dict(xmlname):
    """Build a dict for xml file

    Args:
        xmlname (str): xml file name

    Returns:
        dict: key is a tuple of chain ID, index, residue name, value is the label
    """
    d = {} # save xml dict
    dom = parse(xmlname)
    entries = dom.getElementsByTagName('ModelledSubgroup')
    for entry in entries:
        chain = entry.getAttribute('chain')
        if RepresentsInt(chain):
            chain = str(int(chain))
        index = entry.getAttribute('resnum')
        residueID = entry.getAttribute('resname')
        count_errors = 0
        rama = entry.getAttribute('rama')
        if rama == 'OUTLIER':
            count_errors += 1
        rota = entry.getAttribute('rota')
        if rota == 'OUTLIER':
            count_errors += 1
        clash = entry.getElementsByTagName('clash')
        if len(clash) > 0:
            count_errors += 1
        angle_outlier = entry.getElementsByTagName('angle-outlier')
        if len(angle_outlier) > 0:
            count_errors += 1
        plane_outlier = entry.getElementsByTagName('plane-outlier')
        if len(plane_outlier) > 0:
            count_errors += 1
        d[(chain, index, residueID)] = count_errors
    return d

def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

import pandas as pd
def main():
    table = pd.read_csv(sys.argv[1], header = None)

    csvdict = csv_dict()

    result_list = []
    for proteinID in csvdict:
        print('Processing '+proteinID)
        gzname = proteinID+'_validation.xml.gz'
        xmlname = proteinID+'_validation.xml'

        label_list = []

        # download xml file
        if not file_exist(xmlname):
            if xml_download(proteinID):
                unpack_gz(gzname, xmlname)
            else:
                print('-----------'+xmlname+' not exist ...')
                label_list.append('Not_exist')
                label_list = label_list*len(csvdict[proteinID])
                result_list += label_list
                continue

        # get xml dict
        xmldict = xml_dict(xmlname)

        for ele in csvdict[proteinID]:
            key = (ele[0], ele[1], ele[2]) # chain ID, index, residue name
            label = xmldict[key]
            label_list.append(label)

        result_list += label_list

    table.insert(4, 'validation', value = result_list)
    table.to_csv('output.csv', index=False, header=False)

if __name__ == '__main__':
    main()
