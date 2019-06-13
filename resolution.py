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
def main():
    source_folder = '/Users/lchen/Downloads/EM_2018/'
    target_folder = "temp"
    files = get_files(source_folder)
    for f in files:
        name = f[:4]
        gzname = name+'_validation.xml.gz'
        xmlname = name+'_validation.xml'

        if not file_exist(xmlname):
            if xml_download(name): # download to the current directory
                unpack_gz(gzname, xmlname)
        else:
            print(xmlname+' exist ...')

        resolution = get_resolution(xmlname)

        low_limit = 0
        up_limit = 4

        if low_limit <= resolution < up_limit:
            shutil.move(join(source_folder, f), join(target_folder, f))

if __name__ == '__main__':
    main()
