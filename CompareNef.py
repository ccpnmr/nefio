"""
compareNef - a series of routines to compare the contents of two Nef files

Command Line Usage:
  compareNef for execution from the command line with a suitable script
  An example can be found in AnalysisV3/bin/compareNef:

      #!/usr/bin/env sh
      export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
      export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
      export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
      ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/compareNef.py $*

  Usage:  compareNef inFile1, inFile2     compare two Nef files
                                          print results to the screen
          compareNef /?                   simple instructions

Details of the contents of Nef files can be found in GenericStarParser
The general structure of a Nef file is:

::

    DataExtent
      DataBlock
        Item
        Loop
        SaveFrame
          Item
          Loop

DataExtent, DataBlock and SaveFrame are Python OrderedDict with an additional 'name' attribute
DataBlocks and SaveFrames are entered in their container using their name as the key.

Loop is an object with a 'columns' list, a 'data' list-of-row-OrderedDict, and a name attribute
set equal to the name of the first column. A loop is entered in its container under each
column name, so that e.g. aSaveFrame['_Loopx.loopcol1'] and aSaveFrame['_Loopx.loopcol2'] both
exist and both correspond to the same loop object.

Items are entered as a string key - string value pair.


Module Contents
===============

compareNef contains the following routines:
  compareNefFiles      compare two Nef files and return a compare object of the form:

  ::

      [] list of NefItems of the form:

          inWhich     a flag labelling which file the item was found in
                      1 = found in the first file, 2 = found on the second file, 3 = common to both
          []          list of strings containing the comparison information

  compareDataExtents    compare two DataExtent objects and return a compare list as above
  compareDataBlocks     compare two DataBlock objects and return a compare list as above
  compareSaveFrames     compare two SaveFrame objects and return a compare list as above
  compareLoops          compare two Loop objects and return a compare list as above

"""

#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-04-07 11:40:47 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
# import time
import copy
import sys
from ccpn.util.nef import GenericStarParser, StarIo
from ccpn.util import Path
# import unittest
# import contextlib
# import difflib
# import json
# from collections import OrderedDict

TEST_FILE_PATH = os.path.join(Path.getTopDirectory(), 'internal', 'data', 'starExamples')

DATAEXTENT = ''
DATABLOCK = ''
SAVEFRAME = ''
LOOP = ''
COLUMN = ''

# DATAEXTENT = 'dataExtent:'
# DATABLOCK = 'dataBlock:'
# SAVEFRAME = 'saveFrame:'
# LOOP = 'Loop:'
# COLUMN = 'Column'

#=========================================================================================
# nefItem
#=========================================================================================

class nefItem():
  """
  Holds the contents of a single Nef comparison
  inWhich   a flag labelling which file the item was found in
            1 = found in the first file, 2 = found on the second file, 3 = common to both
  list      a list of strings containing the comparison information
  """
  def __init__(self, cItem=None):
    self.inWhich = None
    self.list = []

    if cItem is not None:
      self.list = copy.deepcopy(cItem.list)
      self.inWhich = cItem.inWhich

#=========================================================================================
# _loadGeneralFile
#=========================================================================================

def _loadGeneralFile(path=None):
  """
  Load a file with the given pathname and return a dict of the contents
  :return entry:dict
  """
  usePath = path if path.startswith('/') else os.path.join(TEST_FILE_PATH, path)
  # t0 = time.time()
  entry = StarIo.parseNefFile(usePath)  # 'lenient')
  # print("Parsing time %s for %s" % (time.time() - t0, path))
  print (' %s' % path)
  return entry

#=========================================================================================
# printFile
#=========================================================================================

def printFile(thisFile):
  """
  Print a file to the screen
  """
  print ('~'*80)
  print(thisFile)
  for i, val in enumerate(thisFile):
    print (i, thisFile[val])
    sub = thisFile[val]
    for j, valj in enumerate(sub):
      print ('  ', j, sub[valj])
      if j > 3:
        break

      sub2 = sub[valj]
      for k, valk in enumerate(sub2):
        loopType = sub2[valk]
        if isinstance(loopType, GenericStarParser.Loop):
          print ('    ', k, 'LOOP', loopType)
        else:
          print('    ', k, loopType)

        if k > 3:
          break

#=========================================================================================
# addToList
#=========================================================================================

def addToList(inList, cItem=nefItem(), nefList=[]) -> list:
  """
  Append cItem to the compare list
  :param inList: a list of items to add to the compare list
  :return: list of type nefItem
  """
  if len(inList) > 0:
    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append('contains --> '+', '.join(inList))
    nefList.append(nefItem(cItem=cItem3))

  return nefList

#=========================================================================================
# compareLoops
#=========================================================================================

def compareLoops(loop1:GenericStarParser.Loop
                , loop2:GenericStarParser.Loop
                , cItem=nefItem()
                , nefList=[]) -> list:
  """
  Compare two Loops
  :param loop1: name of the first Loop object
  :param loop2: name of the second Loop object
  :return: list of type nefItem
  """
  lSet = [bl for bl in loop1.columns]
  rSet = [bl for bl in loop2.columns]
  inLeft = set(lSet).difference(rSet)
  dSet = set(lSet).intersection(rSet)
  inRight = set(rSet).difference(lSet)

  cItem1 = copy.deepcopy(cItem)
  cItem1.list.append(LOOP+loop1.name)
  cItem1.inWhich = 1
  addToList(inLeft, cItem=cItem1, nefList=nefList)

  cItem2 = copy.deepcopy(cItem)
  cItem2.list.append(LOOP+loop2.name)
  cItem2.inWhich = 2
  addToList(inRight, cItem=cItem2, nefList=nefList)

  if loop1.data and loop2.data:
    rowRange = min(len(loop1.data), len(loop2.data))
    if len(loop1.data) == len(loop2.data):        # simple compare, same length tables
      for compName in dSet:
        for rowIndex in range(rowRange):
          if loop1.data[rowIndex][compName] != loop2.data[rowIndex][compName]:
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(LOOP+loop1.name)
            cItem3.list.append('<Column>: '+compName+'  <rowIndex>: '\
                        +str(rowIndex)+'  -->  '\
                        +str(loop1.data[rowIndex][compName])+' != '\
                        +str(loop2.data[rowIndex][compName]))
            cItem3.inWhich = 3
            nefList.append(nefItem(cItem=cItem3))
    else:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop1.name)
      cItem3.list.append('<rowLength>:  '+str(len(loop1.data))+' != '+str(len(loop2.data)))
      cItem3.inWhich = 3
      nefList.append(nefItem(cItem=cItem3))

      #TODO
      # need to add a further test here, could do a diff on the tables which would pick up
      # insertions to the table - this columns would need to be reordered for this to work
      # what if there are a different number of columns?
      # also check for Mandatory items

  else:
    # can't compare non-existant loopdata
    if loop1.data is None:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop1.name)
      cItem3.list.append('<Contains no data>')
      cItem3.inWhich = 1
      nefList.append(nefItem(cItem=cItem3))
    if loop2.data is None:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop2.name)
      cItem3.list.append('<Contains no data>')
      cItem3.inWhich = 2
      nefList.append(nefItem(cItem=cItem3))

  return nefList

#=========================================================================================
# compareSaveFrames
#=========================================================================================

def compareSaveFrames(saveFrame1:GenericStarParser.SaveFrame
                      , saveFrame2:GenericStarParser.SaveFrame
                      , cItem=nefItem()
                      , nefList=[]) -> list:
  """
  Compare two saveFrames, if they have the same name then check their contents
  :param saveFrame1: name of the first SaveFrame object
  :param saveFrame2: name of the second SaveFrame object
  :return: list of type nefItem
  """
  lSet = [' ' if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else saveFrame1[bl].name for bl in saveFrame1]
  rSet = [' ' if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else saveFrame2[bl].name for bl in saveFrame2]
  inLeft = set(lSet).difference(rSet).difference({' '})
  dSet = set(lSet).intersection(rSet).difference({' '})
  inRight = set(rSet).difference(lSet).difference({' '})

  # lVSet = [str(bl)+':'+str(saveFrame1[bl]) if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame1]
  # rVSet = [str(bl)+':'+str(saveFrame2[bl]) if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame2]

  lVSet = [str(bl) if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame1]
  rVSet = [str(bl) if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame2]
  inVLeft = set(lVSet).difference(rVSet).difference({' '})
  dVSet = set(lVSet).intersection(rVSet).difference({' '})
  inVRight = set(rVSet).difference(lVSet).difference({' '})

  # list everything only present in the first saveFrame

  cItem1 = copy.deepcopy(cItem)
  cItem1.list.append(SAVEFRAME+saveFrame1.name)
  cItem1.inWhich = 1
  addToList(inLeft, cItem=cItem1, nefList=nefList)
  addToList(inVLeft, cItem=cItem1, nefList=nefList)

  # list everything only present in the second saveFrame

  cItem2 = copy.deepcopy(cItem)
  cItem2.list.append(SAVEFRAME+saveFrame2.name)
  cItem2.inWhich = 2
  addToList(inRight, cItem=cItem2, nefList=nefList)
  addToList(inVRight, cItem=cItem2, nefList=nefList)

  # compare the common items

  cItem3 = copy.deepcopy(cItem)
  cItem3.list.append(SAVEFRAME+saveFrame1.name)
  cItem3.inWhich = 3
  for compName in dSet:
    # compare the loop items of the matching saveFrames

    compareLoops(saveFrame1[compName], saveFrame2[compName], cItem=copy.deepcopy(cItem3), nefList=nefList)

  for compName in dVSet:
    # compare the other items in the saveFrames
    #   mandatory/optional parameters

    if saveFrame1[compName] != saveFrame2[compName]:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(SAVEFRAME+saveFrame2.name)
      cItem3.list.append('Value:  '+compName+'  -->  '\
                  +str(saveFrame1[compName])+' != '\
                  +str(saveFrame2[compName]))
      cItem3.inWhich = 3
      nefList.append(nefItem(cItem=cItem3))

  return nefList

#=========================================================================================
# compareDataBlocks
#=========================================================================================

def compareDataBlocks(dataBlock1:GenericStarParser.DataBlock
                      , dataBlock2:GenericStarParser.DataBlock
                      , cItem=nefItem()
                      , nefList=[]) -> list:
  """
  Compare two dataBlocks, if they have the same name then check their contents
  :param dataBlock1: name of the first DataBlock object
  :param dataBlock2: name of the second DataBlock object
  :return: list of type nefItem
  """
  lSet = [dataBlock1[bl].name for bl in dataBlock1]
  rSet = [dataBlock2[bl].name for bl in dataBlock2]
  inLeft = set(lSet).difference(rSet)
  dSet = set(lSet).intersection(rSet)
  inRight = set(rSet).difference(lSet)

  # list everything only present in the first DataBlock

  cItem1 = copy.deepcopy(cItem)
  cItem1.list.append(DATABLOCK+dataBlock1.name)
  cItem1.inWhich = 1
  addToList(inLeft, cItem=cItem1, nefList=nefList)

  # list everything only present in the second DataBlock

  cItem2 = copy.deepcopy(cItem)
  cItem2.list.append(DATABLOCK+dataBlock2.name)
  cItem2.inWhich = 2
  addToList(inRight, cItem=cItem2, nefList=nefList)

  # compare the common items - strictly there should only be one DataBlock

  cItem3 = copy.deepcopy(cItem)
  cItem3.list.append(DATABLOCK+dataBlock1.name)
  cItem3.inWhich = 3
  for compName in dSet:
    compareSaveFrames(dataBlock1[compName], dataBlock2[compName], cItem=copy.deepcopy(cItem3), nefList=nefList)

  return nefList

#=========================================================================================
# compareDataExtents
#=========================================================================================

def compareDataExtents(dataExt1:GenericStarParser.DataExtent
                      , dataExt2:GenericStarParser.DataExtent
                      , cItem=nefItem()
                      , nefList=[]) -> list:
  """
  Compare two dataExtents, if they have the same name then check their contents
  :param dataExt1: name of the first DataExtent object
  :param dataExt2: name of the second DataExtent object
  :return: list of type nefItem
  """
  lSet = [dataExt1[bl].name for bl in dataExt1]
  rSet = [dataExt2[bl].name for bl in dataExt2]
  inLeft = set(lSet).difference(rSet)
  dSet = set(lSet).intersection(rSet)
  inRight = set(rSet).difference(lSet)

  # list everything only present in the first DataExtent

  cItem1 = copy.deepcopy(cItem)
  cItem1.list = [DATAEXTENT+dataExt1.name]
  cItem1.inWhich = 1                                # left
  addToList(inLeft, cItem=cItem1, nefList=nefList)

  # list everything only present in the second DataExtent

  cItem2 = copy.deepcopy(cItem)
  cItem2.list = [DATAEXTENT+dataExt2.name]
  cItem2.inWhich = 2                                # right
  addToList(inRight, cItem=cItem2, nefList=nefList)

  # compare the common items - strictly there should only be one DataExtent

  cItem3 = copy.deepcopy(cItem)
  cItem3.list = [DATAEXTENT+dataExt1.name]
  cItem3.inWhich = 3                                # both
  for compName in dSet:
    compareDataBlocks(dataExt1[compName], dataExt2[compName], cItem=copy.deepcopy(cItem3), nefList=nefList)

  return nefList

#=========================================================================================
# compareFiles
#=========================================================================================

def compareNefFiles(inFile1, inFile2, cItem=nefItem(), nefList=[]):
  """
  Compare two Nef files and return comparison as a nefItem list
  :param inFile1: name of the first file
  :param inFile2: name of the second file
  :return: list of type nefItem
  """
  if not os.path.isfile(inFile1):
    print('File Error:', sys.argv[1])
  elif not os.path.isfile(inFile2):
    print('File Error:', sys.argv[2])
  else:
    try:
      NefData1 = _loadGeneralFile(path=sys.argv[1])
    except Exception as e:
      print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
      return None

    try:
      NefData2 = _loadGeneralFile(path=sys.argv[2])
    except Exception as e:
      print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
      return None

    nefList = compareDataExtents(NefData1, NefData2)  # , cItem=cItem, nefList=nefList)

  return nefList

#=========================================================================================
# __main__
#=========================================================================================

if __name__ == '__main__':
  """
  Load two files and compare

  compareNef for execution from the command line with a suitable script
  An example can be found in AnalysisV3/bin/compareNef
  
  Usage:  compareNef inFile1, inFile2     compare two Nef files
          compareNef /?                   simple instructions          
  """

  if len(sys.argv) == 2 and sys.argv[1] == '/?':
    print ('Compare contents of Nef files:')
    print ('  usage: compareNef inFile1 inFile2')
    print ('')
    print ('  Searches through all saveFrames and Loops within the files.')
    print ('  Comparisons are made for all data structures that have the same name.')
    print ('  Differences for items within a column are listed in the form:')
    print ('    dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2')
    print ('  saveFrames, Loops, columns present only in one file are listed.')
  else:
    if len(sys.argv) != 3:
      print ('Incorrect number of arguments')
    else:
      inFile1 = sys.argv[1]
      inFile2 = sys.argv[2]

      print ()
      print ('Loading...')
      nefList = compareNefFiles(inFile1, inFile2)

        # if not os.path.isfile(inFile1):
        #   print ('File Error:', sys.argv[1])
        # elif not os.path.isfile(inFile2):
        #   print('File Error:', sys.argv[2])
        # else:
        # should be okay to load the files from here

        # print ()
        # print ('Loading...')
        #
        # try:
        #   NefData1 = _loadGeneralFile(path=sys.argv[1])
        # except Exception as e:
        #   print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
        #   sys.exit()
        #
        # try:
        #   NefData2 = _loadGeneralFile(path=sys.argv[2])
        # except Exception as e:
        #   print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
        #   sys.exit()
        #
        # # NefData1 = _loadGeneralFile(path='/Users/ejb66/Downloads/cyana-3.98/demo/basic/demo.nef')
        # # NefData2 = _loadGeneralFile(path='/Users/ejb66/Downloads/cyana-3.98/demo/basic/demo2.nef')
        # # print ('  # CCPN_2l9r_Paris_155.nef')
        # # NefData3 = _loadGeneralFile(path='CCPN_2l9r_Paris_155.nef')
        #
        # cItem = nefItem()
        # nefList = []
        # # compareDataExtents(NefData1, NefData2, cItem=cItem, nefList=nefList)
        # nefList = compareDataExtents(NefData1, NefData2)    #, cItem=cItem, nefList=nefList)

      for cCount, cc in enumerate(nefList):
        if cc.inWhich == 1:
          print ('inFile1: '+ ':'.join(cc.list[:]))
        elif cc.inWhich == 2:
          print ('inFile2: '+ ':'.join(cc.list[:]))
        elif cc.inWhich == 3:
          print (':'.join(cc.list[:]))
