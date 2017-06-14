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

    Usage:  compareNef inFile1 inFile2            compare two Nef files
                                                  print results to the screen
            compareNef -b inDir1 inDir2 outDir    compare Nef files common to directories
                                                  inDir1 and inDir2. Write output *.txt for each
                                                  file into the outDir directory.
            compareNef /?                         simple instructions

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
                                    the string value can be a dictionary

Module Contents
===============

compareNef contains the following routines:
  compareNefFiles      compare two Nef files and return a compare object

    Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
    Comparisons are made for all data structures that have the same name.
    Differences for items within a column are listed in the form:
      dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2

    dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed, in the form:
      dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                        parameter2
                                                        ...
                                                        parameterN

  A compare object is a list of nefItems of the form:

  ::

      NefItem
        inWhich         a flag labelling which file the item was found in
                        1 = found in first file, 2 = found in second file, 3 = common to both
        List
          Item          multiple strings containing the comparison tree
          (,List)       the last item of which may be a list of items common to the tree

  e.g., for parameters present in the first file:
        [
          inWhich=1
          list=[dataExtent1, dataBlock1, saveFrame1, Loop1, [parameter1, parameter2, parameter3]]
        ]

  compareDataExtents    compare two DataExtent objects and return a compare list as above
  compareDataBlocks     compare two DataBlock objects and return a compare list as above
  compareSaveFrames     compare two SaveFrame objects and return a compare list as above
  compareLoops          compare two Loop objects and return a compare list as above

  printCompareList      print the compare list to the screen
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
import unittest
from typing import Optional
# import json
# from collections import OrderedDict
# import contextlib
# import difflib
from ast import literal_eval
from os import listdir
from os.path import isfile, join

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
# sizeNefList
#=========================================================================================

def sizeNefList(nefList, whichType=0) -> int:
  """
  List only those items that are of type whichType
  :param nefList: list to print
  :param whichType: type to print
  """
  count=0
  if nefList is not None:
    for cCount, cc in enumerate(nefList):
      if cc.inWhich == whichType:
        count += 1
  return count

#=========================================================================================
# printWhichList
#=========================================================================================

def printWhichList(nefList, whichType=0):
  """
  List only those items that are of type whichType
  :param nefList: list to print
  :param whichType: type to print
  """
  for cCount, cc in enumerate(nefList):
    if cc.inWhich == whichType:

      if isinstance(cc.list[-1], str):
        print('  ' + ':'.join(cc.list[:]))
      else:
        outStr = '  '+':'.join(cc.list[:-1])+': contains --> '
        lineTab = '\n'+' '*len(outStr)
        print (outStr+lineTab.join(cc.list[-1]))

#=========================================================================================
# printCompareList
#=========================================================================================

def printCompareList(nefList, inFile1, inFile2):
  """
  Print the contents of the nef compare list to the screen

  Output is in three parts:
    - items that are present only in the first file
    - items that are only in the second file
    - differences between objects that are common in both files

  :param nefList: list to print
  :param inFile1: name of the first file
  :param inFile2: name of the second file
  """

  if not isinstance(inFile1, str):
    print ('TypeError: inFile1 must be a string.')
    return
  if not isinstance(inFile2, str):
    print ('TypeError: inFile2 must be a string.')
    return

  if sizeNefList(nefList, whichType=1) > 0:
    print ('\nItems that are only present in '+inFile1+':')
    printWhichList(nefList, 1)

  if sizeNefList(nefList, whichType=2) > 0:
    print ('\nItems that are only present in '+inFile2+':')
    printWhichList(nefList, 2)

  if sizeNefList(nefList, whichType=3) > 0:
    print ('\nItems that are present in both files:')
    printWhichList(nefList, 3)

#=========================================================================================
# addToList
#=========================================================================================

def addToList(inList, cItem, nefList) -> list:
  """
  Append cItem to the compare list
  Currently adds one cItem with a list as the last element

  :param inList: a list of items to add to the end of cItem
  :param cItem: object containing the current tree to add to the list
  :param nefList: current list of comparisons
  :return: list of type nefItem
  """
  if len(inList) > 0:
    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append(list(inList))              # nest the list within the cItem
    nefList.append(nefItem(cItem=cItem3))

  # if len(inList) > 0:
  #   cItem3 = copy.deepcopy(cItem)
  #   lineTab = len('contains --> ')+len(cItem.list)+1
  #   for cc in cItem.list:
  #     lineTab += len(cc)
  #   preStr = '\n'+' '*lineTab
  #   cItem3.list.append('contains --> '+preStr.join(inList))
  #   nefList.append(nefItem(cItem=cItem3))

  # for eachItem in inList:
  #   cItem3 = copy.deepcopy(cItem)
  #   cItem3.list.append('contains --> '+eachItem)
  #   nefList.append(nefItem(cItem=cItem3))

  return nefList

#=========================================================================================
# compareLoops
#=========================================================================================

def compareLoops(loop1:GenericStarParser.Loop
                , loop2:GenericStarParser.Loop
                , cItem=None
                , nefList=None) -> list:
  """
  Compare two Loops
  :param loop1: name of the first Loop object
  :param loop2: name of the second Loop object
  :return: list of type nefItem
  """
  if cItem is None:
    cItem = nefItem()
  if nefList is None:
    nefList = []

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

    if len(loop1.data) != len(loop2.data):        # simple compare, same length tables
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop1.name)
      cItem3.list.append(' <rowLength>:  '+str(len(loop1.data))+' != '+str(len(loop2.data)))
      cItem3.inWhich = 3
      nefList.append(nefItem(cItem=cItem3))

    # carry on and compare the common table

    for compName in dSet:
      for rowIndex in range(rowRange):

        loopValue1 = loop1.data[rowIndex][compName]
        loopValue2 = loop2.data[rowIndex][compName]

        if loopValue1 != loopValue2:

          # The value_strings are different
          # Check to see if they are dictionaries
          # and compare contents

          try:
            # these 2 lines will crash as a ValueError: malformed string if cannot
            # be evaluated as:
            #   - strings, numbers, tuples, lists, dicts, booleans, and None.

            loopValue1 = literal_eval(loopValue1)
            loopValue2 = literal_eval(loopValue2)
            assert isinstance(loopValue1, dict)
            assert isinstance(loopValue2, dict)

            if isinstance(loopValue1, dict) and isinstance(loopValue2, dict):
              if loopValue1 != loopValue2:

                # may need a deeper compare of inserted dictionaries in here
                cItem3 = copy.deepcopy(cItem)
                cItem3.list.append(LOOP + loop1.name)
                cItem3.list.append(' <Column>: ' + compName + '  <rowIndex>: ' \
                                   + str(rowIndex) + '  -->  ' \
                                   + str(loopValue1) + ' != ' \
                                   + str(loopValue2))
                cItem3.inWhich = 3
                nefList.append(nefItem(cItem=cItem3))
            else:
              # not both dicts so compare as normal strings

              cItem3 = copy.deepcopy(cItem)
              cItem3.list.append(LOOP+loop1.name)
              cItem3.list.append(' <Column>: '+compName+'  <rowIndex>: '\
                          +str(rowIndex)+'  -->  '\
                          +str(loopValue1)+' != '\
                          +str(loopValue2))
              cItem3.inWhich = 3
              nefList.append(nefItem(cItem=cItem3))

          except (SyntaxError, ValueError, AssertionError):
            # loopvalues cannot be converted to proper values
            # need to check that comments are being loaded correctly

            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(LOOP + loop1.name)
            cItem3.list.append(' <Column>: ' + compName + '  <rowIndex>: ' \
                               + str(rowIndex) + '  -->  ' \
                               + str(loopValue1) + ' != ' \
                               + str(loopValue2))
            cItem3.inWhich = 3
            nefList.append(nefItem(cItem=cItem3))

        else:

          # nothing for the minute as identical already but may want to keep a log
          pass

    # else:
    #   cItem3 = copy.deepcopy(cItem)
    #   cItem3.list.append(LOOP+loop1.name)
    #   cItem3.list.append(' <rowLength>:  '+str(len(loop1.data))+' != '+str(len(loop2.data)))
    #   cItem3.inWhich = 3
    #   nefList.append(nefItem(cItem=cItem3))

      #TODO
      # need to add a further test here, could do a diff on the tables which would pick up
      # insertions to the table - this columns would need to be reordered for this to work
      # what if there are a different number of columns?
      # also check for Mandatory items

  else:
    # can't compare non-existent loopdata
    if loop1.data is None:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop1.name)
      cItem3.list.append(' <Contains no data>')
      cItem3.inWhich = 1
      nefList.append(nefItem(cItem=cItem3))
    if loop2.data is None:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(LOOP+loop2.name)
      cItem3.list.append(' <Contains no data>')
      cItem3.inWhich = 2
      nefList.append(nefItem(cItem=cItem3))

  return nefList

#=========================================================================================
# compareSaveFrames
#=========================================================================================

def compareSaveFrames(saveFrame1:GenericStarParser.SaveFrame
                      , saveFrame2:GenericStarParser.SaveFrame
                      , cItem=None
                      , nefList=None) -> list:
  """
  Compare two saveFrames, if they have the same name then check their contents
  :param saveFrame1: name of the first SaveFrame object
  :param saveFrame2: name of the second SaveFrame object
  :return: list of type nefItem
  """
  if cItem is None:
    cItem = nefItem()
  if nefList is None:
    nefList = []

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

    compareLoops(saveFrame1[compName], saveFrame2[compName], cItem=cItem3, nefList=nefList)

  for compName in dVSet:
    # compare the other items in the saveFrames
    #   mandatory/optional parameters

    if saveFrame1[compName] != saveFrame2[compName]:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(SAVEFRAME+saveFrame2.name)
      cItem3.list.append(' <Value>:  '+compName+'  -->  '\
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
                      , cItem=None
                      , nefList=None) -> list:
  """
  Compare two dataBlocks, if they have the same name then check their contents
  :param dataBlock1: name of the first DataBlock object
  :param dataBlock2: name of the second DataBlock object
  :return: list of type nefItem
  """
  if cItem is None:
    cItem = nefItem()
  if nefList is None:
    nefList = []

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
    compareSaveFrames(dataBlock1[compName], dataBlock2[compName], cItem=cItem3, nefList=nefList)

  return nefList

#=========================================================================================
# compareDataExtents
#=========================================================================================

def compareDataExtents(dataExt1:GenericStarParser.DataExtent
                      , dataExt2:GenericStarParser.DataExtent
                      , cItem=None
                      , nefList=None) -> list:
  """
  Compare two dataExtents, if they have the same name then check their contents
  :param dataExt1: name of the first DataExtent object
  :param dataExt2: name of the second DataExtent object
  :return: list of type nefItem
  """
  if cItem is None:
    cItem = nefItem()
  if nefList is None:
    nefList = []

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
    compareDataBlocks(dataExt1[compName], dataExt2[compName], cItem=cItem3, nefList=nefList)

  return nefList

#=========================================================================================
# compareFiles
#=========================================================================================

def compareNefFiles(inFile1, inFile2, cItem=None, nefList=None) -> Optional[list]:
  """
  Compare two Nef files and return comparison as a nefItem list
  :param inFile1: name of the first file
  :param inFile2: name of the second file
  :return: list of type nefItem
  """
  if cItem is None:
    cItem = nefItem()
  if nefList is None:
    nefList = []

  if not os.path.isfile(inFile1):
    print('File Error:', inFile1)
  elif not os.path.isfile(inFile2):
    print('File Error:', inFile2)
  else:
    try:
      NefData1 = _loadGeneralFile(path=inFile1)
    except Exception as e:
      print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
      return None

    try:
      NefData2 = _loadGeneralFile(path=inFile2)
    except Exception as e:
      print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
      return None

    compareDataExtents(NefData1, NefData2, cItem=cItem, nefList=nefList)

  return nefList

#=========================================================================================
# Test_Compare_Files
#=========================================================================================

def batchCompareNefFiles(inDir1, inDir2, outDir):
  """
  Batch compare the Nef files common to the two directories
  For each file found, write the compare log to the corresponding .txt file

  :param inDir1:
  :param inDir2:
  :param outDir:
  """
  inFileList = [f for f in listdir(inDir1) if isfile(join(inDir1, f)) and f[-4:] == '.nef']
  outFileList = [f for f in listdir(inDir2) if isfile(join(inDir2, f)) and f[-4:] == '.nef']

  for fl in inFileList:
    if fl in outFileList:

      outFileName = join(outDir, fl[:-4] + '.txt')
      print('Batch processing %s > %s' % (fl, outFileName))

      stdOriginal = sys.stdout

      nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl))
      with open(outFileName, 'w') as outLog:
        sys.stdout = outLog
        print (join(inDir1, fl))
        print (join(inDir2, fl))
        printCompareList(nefList, join(inDir1, fl), join(inDir2, fl))
        sys.stdout = stdOriginal

#=========================================================================================
# Test_Compare_Files
#=========================================================================================

class Test_Compare_Files(unittest.TestCase):
  """
  Test the comparison of nef files and print the results
  """
  def _test_Compare_Files(self):
    """
    Load two files and compare
    """
    # test set one
    inFile1 = '/Users/ejb66/PycharmProjects/AnalysisV3/internal/data/starExamples/Commented_Example.nef'
    inFile2 = '/Users/ejb66/PycharmProjects/AnalysisV3/internal/data/starExamples/Commented_Example_Change.nef'

    print ('\nTEST COMPARISON')
    print ('   file1 = '+inFile1)
    print ('   file2 = '+inFile2)
    print ('Loading...')
    nefList = compareNefFiles(inFile1, inFile2)
    printCompareList(nefList, inFile1, inFile2)

    # test set two
    inFile1 = '/Users/ejb66/Desktop/Temporary/1nk2_docr_extended.ccpn.nef'
    inFile2 = '/Users/ejb66/PycharmProjects/AnalysisV3/internal/data/NEF_test_data/NMRX/1nk2_docr_extended.ccpn.nef'

    print ('~'*80)
    print ('\nTEST COMPARISON')
    print ('   file1 = '+inFile1)
    print ('   file2 = '+inFile2)
    print ('Loading...')
    nefList = compareNefFiles(inFile1, inFile2)
    printCompareList(nefList, inFile1, inFile2)

    inFile1 = '/Users/ejb66/Desktop/Temporary/2mqq_docr_extended.nef'
    inFile2 = '/Users/ejb66/PycharmProjects/AnalysisV3/internal/data/NEF_test_data/NMRX/2mqq_docr_extended.nef'

    print ('~'*80)
    print ('\nTEST COMPARISON')
    print ('   file1 = '+inFile1)
    print ('   file2 = '+inFile2)
    print ('Loading...')
    nefList = compareNefFiles(inFile1, inFile2)
    printCompareList(nefList, inFile1, inFile2)


  def test_Compare_BatchFiles(self):
    """
    Compare the Nef files in two directories
    """
    inDir1 = '/Users/ejb66/Desktop/Temporary'
    inDir2 = '/Users/ejb66/Dropbox/CCPNdocsShared/NefTestData/v3'
    outDir = '/Users/ejb66/Desktop/Temporary'

    batchCompareNefFiles(inDir1, inDir2, outDir)

#=========================================================================================
# __main__
#=========================================================================================

if __name__ == '__main__':
  """
  Command Line Usage:
    compareNef for execution from the command line with a suitable script
    An example can be found in AnalysisV3/bin/compareNef:
  
        #!/usr/bin/env sh
        export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
        export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
        export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
        ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/compareNef.py $*
  
    Usage:  compareNef inFile1 inFile2            compare two Nef files
                                                  print results to the screen
            compareNef -b inDir1 inDir2 outDir    compare Nef files common to directories
                                                  inDir1 and inDir2. Write output *.txt for each
                                                  file into the outDir directory.
            compareNef /?                         simple instructions
                                                  
    Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
    Comparisons are made for all data structures that have the same name.
    Differences for items within a column are listed in the form:
      dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2
      
    dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed.
      dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                        parameter2
                                                        ...
                                                        parameterN
  """

  if len(sys.argv) == 2 and sys.argv[1] == '/?':
    print ('Command Line Usage:')
    print ('  compareNef for execution from the command line with a suitable script')
    print ('  An example can be found in AnalysisV3/bin/compareNef:')
    print ('')
    print ('      #!/usr/bin/env sh')
    print ('      export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"')
    print ('      export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda')
    print ('      export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c')
    print ('      ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/compareNef.py $*')
    print ('')
    print ('  Usage:  compareNef inFile1 inFile2            compare two Nef files')
    print ('                                                print results to the screen')
    print ('          compareNef -b inDir1 inDir2 outDir    compare Nef files common to directories')
    print ('                                                inDir1 and inDir2. Write output *.txt for each')
    print ('                                                file into the outDir directory.')
    print ('          compareNef /?                         simple instructions')
    print ('')
    print ('  Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.')
    print ('  Comparisons are made for all data structures that have the same name.')
    print ('  Differences for items within a column are listed in the form:')
    print ('    dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2')
    print ('')
    print ('  dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed.')
    print ('    dataExtent: dataBlock:saveFrame: Loop: contains --> parameter1')
    print ('                                                        parameter2')
    print ('                                                        ...')
    print ('                                                        parameterN')
  else:
    if len(sys.argv) == 3:          # assume compareNef inFile1 inFile2
      inFile1 = sys.argv[1]
      inFile2 = sys.argv[2]

      print ()
      print ('Loading...')
      nefList = compareNefFiles(inFile1, inFile2)
      printCompareList(nefList, inFile1, inFile2)

    elif len(sys.argv) == 5 and sys.argv[1] == '-b':     # assume compareNef -b inDir1 inDir2 outDir
      inDir1 = sys.argv[2]
      inDir2 = sys.argv[3]
      outDir = sys.argv[4]

      batchCompareNefFiles(inDir1, inDir2, outDir)

    else:
      print ('Incorrect arguments')
