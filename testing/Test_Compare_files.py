"""
Module Documentation here
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = ""
__credits__ = ""
__licence__ = ("")
__reference__ = ("")
#=========================================================================================
# Last code modification:
#=========================================================================================
__modifiedBy__ = "$modifiedBy$"
__dateModified__ = "$dateModified$"
__version__ = "$Revision$"
#=========================================================================================
# Created:
#=========================================================================================
__author__ = "$Author$"
__date__ = "$Date$"
#=========================================================================================
# Start of code
#=========================================================================================

import os
import time
import copy
from .. import GenericStarParser, StarIo
from ccpn.util import Path
import unittest
import contextlib
import difflib
import json
from collections import OrderedDict

TEST_FILE_PATH = os.path.join(Path.getTopDirectory(), 'internal', 'data', 'starExamples')

class compareItem():
  def __init__(self, inDE=None, inDB=None, inSF=None, inLP=None, inCL=None, inInd=None, inStr=None, inWhich=None, cItem=None):
    # self.dExtend = inDE
    # self.dBlock = inDB
    # self.sFrame = inSF
    # self.loop = inLP
    # self.column = inCL
    # self.index = inInd
    # self.compStr = inStr
    self.inWhich = None

    if cItem is not None:
      self.list = copy.deepcopy(cItem.list)
      self.inWhich = cItem.inWhich

class compareFiles():
  def __init__(self):
    """
    Initialise an empty list of comparisons
    """
    self.numItems = 0
    self.compareList = []
    self.current = compareItem()
    self.file1 = None
    self.file2 = None

  def addItem(self, inCompare:compareItem):
    self.compareList.append(self.compareItem(inCompare))


class Test_Compare_Files(unittest.TestCase):
  """
  Test the comparison of nef files and show a diff of the results
  """

  #=========================================================================================
  # setUp     Pretest processing
  #=========================================================================================

  def setUp(self):
    """
    Preprocessing for all test cases in the class
    """
    super(Test_Compare_Files, self).setUp()
    self.bigList = []
    # self.bigList.append(compareItem())

  #=========================================================================================
  # tearDown     Posttest processing
  #=========================================================================================

  def tearDown(self):
    """
    Postprocessing for all test cases in the class
    """
    super(Test_Compare_Files, self).tearDown()

  #=========================================================================================
  # _loadGeneralFile
  #=========================================================================================

  def _loadGeneralFile(self, path=None):
    """
    Load a file with the given pathname and return a dict
    :return entry:dict
    """
    usePath = path if path.startswith('/') else os.path.join(TEST_FILE_PATH, path)
    t0 = time.time()
    entry = StarIo.parseNefFile(usePath)  # 'lenient')
    print("Parsing time %s for %s" % (time.time() - t0, path))
    return entry

  #=========================================================================================
  # test_printFile
  #=========================================================================================

  def printFile(self, thisFile):
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
  # compare data blocks
  #=========================================================================================

  def addToList(self, inList, changeList, changeHeader, name, cItem=None):
    if len(inList) > 0:
      changeList.append(changeHeader+name+' : contains --> '+','.join(inList))

    for ll in inList:
      self.bigList.append(compareItem(cItem=cItem))

  def compareLoop(self
                  , loop1:GenericStarParser.Loop
                  , loop2:GenericStarParser.Loop
                  , changeList=[]
                  , changeHeader=''
                  , cItem=None):
    """
    Compare two Loops
    """
    lSet = [bl for bl in loop1.columns]
    rSet = [bl for bl in loop2.columns]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(loop1.name)
    cItem1.inWhich = 1      # left

    self.addToList(inLeft, changeList, '---LEFT  '+changeHeader, loop1.name, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(loop2.name)
    cItem2.inWhich = 2      # left

    self.addToList(inRight, changeList, '---RIGHT '+changeHeader, loop2.name, cItem=cItem2)

    if loop1.data and loop2.data:
      for compName in dSet:
        rowRange = min(len(loop1.data), len(loop2.data))
        for rowIndex in range(rowRange):
          if loop1.data[rowIndex][compName] != loop2.data[rowIndex][compName]:
            changeStr = 'colDiff: '\
                        +changeHeader\
                        +compName+' : '\
                        +str(rowIndex)+' --> '\
                        +str(loop1.data[rowIndex][compName])+' != '\
                        +str(loop2.data[rowIndex][compName])

            changeList.append(changeStr)
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append(loop1.name)
            cItem3.list.append(compName)
            cItem3.list.append(str(rowIndex))
            cItem3.list.append(str(loop1.data[rowIndex][compName]))
            cItem3.list.append(str(loop2.data[rowIndex][compName]))
            cItem3.list.append(compName+' : '\
                        +str(rowIndex)+' --> '\
                        +str(loop1.data[rowIndex][compName])+' != '\
                        +str(loop2.data[rowIndex][compName]))
            cItem3.inWhich = 3
            self.bigList.append(compareItem(cItem=cItem3))


  def compareSaveFrame(self
                        , saveFrame1:GenericStarParser.SaveFrame
                        , saveFrame2:GenericStarParser.SaveFrame
                        , changeList=[]
                        , changeHeader=''
                        , cItem=None):
    """
    Compare two saveFrames, if they have the same name then check their contents
    """
    lSet = [None if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else saveFrame1[bl].name for bl in saveFrame1]
    rSet = [None if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else saveFrame2[bl].name for bl in saveFrame2]
    inLeft = set(lSet).difference(rSet).difference({None})
    dSet = set(lSet).intersection(rSet).difference({None})     # get rid of None
    inRight = set(rSet).difference(lSet).difference({None})

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(saveFrame1.name)
    cItem1.inWhich = 1      # left

    self.addToList(inLeft, changeList, '--LEFT  '+changeHeader, saveFrame1.name, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(saveFrame2.name)
    cItem2.inWhich = 2      # left

    self.addToList(inRight, changeList, '--RIGHT '+changeHeader, saveFrame2.name, cItem=cItem2)

    for compName in dSet:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(saveFrame1.name)
      cItem3.list.append(compName)
      cItem3.inWhich = 3  # right

      self.compareLoop(saveFrame1[compName], saveFrame2[compName], changeList, changeHeader+ '\n      ' + compName + ':', cItem=copy.deepcopy(cItem3))

  def compareDataBlock(self
                        , dataBlock1:GenericStarParser.DataBlock
                        , dataBlock2:GenericStarParser.DataBlock
                        , changeList=[]
                        , changeHeader=''
                        , cItem=None):
    """
    Compare two dataBlocks, if they have the same name then check their contents
    """
    lSet = [dataBlock1[bl].name for bl in dataBlock1]
    rSet = [dataBlock2[bl].name for bl in dataBlock2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append(dataBlock1.name)
    cItem1.inWhich = 1      # left

    self.addToList(inLeft, changeList, '-LEFT  '+changeHeader, dataBlock1.name, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append(dataBlock2.name)
    cItem2.inWhich = 2      # right

    self.addToList(inRight, changeList, '-RIGHT '+changeHeader, dataBlock2.name, cItem=cItem2)

    for compName in dSet:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append(dataBlock1.name)
      cItem3.list.append(compName)
      cItem3.inWhich = 3  # right

      self.compareSaveFrame(dataBlock1[compName], dataBlock2[compName], changeList, changeHeader+ '\n    ' + compName + ':', cItem=copy.deepcopy(cItem3))

  def compareDataExtent(self
                        , dataExt1:GenericStarParser.DataExtent
                        , dataExt2:GenericStarParser.DataExtent
                        , changeList=[]
                        , changeHeader=''
                        , cItem=None):
    """
    Compare two dataExtents, if they have the same name then check their contents
    """
    lSet = [dataExt1[bl].name for bl in dataExt1]
    rSet = [dataExt2[bl].name for bl in dataExt2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    cItem1 = copy.deepcopy(cItem)
    cItem1.list = [dataExt1.name]
    cItem1.inWhich = 1      # left

    self.addToList(inLeft, changeList, 'LEFT  '+changeHeader, dataExt1.name, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list = [dataExt2.name]
    cItem2.inWhich = 2      # right

    self.addToList(inRight, changeList, 'RIGHT '+changeHeader, dataExt2.name, cItem=cItem2)

    for compName in dSet:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list = [dataExt1.name]
      cItem3.list.append(compName)
      cItem3.inWhich = 3  # right

      self.compareDataBlock(dataExt1[compName], dataExt2[compName], changeList, changeHeader+ '\n  ' + compName + ':', cItem=copy.deepcopy(cItem3))

  #=========================================================================================
  # test_Compare_Files
  #=========================================================================================

  def test_Compare_Files(self):
    """
    Load two files and compare
    """
    print ('Loading...')
    print ('  # Commented_Example.nef')
    NefData1 = self._loadGeneralFile(path='Commented_Example.nef')
    print ('  # Commented_Example_Change.nef')
    NefData2 = self._loadGeneralFile(path='Commented_Example_Change.nef')

    # print ('  # CCPN_2l9r_Paris_155.nef')
    # NefData3 = self._loadGeneralFile(path='CCPN_2l9r_Paris_155.nef')

    print ('~'*80)
    compareSet = []
    cItem = compareItem()
    self.compareDataExtent(NefData1, NefData2, compareSet, cItem=cItem)
    if len(compareSet) > 0:
      for x in compareSet:
        print (x)