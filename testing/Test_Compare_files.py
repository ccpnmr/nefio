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
from .. import GenericStarParser, StarIo
from ccpn.util import Path
import unittest
import contextlib
import difflib
import json
from collections import OrderedDict

TEST_FILE_PATH = os.path.join(Path.getTopDirectory(), 'internal', 'data', 'starExamples')

class compareItem():
  def __init__(self, inDE=None, inDB=None, inSF=None, inLP=None, inCL=None, inStr=None):
    self.dExtend = inDE
    self.dBlock = inDB
    self.sFrame = inSF
    self.loop = inLP
    self.column = inCL
    self.compStr = inStr
    self.inWhich = 0

class compareFiles():
  def __init__(self):
    """
    Initialise an empty list of comparisons
    """
    self.numItems = 0
    self.compareList = []
    self.current = compareItem()

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

  def addToList(self, inList, changeList, changeHeader, name):
    if len(inList) > 0:
      changeList.append(changeHeader+name+' : contains --> '+','.join(inList))

  def compareLoop(self
                  , loop1:GenericStarParser.Loop
                  , loop2:GenericStarParser.Loop
                  , changeList=[]
                  , changeHeader=''):
    """
    Compare two Loops
    """
    lSet = [bl for bl in loop1.columns]
    rSet = [bl for bl in loop2.columns]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    self.addToList(inLeft, changeList, '---LEFT  '+changeHeader, loop1.name)
    self.addToList(inRight, changeList, '---RIGHT '+changeHeader, loop2.name)

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

  def compareSaveFrame(self
                        , saveFrame1:GenericStarParser.SaveFrame
                        , saveFrame2:GenericStarParser.SaveFrame
                        , changeList=[]
                        , changeHeader=''):
    """
    Compare two saveFrames, if they have the same name then check their contents
    """
    lSet = [None if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else saveFrame1[bl].name for bl in saveFrame1]
    rSet = [None if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else saveFrame2[bl].name for bl in saveFrame2]
    inLeft = set(lSet).difference(rSet).difference({None})
    dSet = set(lSet).intersection(rSet).difference({None})     # get rid of None
    inRight = set(rSet).difference(lSet).difference({None})

    self.addToList(inLeft, changeList, '--LEFT  '+changeHeader, saveFrame1.name)
    self.addToList(inRight, changeList, '--RIGHT '+changeHeader, saveFrame2.name)

    for compName in dSet:
      self.compareLoop(saveFrame1[compName], saveFrame2[compName], changeList, changeHeader+ '\n      ' + compName + ':')

  def compareDataBlock(self
                        , dataBlock1:GenericStarParser.DataBlock
                        , dataBlock2:GenericStarParser.DataBlock
                        , changeList=[]
                        , changeHeader=''):
    """
    Compare two dataBlocks, if they have the same name then check their contents
    """
    lSet = [dataBlock1[bl].name for bl in dataBlock1]
    rSet = [dataBlock2[bl].name for bl in dataBlock2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    self.addToList(inLeft, changeList, '-LEFT  '+changeHeader, dataBlock1.name)
    self.addToList(inRight, changeList, '-RIGHT '+changeHeader, dataBlock2.name)

    for compName in dSet:
      self.compareSaveFrame(dataBlock1[compName], dataBlock2[compName], changeList, changeHeader+ '\n    ' + compName + ':')

  def compareDataExtent(self
                        , dataExt1:GenericStarParser.DataExtent
                        , dataExt2:GenericStarParser.DataExtent
                        , changeList=[]
                        , changeHeader=''):
    """
    Compare two dataExtents, if they have the same name then check their contents
    """
    lSet = [dataExt1[bl].name for bl in dataExt1]
    rSet = [dataExt2[bl].name for bl in dataExt2]
    inLeft = set(lSet).difference(rSet)
    dSet = set(lSet).intersection(rSet)
    inRight = set(rSet).difference(lSet)

    self.addToList(inLeft, changeList, 'LEFT  '+changeHeader, dataExt1.name)
    self.addToList(inRight, changeList, 'RIGHT '+changeHeader, dataExt2.name)

    for compName in dSet:
      self.compareDataBlock(dataExt1[compName], dataExt2[compName], changeList, changeHeader+ '\n  ' + compName + ':')

  #=========================================================================================
  # test_Compare_Files
  #=========================================================================================

  def test_Compare_Files(self):
    """
    Load two files and compare
    """
    compList = compareFiles()

    print ('Loading...')
    print ('  # Commented_Example.nef')
    NefData1 = self._loadGeneralFile(path='Commented_Example.nef')
    print ('  # Commented_Example_Change.nef')
    NefData2 = self._loadGeneralFile(path='Commented_Example_Change.nef')

    # print ('  # CCPN_2l9r_Paris_155.nef')
    # NefData3 = self._loadGeneralFile(path='CCPN_2l9r_Paris_155.nef')

    print ('~'*80)
    compareSet = []
    self.compareDataExtent(NefData1, NefData2, compareSet)
    if len(compareSet) > 0:
      for x in compareSet:
        print (x)