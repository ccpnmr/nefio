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
# from .. import GenericStarParser, StarIo
from ccpn.util.nef import GenericStarParser, StarIo
from ccpn.util import Path
import unittest
import contextlib
import difflib
import json
from collections import OrderedDict

TEST_FILE_PATH = os.path.join(Path.getTopDirectory(), 'internal', 'data', 'starExamples')

class compareItem():
  def __init__(self, cItem=None):
    self.inWhich = None

    if cItem is not None:
      self.list = copy.deepcopy(cItem.list)
      self.inWhich = cItem.inWhich

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

  def addToList(self, inList, cItem=None):
    if len(inList) > 0:
      cItem3 = copy.deepcopy(cItem)
      cItem3.list.append('contains --> '+', '.join(inList))
      self.bigList.append(compareItem(cItem=cItem3))

  def compareLoop(self
                  , loop1:GenericStarParser.Loop
                  , loop2:GenericStarParser.Loop
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
    cItem1.list.append('Loop:'+loop1.name)
    cItem1.inWhich = 1
    self.addToList(inLeft, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append('Loop:'+loop2.name)
    cItem2.inWhich = 2
    self.addToList(inRight, cItem=cItem2)

    if loop1.data and loop2.data:
      for compName in dSet:
        rowRange = min(len(loop1.data), len(loop2.data))
        for rowIndex in range(rowRange):
          if loop1.data[rowIndex][compName] != loop2.data[rowIndex][compName]:
            cItem3 = copy.deepcopy(cItem)
            cItem3.list.append('Loop:'+loop1.name)
            cItem3.list.append('<Column>: '+compName+'  <rowIndex>: '\
                        +str(rowIndex)+'  -->  '\
                        +str(loop1.data[rowIndex][compName])+' != '\
                        +str(loop2.data[rowIndex][compName]))
            cItem3.inWhich = 3
            self.bigList.append(compareItem(cItem=cItem3))

  def compareSaveFrame(self
                        , saveFrame1:GenericStarParser.SaveFrame
                        , saveFrame2:GenericStarParser.SaveFrame
                        , cItem=None):
    """
    Compare two saveFrames, if they have the same name then check their contents
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

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append('SaveFrame:'+saveFrame1.name)
    cItem1.inWhich = 1
    self.addToList(inLeft, cItem=cItem1)
    self.addToList(inVLeft, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append('SaveFrame:'+saveFrame2.name)
    cItem2.inWhich = 2
    self.addToList(inRight, cItem=cItem2)
    self.addToList(inVRight, cItem=cItem2)

    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append('SaveFrame:'+saveFrame1.name)
    cItem3.inWhich = 3
    for compName in dSet:
      self.compareLoop(saveFrame1[compName], saveFrame2[compName], cItem=copy.deepcopy(cItem3))

    for compName in dVSet:
      if saveFrame1[compName] != saveFrame2[compName]:
        cItem3 = copy.deepcopy(cItem)
        cItem3.list.append('SaveFrame:' + saveFrame2.name)
        cItem3.list.append('Value:  '+compName+'  -->  '\
                    +str(saveFrame1[compName])+' != '\
                    +str(saveFrame2[compName]))
        cItem3.inWhich = 3
        self.bigList.append(compareItem(cItem=cItem3))

  def compareDataBlock(self
                        , dataBlock1:GenericStarParser.DataBlock
                        , dataBlock2:GenericStarParser.DataBlock
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
    cItem1.list.append('DataBlock:'+dataBlock1.name)
    cItem1.inWhich = 1
    self.addToList(inLeft, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append('DataBlock:'+dataBlock2.name)
    cItem2.inWhich = 2
    self.addToList(inRight, cItem=cItem2)

    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append('DataBlock:' + dataBlock1.name)
    cItem3.inWhich = 3
    for compName in dSet:
      self.compareSaveFrame(dataBlock1[compName], dataBlock2[compName], cItem=copy.deepcopy(cItem3))

  def compareDataExtent(self
                        , dataExt1:GenericStarParser.DataExtent
                        , dataExt2:GenericStarParser.DataExtent
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
    cItem1.list = ['DataExtent:'+dataExt1.name]
    cItem1.inWhich = 1                                # left
    self.addToList(inLeft, cItem=cItem1)

    cItem2 = copy.deepcopy(cItem)
    cItem2.list = ['DataExtent:'+dataExt2.name]
    cItem2.inWhich = 2                                # right
    self.addToList(inRight, cItem=cItem2)

    cItem3 = copy.deepcopy(cItem)
    cItem3.list = ['DataExtent:' + dataExt1.name]
    cItem3.inWhich = 3                                # both
    for compName in dSet:
      self.compareDataBlock(dataExt1[compName], dataExt2[compName], cItem=copy.deepcopy(cItem3))

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

    cItem = compareItem()

    #                       file1='Commented_Example.nef'
    #                     , file2='Commented_Example_Change.nef'
    #                     , NefData1=NefData1
    #                     , NefData2=NefData2)

    self.compareDataExtent(NefData1, NefData2, cItem=cItem)

    for cCount, cc in enumerate(self.bigList):
      print ('~'*80)
      print(str(cCount) + ': ' + 'difference')
      if cc.inWhich == 1:
        print ('Present in LEFT')
      elif cc.inWhich == 2:
        print ('Present in RIGHT')
      elif cc.inWhich == 3:
        print ('Present in Both')

      for ind, ll in enumerate(cc.list):
        print ('.  '*ind, ll)
