"""
  This has now been moved to file ../CompareNef.py
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license",
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license",
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: CCPN $"
__dateModified__ = "$dateModified: 2017-07-07 16:33:02 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b3 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
import time
import copy
import sys
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
      rowRange = min(len(loop1.data), len(loop2.data))
      if len(loop1.data) == len(loop2.data):        # simple compare, same length tables
        for compName in dSet:
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
      else:
        cItem3 = copy.deepcopy(cItem)
        cItem3.list.append('Loop:' + loop1.name)
        cItem3.list.append('<rowLength>:  '+str(len(loop1.data))+' != '+str(len(loop2.data)))
        cItem3.inWhich = 3
        self.bigList.append(compareItem(cItem=cItem3))

        #TODO
        # need to add a further test here, could do a diff on the tables which would pick up
        # insertions to the table - this columns would need to be reordered for this to work
        # what if there are a different number of columns?
        # also check for Mandatory items
        #
    else:
      # can't compare non-existant loopdata
      if loop1.data is None:
        cItem3 = copy.deepcopy(cItem)
        cItem3.list.append('Loop:' + loop1.name)
        cItem3.list.append('<Contains no data>')
        cItem3.inWhich = 1
        self.bigList.append(compareItem(cItem=cItem3))
      if loop2.data is None:
        cItem3 = copy.deepcopy(cItem)
        cItem3.list.append('Loop:' + loop2.name)
        cItem3.list.append('<Contains no data>')
        cItem3.inWhich = 2
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

    # list everything only present in the first saveFrame

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append('SaveFrame:'+saveFrame1.name)
    cItem1.inWhich = 1
    self.addToList(inLeft, cItem=cItem1)
    self.addToList(inVLeft, cItem=cItem1)

    # list everything only present in the second saveFrame

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append('SaveFrame:'+saveFrame2.name)
    cItem2.inWhich = 2
    self.addToList(inRight, cItem=cItem2)
    self.addToList(inVRight, cItem=cItem2)

    # compare the common items

    cItem3 = copy.deepcopy(cItem)
    cItem3.list.append('SaveFrame:'+saveFrame1.name)
    cItem3.inWhich = 3
    for compName in dSet:
      # compare the loop items of the matching saveFrames

      self.compareLoop(saveFrame1[compName], saveFrame2[compName], cItem=copy.deepcopy(cItem3))

    for compName in dVSet:
      # compare the other items in the saveFrames
      #   mandatory/optional parameters

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

    # list everything only present in the first DataBlock

    cItem1 = copy.deepcopy(cItem)
    cItem1.list.append('DataBlock:'+dataBlock1.name)
    cItem1.inWhich = 1
    self.addToList(inLeft, cItem=cItem1)

    # list everything only present in the second DataBlock

    cItem2 = copy.deepcopy(cItem)
    cItem2.list.append('DataBlock:'+dataBlock2.name)
    cItem2.inWhich = 2
    self.addToList(inRight, cItem=cItem2)

    # compare the common items - strictly there should only be one DataBlock

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

    # list everything only present in the first DataExtent

    cItem1 = copy.deepcopy(cItem)
    cItem1.list = ['DataExtent:'+dataExt1.name]
    cItem1.inWhich = 1                                # left
    self.addToList(inLeft, cItem=cItem1)

    # list everything only present in the second DataExtent

    cItem2 = copy.deepcopy(cItem)
    cItem2.list = ['DataExtent:'+dataExt2.name]
    cItem2.inWhich = 2                                # right
    self.addToList(inRight, cItem=cItem2)

    # compare the common items - strictly there should only be one DataExtent

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

    # NefData1 = self._loadGeneralFile(path='/Users/ejb66/Downloads/cyana-3.98/demo/basic/demo.nef')
    # NefData2 = self._loadGeneralFile(path='/Users/ejb66/Downloads/cyana-3.98/demo/basic/demo2.nef')

    # print ('  # CCPN_2l9r_Paris_155.nef')
    # NefData3 = self._loadGeneralFile(path='CCPN_2l9r_Paris_155.nef')

    cItem = compareItem()

    #                       file1='Commented_Example.nef'
    #                     , file2='Commented_Example_Change.nef'
    #                     , NefData1=NefData1
    #                     , NefData2=NefData2)

    self.compareDataExtent(NefData1, NefData2, cItem=cItem)

    for cCount, cc in enumerate(self.bigList):
      # print ('~'*80)
      # print(str(cCount) + ': ' + 'difference')
      # if cc.inWhich == 1:
      #   print ('Present in First file\n')
      # elif cc.inWhich == 2:
      #   print ('Present in Second file\n')
      # elif cc.inWhich == 3:
      #   print ('Present in Both\n')

      # for ind, ll in enumerate(cc.list):
      #   print ('.  '*ind, ll)

      print (cc.list)
