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
  # test_Compare_Files
  #=========================================================================================

  def test_Compare_Files(self):
    """
    Load two files and compare
    """
    print ('Loading...')
    print ('  # Commented_Example.nef')
    file1 = self._loadGeneralFile(path='Commented_Example.nef')
    print ('  # CCPN_H1GI.nef')
    file2 = self._loadGeneralFile(path='CCPN_H1GI.nef')

    textFile1 = file1.toString(indent=' ')
    textFile2 = file2.toString(indent=' ')

    print ('~'*80)

    self.printFile(file1)
    self.printFile(file2)

    # compare name of root
    # compare names of saveframes
    # compare names of loops within saveframes

    print ('='*80)

    if file1.name == file2.name:
      print ('DE:  Same DataExtent name')
    else:
      print('DE:  Different DataExtent name')

    for i, val in enumerate(file1):
      left = file1[val]
      for j, valj in enumerate(file1):
        right = file1[valj]
        if left.name == right.name:
          print ('DB:    ', left.name, 'in both - do more comparing')

          for x, valx in enumerate(left):
            lleft = left[valx]
            for y, valy in enumerate(right):
              rright = right[valy]

              if lleft.name == rright.name:
                print ('SF:      ', lleft.name, 'in both - do more comparing')

        else:
          print ('DB:    ', right.name, 'not in file2')