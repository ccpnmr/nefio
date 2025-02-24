"""
nef - nef handling routines; a series of routines to compare/verify Nef files

Command Line Usage:
  nef for execution from the command line with a suitable script
  An example can be found in AnalysisV3/bin/nef:

      #!/usr/bin/env sh
      export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
      export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
      export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
      ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/nef.py $*

  Usage:  nef [options]

  optional arguments:
    -h, --help              show this help message
    -H, --Help              Show detailed help

    --compare               Compare Nef files: with the following options

        -i, --ignoreblockname   Ignore the blockname when comparing two Nef files
                                May be required when converting Nef files through
                                different applications.
                                May be used with -f and -b

        -f file1 file2, --files file1 file2
                                Compare two Nef files and print the results to the
                                screen

        -d dir1 dir2, --dirs dir1 dir2
                                compare Nef files common to directories
                                dir1 and dir2. Write output *.txt for each
                                file into the output directory specified below

        -o outDir, --outdir     Specify output directory for batch processing

        -s, --screen            Output batch processing to screen, default is to .txt files
                                may be used with -d

        -r, --replace           Replace existing .txt files. If false then files are
                                appended with '(n)' before the extension, where n is
                                the next available number

        -c, --create            Automatically create directories as required

        -I, --ignorecase        Ignore case when comparing items

        --same                  output similarities between Nef files
                                default is differences

        -a, --almostequal       Consider float/complex numbers to be equal if within the
                                relative tolerance

        -p, --places            Specify the number of decimal places for the relative
                                tolerance

    --verify                Verify Nef files

                            Can be used with switches: -f, -d

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

nef.py contains the following routines:
  compareNefFiles      compare two Nef files and return a comparison object

    Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
    Comparisons are made for all data structures that have the same name.
    Differences for items within a column are listed in the form:
      dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2

    dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed, in the form:
      dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                        parameter2
                                                        ...
                                                        parameterN

  A comparison object is a list of nefItems of the form:

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

  compareDataExtents    compare two DataExtent objects and return a comparison list as above
  compareDataBlocks     compare two DataBlock objects and return a comparison list as above
  compareSaveFrames     compare two SaveFrame objects and return a comparison list as above
  compareLoops          compare two Loop objects and return a comparison list as above

  compareNefFiles       compare two Nef files and return a comparison list as above
  batchCompareNefFiles  compare two directories of Nef files.
                        Nef Files common to specified directories are compared and the comparison
                        lists are written to the third directory as .txt

  printCompareList      print the comparison list to the screen
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2020"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2020-06-25 10:46:01 +0100 (Thu, June 25, 2020) $"
__version__ = "$Revision: 3.0.1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2020-05-14 16:11:22 +0000 (Thu, May 14, 2020) $"
__date__ = "$Date: 2020-06-24 13:57:48 +0000 (Wed, June 24, 2020) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
import copy
import sys


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is a fix to get the import to work when running as a standalone
# when importing into your own code, with PYTHON_PATH defined it can be safely removed

def import_parents(level=1):
    global __package__

    import sys
    from os import path
    import importlib

    # pathlib does all this a lot nicer, but don't think it's in python2.7
    top = parent = path.dirname(path.abspath(__file__))
    package = []
    for t in range(level):
        package.insert(0, os.path.basename(top))
        top = path.dirname(top)

    sys.path.append(str(top))
    try:
        sys.path.remove(str(parent))
    except ValueError:  # already removed
        pass

    __package__ = str('.'.join(package))
    importlib.import_module(__package__)


if __name__ == '__main__' and __package__ is None:
    import_parents(level=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import re
import unittest
from . import GenericStarParser, StarIo
from .SafeOpen import safeOpen
from os import listdir
from os.path import isfile, join
from enum import Enum
from collections.abc import Iterable
from collections import OrderedDict
from math import isclose
from cmath import isclose as cisclose
from ast import literal_eval

try:
    # Python 3
    from itertools import zip_longest
except:
    # python 2.7
    from itertools import izip_longest as zip_longest

EXCLUSIVEGROUP = ['compare', 'verify']
CONVERTTOSTRINGS = (int, float, complex, bool, list, tuple, dict, set, frozenset, OrderedDict, type(None))


class NEFOPTIONS(Enum):
    COMPARE = EXCLUSIVEGROUP[0]
    VERIFY = EXCLUSIVEGROUP[1]


class whichTypes(Enum):
    NONE = 0
    LEFT = 1
    RIGHT = 2
    BOTH = 3


def showMessage(msg, *args, **kwds):
    """Show a warning message
    """
    # to be subclassed as required
    print('Warning: {}'.format(msg))


def showError(msg, *args, **kwds):
    """Show an error message
    """
    # to be subclassed as required
    print('Error: {}'.format(msg))


def printOutput(*args, **kwds):
    """Output a message
    """
    # to be subclassed as required
    print(*args, **kwds)


def defineArguments():
    """Define the arguments of the program

    :return argparse instance
    """
    import argparse

    def _checkInt(value):
        try:
            intValue = int(value)
        except Exception as es:
            raise argparse.ArgumentTypeError("{} must be an int {}".format(value, es))

        if intValue < 0:
            raise argparse.ArgumentTypeError("{} must be positive int".format(value))

        return intValue

    parser = argparse.ArgumentParser(prog='compareNef',
                                     usage='%(prog)s [options]',
                                     description='Compare the contents of Nef files')

    parser.add_argument('-H', '--Help', dest='help', action='store_true', default=False, help='Show detailed help')
    parser.add_argument('-i', '--ignoreblockname', dest='ignoreBlockName', action='store_true', default=False,
                        help='Ignore the blockname when comparing two Nef files')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--files', dest='inFiles', nargs='*', default=None,
                       help='List of files for compare|verify')
    group.add_argument('-d', '--dirs', dest='batchDirs', nargs='*', default=None,
                       help='List of directories for compare|verify')
    parser.add_argument('-o', '--outdir', dest='outDir', nargs=1, default=None,
                        help='Output directory for batch compare')

    parser.add_argument('-s', '--screen', dest='screen', action='store_true', default=False, help='Output batch processing to screen')
    parser.add_argument('-r', '--replace', dest='replaceExisting', action='store_true', default=False,
                        help='Replace existing .txt files. If false, new files are appended with "(n)"')
    parser.add_argument('-c', '--create', dest='createDirs', action='store_true', default=False,
                        help='Create directories as required')
    parser.add_argument('-I', '--ignorecase', dest='ignoreCase', action='store_true', default=False,
                        help='Ignore case when comparing items')

    parser.add_argument('--same', dest='identical', action='store_true', default=False,
                        help='Output similarities between Nef files; default is differences')

    parser.add_argument('-a', '--almostequal', dest='almostEqual', action='store_true', default=True,
                        help='Consider float/complex values as equal if within tolerance')
    parser.add_argument('-p', '--places', dest='places', nargs=1, default=10, type=int, choices=range(1, 16),
                        help='Specify number of decimal places for relative tolerance')

    parser.add_argument('-m', '--maxrows', dest='maxRows', default=None, type=_checkInt,
                        help='Specify the maximum number of rows to show/print in each loop/saveframe')

    group = parser.add_mutually_exclusive_group()
    for nefItem in NEFOPTIONS:
        group.add_argument('--{}'.format(nefItem.value), dest='nefOption', action='store_const', const=nefItem,
                           help='{} Nef files'.format(nefItem.value.capitalize()))
    group.set_defaults(nefOption=NEFOPTIONS.COMPARE)

    return parser


#=========================================================================================
# compareItem
#=========================================================================================

class compareItem(object):
    """Holds the details of a compared loop/saveFrame item at a particular row/column (if required)
    """

    def __init__(self, attribute=None, row=None, column=None, thisValue=None, compareValue=None):
        self.attribute = attribute
        self.row = row
        self.column = column
        self.thisValue = thisValue
        self.compareValue = compareValue


#=========================================================================================
# nefItem
#=========================================================================================

class nefItem(object):
    """Holds the contents of a single Nef comparison
    inWhich   a flag labelling which file the item was found in
              1 = found in the first file, 2 = found on the second file, 3 = common to both
    list      a list of strings containing the comparison information
    """

    def __init__(self, cItem=None):
        self.inWhich = whichTypes.NONE
        self.strList = []
        self.objList = []
        self.compareList = []
        self.differenceList = []
        self.warningList = []
        self.errorList = []
        self.thisObj = None
        self.compareObj = None
        self._identical = False


#=========================================================================================
# _loadGeneralFile
#=========================================================================================

def _loadGeneralFile(path=None):
    """Load a file with the given pathname and return a dict of the contents

    :return entry:dict
    """
    usePath = path if path.startswith('/') else os.path.join(os.getcwd(), path)
    entry = StarIo.parseNefFile(usePath)  # 'lenient')
    printOutput(' %s' % path)
    return entry


#=========================================================================================
# printFile
#=========================================================================================

def printFile(thisFile):
    """Print a file to the screen
    """
    printOutput('~' * 80)
    printOutput(thisFile)
    for i, val in enumerate(thisFile):
        printOutput(i, thisFile[val])
        sub = thisFile[val]
        for j, valj in enumerate(sub):
            printOutput('  ', j, sub[valj])
            if j > 3:
                break

            sub2 = sub[valj]
            for k, valk in enumerate(sub2):
                loopType = sub2[valk]
                if isinstance(loopType, GenericStarParser.Loop):
                    printOutput('    ', k, 'LOOP', loopType)
                else:
                    printOutput('    ', k, loopType)

                if k > 3:
                    break


#=========================================================================================
# sizeNefList
#=========================================================================================

def sizeNefList(nefList, whichType=whichTypes.NONE):
    """List only those items that are of type whichType

    :param nefList: list to print
    :param whichType: type to print
    """
    count = 0
    if nefList is not None:
        for cCount, nefItem in enumerate(nefList):
            if nefItem.inWhich == whichType:
                count += 1
    return count


#=========================================================================================
# printWhichList
#=========================================================================================

def printWhichList(nefList, options, whichType=whichTypes.NONE):
    """List only those items that are of type whichType

    :param nefList: list to print
    :param whichType: type to print
    """

    def _remainingRows(thisList):
        """Output the remaining number of elements on the list
        """
        if maxRows is not None and maxRows < len(thisList):
            print('{} ... {} more row{}'.format(lineLeader, len(thisList) - maxRows, 's' if (len(thisList) - maxRows) > 1 else ''))

    maxRows = options.maxRows
    for cCount, nefItem in enumerate(nefList):
        if nefItem.inWhich == whichType:

            outStr = '  ' + ':'.join(obj.name for obj in nefItem.objList) + ': contains --> '
            lineTab = ' ' * len(outStr)
            lineLeader = outStr

            if isinstance(nefItem.thisObj, GenericStarParser.Loop):
                for warn in nefItem.warningList[:maxRows]:
                    print('{} {}'.format(lineLeader, warn))
                    lineLeader = lineTab
                _remainingRows(nefItem.warningList)

                for error in nefItem.errorList[:maxRows]:
                    print('{} {}'.format(lineLeader, error))
                    lineLeader = lineTab
                _remainingRows(nefItem.errorList)

                for compareObj in nefItem.compareList[:maxRows]:
                    symbol = ' == ' if nefItem._identical else ' != '

                    dataStr = '{} <Col>: {} <Row:> {} --> {} {} {}'.format(lineLeader,
                                                                           compareObj.attribute,
                                                                           compareObj.row,
                                                                           compareObj.thisValue,
                                                                           symbol,
                                                                           compareObj.compareValue)
                    printOutput(dataStr)
                    lineLeader = lineTab
                _remainingRows(nefItem.compareList)

            if isinstance(nefItem.thisObj, GenericStarParser.SaveFrame):
                for warn in nefItem.warningList[:maxRows]:
                    print('{} {}'.format(lineLeader, warn))
                    lineLeader = lineTab
                _remainingRows(nefItem.warningList)

                for error in nefItem.errorList[:maxRows]:
                    print('{} {}'.format(lineLeader, error))
                    lineLeader = lineTab
                _remainingRows(nefItem.errorList)

                for compareObj in nefItem.compareList[:maxRows]:
                    symbol = ' == ' if nefItem._identical else ' != '

                    dataStr = '{} <Value>: {} --> {} {} {}'.format(lineLeader,
                                                                   compareObj.attribute,
                                                                   compareObj.thisValue,
                                                                   symbol,
                                                                   compareObj.compareValue)
                    printOutput(dataStr)
                    lineLeader = lineTab
                _remainingRows(nefItem.compareList)

            for diffObj in nefItem.differenceList[:maxRows]:
                dataStr = '{} {}'.format(lineLeader, diffObj.attribute)
                printOutput(dataStr)
                lineLeader = lineTab
            _remainingRows(nefItem.differenceList)


#=========================================================================================
# printCompareList
#=========================================================================================

def printCompareList(nefList, inFile1, inFile2, options):
    """Print the contents of the nef compare list to the screen

    Output is in three parts:
      - items that are present only in the first file
      - items that are only in the second file
      - differences between objects that are common in both files

    :param nefList: list to print
    :param inFile1: name of the first file
    :param inFile2: name of the second file
    """

    if not isinstance(inFile1, str):
        showError('TypeError: inFile1 must be a string.')
        return
    if not isinstance(inFile2, str):
        showError('TypeError: inFile2 must be a string.')
        return

    # print the items that are only present in the first nefFile
    if sizeNefList(nefList, whichType=whichTypes.LEFT) > 0:
        printOutput('\nItems that are only present in ' + inFile1 + ':')
        printWhichList(nefList, options, whichTypes.LEFT)

    # print the items that are only present in the second nefFile
    if sizeNefList(nefList, whichType=whichTypes.RIGHT) > 0:
        printOutput('\nItems that are only present in ' + inFile2 + ':')
        printWhichList(nefList, options, whichTypes.RIGHT)

    # print the common items
    if sizeNefList(nefList, whichType=whichTypes.BOTH) > 0:
        printOutput('\nItems that are present in both files:')
        printWhichList(nefList, options, whichTypes.BOTH)


#=========================================================================================
# _filterName
#=========================================================================================

def _filterName(inName):
    """Remove rogue `n` quotes from the names.
    (This is currently only a test)

    :param inName:
    :return:
    """
    # ejb - need to remove the rogue `n` at the beginning of the name if it exists
    #       as it is passed into the namespace and gets added iteratively every save
    #       next three lines remove all occurrences of `n` from name
    regex = u'\`\d*`+?'
    return re.sub(regex, '', inName)  # substitute with ''


#=========================================================================================
# _createAttributeList
#=========================================================================================

def _createAttributeList(cItem, nefObject, inList, nefList):
    """Create a new attribute list and add to the compare list
    Currently adds one cItem with a list as the last element

    :param inList: a list of items to add to the end of cItem
    :param cItem: object containing the current tree to add to the list
    :param nefList: current list of comparisons
    :return: list of type nefItem
    """
    if len(inList) > 0:
        newItem = nefItem()
        newItem.objList = cItem.objList.copy()
        newItem.thisObj = nefObject
        newItem.inWhich = cItem.inWhich
        newItem.differenceList = [compareItem(attribute=str(item)) for item in inList]
        nefList.append(newItem)

        return newItem


#=========================================================================================
# _compareObjects
#=========================================================================================

def _compareObjects(obj1, obj2, options):
    """Compare the values of two objects
    Objects may be nested objects.
    Dicts are compared by keys
    Strings are considerd equal if lowercase values are the same of options.ignoreCase = True
    Floats/complex are considered equal if values are within the a relative tolerance as defined by
    a number of decimal places
    """

    if not type(obj1) in CONVERTTOSTRINGS:
        try:
            obj1 = literal_eval(obj1)
        except Exception as es:
            obj1 = str(obj1)
    if not type(obj2) in CONVERTTOSTRINGS:
        try:
            obj2 = literal_eval(obj2)
        except Exception as es:
            obj2 = str(obj2)

    if isinstance(obj1, Iterable) and isinstance(obj2, Iterable):
        if type(obj1) != type(obj2):
            # print('  False type >>', obj1, obj2, type(obj1), type(obj2))
            return False

        if len(obj1) != len(obj2):
            # print('  False len >>', obj1, obj2)
            return False

        # if dicts then compare keys/values
        if isinstance(obj1, dict):  # and isinstance(obj2, dict):       # shouldn't need to test both
            # compare dict values
            for d1 in obj1:
                if d1 in obj2:
                    compare = _compareObjects(obj1[d1], obj2[d1], options)
                    if not compare:
                        # print('  False bad dict item >>', obj1, obj2)
                        return False
                else:
                    # print('  False bad dict key >>', obj1, obj2)
                    return False

        elif isinstance(obj1, str):  # and isinstance(obj2, str):
            if options.ignoreCase:
                if obj1.lower() != obj2.lower():
                    # print('  False string case >>', obj1, obj2)
                    return False
            elif obj1 != obj2:
                # print('  False string >>', obj1, obj2)
                return False

        else:
            # compare values
            for s1, s2 in zip(obj1, obj2):
                compare = _compareObjects(s1, s2, options)
                if not compare:
                    # print('  False iterable item >>', obj1, obj2)
                    return False

    else:
        if isinstance(obj1, (int, float)) and isinstance(obj2, (int, float)):
            if options.almostEqual:
                if not isclose(obj1, obj2, rel_tol=pow(10, -options.places)):
                    # print('  False float tolerance >>', obj1, obj2)
                    return False
            elif obj1 != obj2:
                # print('  False float not equal >>', obj1, obj2)
                return False

        elif isinstance(obj1, complex) and isinstance(obj2, complex):
            if options.almostEqual:
                if not cisclose(obj1, obj2, rel_tol=pow(10, -options.places)):
                    # print('  False complex tolerance >>', obj1, obj2)
                    return False
            elif obj1 != obj2:
                # print('  False complex not equal >>', obj1, obj2)
                return False

        elif obj1 != obj2:
            # print('  False unequal objects >>', obj1, obj2, type(obj1), type(obj2))
            return False

    return True


#=========================================================================================
# compareLoops
#=========================================================================================

def compareLoops(loop1, loop2, options, cItem=None, nefList=None):
    """Compare two Loops

    :param loop1: first Loop object, of type GenericStarParser.Loop
    :param loop2: second Loop object, of type GenericStarParser.Loop
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
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

    cItem1 = _duplicateItem(cItem, loop1, None, inWhich=whichTypes.LEFT)
    cItem1.strList.append(loop1.name)
    cItem1.objList.append(loop1)
    _createAttributeList(cItem1, loop1, inLeft, nefList)

    cItem2 = _duplicateItem(cItem, loop2, None, inWhich=whichTypes.RIGHT)
    cItem2.strList.append(loop2.name)
    cItem2.objList.append(loop2)
    _createAttributeList(cItem2, loop2, inRight, nefList)

    if loop1.data and loop2.data:
        rowRange = max(len(loop1.data), len(loop2.data))
        nefLoopItem = None

        symbol = ' == ' if options.identical else ' != '
        # NOTE:ED - not sure whether to add this
        if len(loop1.data) != len(loop2.data):  # simple compare, same length tables - should use longest

            nefLoopItem = _createNewItem(cItem, loop1, nefList, options, inWhich=whichTypes.BOTH)
            nefLoopItem.warningList.append('<rowLength>:  {} {} {}'.format(len(loop1.data),
                                                                           symbol, len(loop2.data)))

        # carry on and compare the common table
        for compName in dSet:
            for rowIndex in range(rowRange):

                loopValue1 = loop1.data[rowIndex][compName] if rowIndex < len(loop1.data) else None
                loopValue2 = loop2.data[rowIndex][compName] if rowIndex < len(loop2.data) else None

                if _compareObjects(loopValue1, loopValue2, options) == options.identical:
                    if not nefLoopItem:
                        nefLoopItem = _createLoopItem(cItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, options, inWhich=whichTypes.BOTH)
                    else:
                        _addLoopItem(nefLoopItem, compName, loop1, loopValue1, loopValue2, nefList, rowIndex, options, inWhich=whichTypes.BOTH)

        #TODO
        # need to add a further test here, could do a diff on the tables which would pick up
        # insertions to the table - the columns would need to be reordered for this to work
        # what if there are a different number of columns?
        # also check for Mandatory items

    else:
        # NOTE:ED - not sure whether to add this
        # can't compare non-existent loopdata
        if loop1.data is None:
            newItem = _createNewItem(cItem, loop1, nefList, options, inWhich=whichTypes.LEFT)
            newItem.warningList.append('<Contains no data>')

        if loop2.data is None:
            newItem = _createNewItem(cItem, loop1, nefList, options, inWhich=whichTypes.RIGHT)
            newItem.warningList.append('<Contains no data>')

    return nefList


#=========================================================================================
# _createNewItem
#=========================================================================================

def _createNewItem(cItem, obj, nefList, options, inWhich):
    """Create a new item in the nefList to hold the current loop
    """
    # create a new item - keeping history of objects, could be loop/saveFrame/dataBock/dataExtent
    newItem = nefItem()
    newItem.strList = cItem.strList.copy()
    newItem.objList = cItem.objList.copy()
    newItem.strList.append(obj.name)
    newItem.objList.append(obj)
    newItem.thisObj = obj
    newItem.inWhich = inWhich
    newItem._identical = options.identical
    nefList.append(newItem)
    return newItem


#=========================================================================================
# _createLoopItem to the NefList or append to existing
#=========================================================================================

def _createLoopItem(cItem, compName, loop, loopValue1, loopValue2, nefList, rowIndex, options, inWhich):
    """Create a new loop item and set the compare item
    """
    # create a new item
    newItem = _createNewItem(cItem, loop, nefList, options, inWhich)
    newItem.inWhich = inWhich
    newItem.compareList = [compareItem(attribute=compName,
                                       row=rowIndex,
                                       column=compName,
                                       thisValue=loopValue1,
                                       compareValue=loopValue2)]
    newItem._identical = options.identical

    return newItem


#=========================================================================================
# _addLoopItem to the NefList or append to existing
#=========================================================================================

def _addLoopItem(cItem, compName, loop, loopValue1, loopValue2, nefList, rowIndex, options, inWhich):
    """Check the list of already added items and append to the end OR create a new item
    """
    cItem.compareList.append(compareItem(attribute=compName,
                                         row=rowIndex,
                                         column=compName,
                                         thisValue=loopValue1,
                                         compareValue=loopValue2))


#=========================================================================================
# _createSaveFrameItem to the NefList or append to existing
#=========================================================================================

def _createSaveFrameItem(cItem, compName, saveFrame, saveFrameValue1, saveFrameValue2, nefList, options, inWhich):
    """Create a new saveFrame item and set the compare item
    """
    # create a new item
    newItem = _createNewItem(cItem, saveFrame, nefList, options, inWhich)
    newItem.inWhich = inWhich
    newItem.compareList = [compareItem(attribute=compName,
                                       thisValue=saveFrameValue1,
                                       compareValue=saveFrameValue2)]
    newItem._identical = options.identical

    return newItem


#=========================================================================================
# _addSaveFrameItem to the NefList or append to existing
#=========================================================================================

def _addSaveFrameItem(cItem, compName, saveFrame, saveFrameValue1, saveFrameValue2, nefList, options, inWhich):
    """Check the list of already added items and append to the end
    """
    cItem.compareList.append(compareItem(attribute=compName,
                                         thisValue=saveFrameValue1,
                                         compareValue=saveFrameValue2))


#=========================================================================================
# compareSaveFrames
#=========================================================================================

def compareSaveFrames(saveFrame1, saveFrame2, options, cItem=None, nefList=None):
    """Compare two saveFrames, if they have the same name then check their contents

    :param saveFrame1: first SaveFrame object, of type GenericStarParser.SaveFrame
    :param saveFrame2: second SaveFrame object, of type GenericStarParser.SaveFrame
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
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

    lVSet = [str(bl) if not isinstance(saveFrame1[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame1]
    rVSet = [str(bl) if not isinstance(saveFrame2[bl], GenericStarParser.Loop) else ' ' for bl in saveFrame2]
    inVLeft = set(lVSet).difference(rVSet).difference({' '})
    dVSet = set(lVSet).intersection(rVSet).difference({' '})
    inVRight = set(rVSet).difference(lVSet).difference({' '})

    # list everything only present in the first saveFrame

    cItem1 = _duplicateItem(cItem, saveFrame1, None, whichTypes.LEFT)
    _createAttributeList(cItem1, saveFrame1, inLeft, nefList)
    _createAttributeList(cItem1, saveFrame1, inVLeft, nefList)

    # list everything only present in the second saveFrame

    cItem2 = _duplicateItem(cItem, saveFrame2, None, whichTypes.RIGHT)
    _createAttributeList(cItem2, saveFrame2, inRight, nefList)
    _createAttributeList(cItem2, saveFrame2, inVRight, nefList)

    # compare the common items

    cItem3 = _duplicateItem(cItem, saveFrame1, saveFrame2, whichTypes.BOTH)
    for compName in dSet:
        # compare the loop items of the matching saveFrames

        compareLoops(saveFrame1[compName], saveFrame2[compName], options, cItem=cItem3, nefList=nefList)

    nefLoopItem = None
    for compName in dVSet:
        if _compareObjects(saveFrame1[compName], saveFrame2[compName], options) == options.identical:
            # need to make sure these go in the same result nefItem
            # i.e. keep first item object
            if not nefLoopItem:
                nefLoopItem = _createSaveFrameItem(cItem, compName, saveFrame2, saveFrame1[compName], saveFrame2[compName], nefList,
                                                   options, inWhich=whichTypes.BOTH)
            else:
                _addSaveFrameItem(nefLoopItem, compName, saveFrame2, saveFrame1[compName], saveFrame2[compName], nefList,
                                  options, inWhich=whichTypes.BOTH)
    return nefList


#=========================================================================================
# compareDataBlocks
#=========================================================================================

def _duplicateItem(cItem, thisObj, compareObj, inWhich):
    """Create a duplicate nefItem
    """
    newItem = nefItem()
    newItem.strList = cItem.strList.copy()
    newItem.objList = cItem.objList.copy()
    newItem.strList.append(thisObj.name)
    newItem.objList.append(thisObj)
    newItem.thisObj = thisObj
    newItem.compareObj = compareObj
    newItem.inWhich = inWhich

    return newItem


def compareDataBlocks(dataBlock1, dataBlock2, options, cItem=None, nefList=None):
    """Compare two dataBlocks, if they have the same name then check their contents

    :param dataBlock1: first DataBlock object, of type GenericStarParser.DataBlock
    :param dataBlock2: second DataBlock object, of type GenericStarParser.DataBlock
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
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

    cItem1 = _duplicateItem(cItem, dataBlock1, None, whichTypes.LEFT)
    _createAttributeList(cItem1, dataBlock1, inLeft, nefList)

    # list everything only present in the second DataBlock

    cItem2 = _duplicateItem(cItem, dataBlock2, None, whichTypes.RIGHT)
    _createAttributeList(cItem2, dataBlock2, inRight, nefList)

    # compare the common items - strictly there should only be one DataBlock

    cItem3 = _duplicateItem(cItem, dataBlock1, dataBlock2, whichTypes.BOTH)
    for compName in dSet:
        compareSaveFrames(dataBlock1[compName], dataBlock2[compName], options, cItem=cItem3, nefList=nefList)

    return nefList


#=========================================================================================
# compareDataExtents
#=========================================================================================

def compareDataExtents(dataExt1, dataExt2, options, cItem=None, nefList=None):
    """Compare two dataExtents, if they have the same name then check their contents

    :param dataExt1: first DataExtent object, of type GenericStarParser.DataExtent
    :param dataExt2: second DataExtent object, of type GenericStarParser.DataExtent
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
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

    cItem1 = _duplicateItem(cItem, dataExt1, None, whichTypes.LEFT)
    _createAttributeList(cItem1, dataExt1, inLeft, nefList)

    # list everything only present in the second DataExtent

    cItem2 = _duplicateItem(cItem, dataExt2, None, whichTypes.RIGHT)
    _createAttributeList(cItem2, dataExt2, inRight, nefList)

    # compare the common items - strictly there should only be one DataExtent

    cItem3 = _duplicateItem(cItem, dataExt1, dataExt2, whichTypes.BOTH)
    for compName in dSet:
        compareDataBlocks(dataExt1[compName], dataExt2[compName], options, cItem=cItem3, nefList=nefList)

    return nefList


#=========================================================================================
# compareNefFiles
#=========================================================================================

def compareNefFiles(inFile1, inFile2, options, cItem=None, nefList=None):
    """Compare two Nef files and return comparison as a nefItem list

    :param inFile1: name of the first file
    :param inFile2: name of the second file
    :param options: nameSpace holding the commandLineArguments
    :param cItem: list of str describing differences between nefItems
    :param nefList: input of nefItems
    :return: list of type nefItem
    """
    if cItem is None:
        cItem = nefItem()
    if nefList is None:
        nefList = []

    if not os.path.isfile(inFile1):
        showError('File Error:', inFile1)
    elif not os.path.isfile(inFile2):
        showError('File Error:', inFile2)
    else:
        try:
            NefData1 = _loadGeneralFile(path=inFile1)
        except Exception as e:
            showError('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
            return None

        try:
            NefData2 = _loadGeneralFile(path=inFile2)
        except Exception as e:
            showError('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e), e)
            return None

        if options.ignoreBlockName is False:
            compareDataExtents(NefData1, NefData2, options, cItem=cItem, nefList=nefList)
        else:

            # assumes that there is only one block in a file
            # but this may change

            compList1 = [cn for cn in NefData1]
            compList2 = [cn for cn in NefData2]
            compareDataBlocks(NefData1[compList1[0]], NefData2[compList2[0]], options, cItem=cItem, nefList=nefList)

    return nefList


#=========================================================================================
# batchCompareNefFiles
#=========================================================================================

def batchCompareNefFiles(inDir1, inDir2, outDir, options):
    """Batch compare the Nef files common to the two directories
    For each file found, write the compare log to the corresponding .txt file

    :param inDir1:
    :param inDir2:
    :param outDir:
    :param options: nameSpace holding the commandLineArguments
    """
    inFileList = [f for f in listdir(inDir1) if isfile(join(inDir1, f)) and f[-4:] == '.nef']
    outFileList = [f for f in listdir(inDir2) if isfile(join(inDir2, f)) and f[-4:] == '.nef']

    if not options.screen:
        if options.createDirs is True and not os.path.exists(outDir):
            os.mkdir(outDir)
        if not (os.path.exists(outDir) and os.path.isdir(outDir)):
            showError('Error: No such directory:', str(outDir))
            return

    commonFiles = set(listdir(inDir1)) & set(listdir(inDir2))
    if not commonFiles:
        # if no files found then write message to the screen or log.tx in the out folder
        if options.screen is True:
            # strip the .nef from the end
            printOutput('inDir1: %s' % str(inDir1))
            printOutput('inDir2: %s' % str(inDir2))
            printOutput('No common files found')
        else:
            outFileName = join(outDir, 'log.txt')
            if options.replaceExisting is False:
                with safeOpen(outFileName, 'w') as (outLog, safeFileName):
                    outLog.write('inDir1: %s\n' % str(inDir1))
                    outLog.write('inDir2: %s\n' % str(inDir2))
                    outLog.write('No common files found')
            else:
                with open(outFileName, 'w') as outLog:
                    sys.stdout = outLog
                    outLog.write('inDir1: %s\n' % str(inDir1))
                    outLog.write('inDir2: %s\n' % str(inDir2))
                    outLog.write('No common files found')
        return

    for fl in inFileList:
        if fl in outFileList:

            if options.screen is True:

                # strip the .nef from the end
                outFileName = join(outDir, fl[:-4] + '.txt')
                printOutput('Batch processing %s > %s' % (fl, outFileName))

                nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)
                printCompareList(nefList, join(inDir1, fl), join(inDir2, fl), options)

            else:
                # strip the .nef from the end
                outFileName = join(outDir, fl[:-4] + '.txt')
                # keep the old output stream
                stdOriginal = sys.stdout

                if options.replaceExisting is False:

                    with safeOpen(outFileName, 'w') as (outLog, safeFileName):
                        sys.stdout = outLog
                        nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)

                        printOutput('Batch processing %s > %s' % (fl, os.path.basename(safeFileName)))
                        printOutput(join(inDir1, fl))
                        printOutput(join(inDir2, fl))
                        printCompareList(nefList, join(inDir1, fl), join(inDir2, fl), options)

                else:
                    with open(outFileName, 'w') as outLog:
                        sys.stdout = outLog
                        nefList = compareNefFiles(join(inDir1, fl), join(inDir2, fl), options)

                        printOutput('Batch processing %s > %s' % (fl, outFileName))
                        printOutput(join(inDir1, fl))
                        printOutput(join(inDir2, fl))
                        printCompareList(nefList, join(inDir1, fl), join(inDir2, fl), options)
                sys.stdout = stdOriginal


#=========================================================================================
# verifyFiles
#=========================================================================================

def verifyFile(file, options):
    """Verify a single file

    :param file:
    :param options:
    :return:
    """
    VALIDATEDICT = '/Users/ejb66/PycharmProjects/Git/NEF/specification/mmcif_nef.dic'

    from ccpn.util.nef import NefImporter as Nef

    # load the file and the validate dict
    _loader = Nef.NefImporter(errorLogging=Nef.el.NEF_STANDARD, hidePrefix=True)
    _loader.loadFile(file)
    _loader.loadValidateDictionary(VALIDATEDICT)

    # validate
    validCheck = _loader.isValid
    print(_loader.validErrorLog)

    # simple test print of saveframes
    print(validCheck)
    # names = _loader.getSaveFrameNames(returnType=Nef.NEF_RETURNALL)
    # for name in names:
    #     print(name)
    #     saveFrame = _loader.getSaveFrame(name)
    #     print(saveFrame)


def verifyFiles(inFiles, options):
    """Verify files

    :param inFiles:
    :param options: nameSpace holding the commandLineArguments
    """
    print('VERIFYING')
    badFiles = [f for f in inFiles if not (os.path.exists(f) and f[-4:] == '.nef')]
    files = [f for f in inFiles if (os.path.exists(f) and f[-4:] == '.nef')]
    if badFiles:
        print('bad')

    for file in files:
        verifyFile(file, options)


#=========================================================================================
# ProcessArguments
#=========================================================================================

def processArguments(options):
    """Process the command line arguments
    """
    if options.help:
        printOutput(_helpText)
    else:
        if options.nefOption == NEFOPTIONS.COMPARE:

            if options.inFiles is not None:

                if len(options.inFiles) == 2:

                    # compare the two files
                    inFile0 = options.inFiles[0]
                    inFile1 = options.inFiles[1]

                    printOutput()
                    printOutput('Loading Nef Files...')
                    nefList = compareNefFiles(inFile0, inFile1, options)
                    printCompareList(nefList, inFile0, inFile1, options)

                elif len(options.inFiles) < 2:
                    showError('too few files specified')
                else:
                    showError('too many files specified')

            elif options.batchDirs is not None:

                if len(options.batchDirs) == 2 and options.outDir:

                    # compare the two directories
                    inDir0 = options.batchDirs[0]
                    inDir1 = options.batchDirs[1]
                    outDir = options.outDir

                    batchCompareNefFiles(inDir0, inDir1, outDir, options)

                else:
                    if len(options.inFiles) < 2:
                        showError('too few directories specified')
                    else:
                        showError('too many directories specified')
                    if not options.outDir:
                        showError('output directory not specified')

            else:
                printOutput('Incorrect arguments, use nef -h')

        elif options.nefOption == NEFOPTIONS.VERIFY:

            # verify options here
            if options.inFiles is not None:
                verifyFiles(inFiles=options.inFiles, options=options)

            elif options.batchDirs is not None:
                pass

            else:
                printOutput('Incorrect arguments, use nef -h')

        else:

            printOutput('Incorrect arguments, use nef -h')


#=========================================================================================
# Test_Compare_Files
#=========================================================================================

class Test_compareFiles(unittest.TestCase):
    """Test the comparison of nef files and print the results
    """

    @unittest.skip
    def test_verifyFiles(self):
        """Load two files and verify
        """
        # define arguments to simulate command line
        parser = defineArguments()

        # set the two files to compare
        inFile1 = os.path.join('.', 'testdata', 'Commented_Example.nef')
        inFile2 = os.path.join('.', 'testdata', 'Commented_Example_Change.nef')

        options = parser.parse_args(('-Ic --verify -f {} {}'.format(inFile1, inFile2)).split())
        processArguments(options)

    @unittest.skip
    def test_verifySingleFile(self):
        """Load single file and verify
        """
        # define arguments to simulate command line
        parser = defineArguments()

        # set the file to compare
        inFile = '/Users/ejb66/Documents/CcpNmrData/nefTestProject.nef'

        options = parser.parse_args(('-Ic --verify -f {}'.format(inFile)).split())
        processArguments(options)

    @unittest.skip
    def test_compareSimilarFiles(self):
        """Load two files and compare
        """
        # define arguments to simulate command line
        parser = defineArguments()
        options = parser.parse_args('-Icf file1 file2'.split())

        # set the two files to compare
        inFile1 = os.path.join('.', 'testdata', 'Commented_Example.nef')
        inFile2 = os.path.join('.', 'testdata', 'Commented_Example_Change.nef')

        print('\nTEST COMPARISON')
        print('   file1 = ' + inFile1)
        print('   file2 = ' + inFile2)
        print('Loading...')

        # load and output results
        nefList = compareNefFiles(inFile1, inFile2, options)
        printCompareList(nefList, inFile1, inFile2, options)

    @unittest.skip
    def test_compareDifferentFiles(self):
        """Load two files and compare
        """
        # define arguments to simulate command line
        parser = defineArguments()
        options = parser.parse_args('-Icf file1 file2'.split())

        # set the two files to compare
        inFile1 = os.path.join('.', 'testdata', 'CCPN_1nk2_docr.nef')
        inFile2 = os.path.join('.', 'testdata', 'CCPN_2kko_docr.nef')

        # inFile1 = '/Users/ejb66/Documents/CcpNmrData/NefTestData_1_1/CCPN_Commented_Example.nef'
        inFile1 = '/Users/ejb66/Documents/CcpNmrData/nefTestProject.nef'
        inFile2 = '/Users/ejb66/Documents/CcpNmrData/nefTestProject2.nef'

        print('\nTEST COMPARISON')
        print('   file1 = ' + inFile1)
        print('   file2 = ' + inFile2)
        print('Loading...')

        # load and output results
        options.identical = False
        options.ignoreBlockName = True
        options.almostEqual = False
        options.maxRows = 8
        nefList = compareNefFiles(inFile1, inFile2, options)
        printCompareList(nefList, inFile1, inFile2, options)

    @unittest.skip
    def test_compareBatchFiles(self):
        """Compare the Nef files in two directories
        """
        # define arguments to simulate command line
        parser = defineArguments()
        options = parser.parse_args('-Icf file1 file2'.split())
        options.createDirs = True
        options.replaceExisting = False

        inDir1 = os.path.join('.', 'testdata', 'testinfolder1')
        inDir2 = os.path.join('.', 'testdata', 'testinfolder2')
        outDir = os.path.join('.', 'testdata', 'testoutfolder')

        batchCompareNefFiles(inDir1, inDir2, outDir, options)

    @unittest.skip
    def test_commandLineParser(self):
        """Test the output from the parser
        """
        # NOTE:ED - need to write some test cases here
        commandLineArguments = parser.parse_args('-Icf file1 file2 file3 -w outDir --verify'.split())

    @unittest.skip
    def test_compareObjects(self):
        """Test the compareObjects method
        """
        # set up a test dict
        testDict1 = {
            "Boolean2"  : True,
            "DictOuter" : {
                "ListSet"    : [[0, {1, 2, 3, 4, 5.0, 'more strings'}],
                                [0, 1000000.0],
                                ['Another string', 0.0]],
                "String1"    : 'this is a string',
                "nestedLists": [[0, 0],
                                [0, 1 + 2.0j],
                                [0, (1, 2, 3, 4, 5, 6), OrderedDict((
                                    ("ListSetInner", [[0, frozenset([1, 2, 3, 4, 5.0, 'more inner strings'])],
                                                      [0, 1000000.0],
                                                      {'Another inner string', 0.0},
                                                      ]),
                                    ("String1Inner", 'this is a inner string'),
                                    ("nestedListsInner", [[0, 0],
                                                          [0, 1 + 2.0j],
                                                          [0, (1, 2, 3, 4, 5, 6)]])
                                    ))
                                 ]]
                },
            "nestedDict": {
                "nestedDictItems": {
                    "floatItem": 1.23,
                    "frozen"   : frozenset({67, 78}),
                    }
                },
            "Boolean1"  : (True, None, False),
            }

        testDict2 = {
            "Boolean2"  : True,
            "DictOuter" : {
                "String1"    : 'this is a string',
                "ListSet"    : [[0, {1, 2, 3, 4, 5.00000001, 'more strings'}],
                                [0, 1000000.0],
                                ['Another string', 0.0]],
                "nestedLists": [[0, 0],
                                [0, 1 + 2.000000001j],
                                [0, (1, 2, 3, 4, 5, 6), OrderedDict((
                                    ("ListSetInner", [[0, frozenset([1, 3, 2, 4, 5.000000001, 'more inner strings'])],
                                                      [0, 1000000.0],
                                                      {'Another inner string', 0.0},
                                                      ]),
                                    ("String1Inner", 'this is a inner string'),
                                    ("nestedListsInner", [[0, 0],
                                                          [0, 1 + 2.000000001j],
                                                          [0, (1, 2, 3, 4, 5, 6)]])
                                    ))
                                 ]]
                },
            "nestedDict": {
                "nestedDictItems": {
                    "floatItem": 1.230000001,
                    "frozen"   : frozenset({78, 67}),
                    }
                },
            "Boolean1"  : (True, None, False),
            }

        options = dict()
        options.identical = False
        options.ignoreCase = True
        options.almostEqual = True
        options.maxRows = 5
        options.places = 8

        options.identical = False
        options.ignoreCase = True  # does not count for dict.keys()
        options.maxRows = 5

        options.almostEqual = True
        options.places = 5
        print('almostEqual:{} dp:{} - {}'.format(options.almostEqual, options.places, _compareObjects(testDict1, testDict2, options)))
        options.places = 10
        print('almostEqual:{} dp:{} - {}'.format(options.almostEqual, options.places, _compareObjects(testDict1, testDict2, options)))

        options.almostEqual = False
        print('almostEqual:{} dp:{} - {}'.format(options.almostEqual, options.places, _compareObjects(testDict1, testDict2, options)))
        print('almostEqual:{} dp:{} - {} (same object)'.format(options.almostEqual, options.places, _compareObjects(testDict1, testDict1, options)))


#=========================================================================================
# __main__
#=========================================================================================

_helpText = """
Command Line Usage:
  nef for execution from the command line with a suitable script
  An example can be found in AnalysisV3/bin/nef:

      #!/usr/bin/env sh
      export CCPNMR_TOP_DIR="$(dirname $(cd $(dirname "$0"); pwd))"
      export ANACONDA3=${CCPNMR_TOP_DIR}/miniconda
      export PYTHONPATH=${CCPNMR_TOP_DIR}/src/python:${CCPNMR_TOP_DIR}/src/c
      ${ANACONDA3}/bin/python ${CCPNMR_TOP_DIR}/src/python/ccpn/util/nef/nef.py $*

  Usage:  nef [options]

  optional arguments:
    -h, --help              show this help message
    -H, --Help              Show detailed help

    --compare               Compare Nef files: with the following options

        -i, --ignoreblockname   Ignore the blockname when comparing two Nef files
                                May be required when converting Nef files through
                                different applications.
                                May be used with -f and -b

        -f file1 file2, --files file1 file2
                                Compare two Nef files and print the results to the
                                screen

        -d dir1 dir2, --dirs dir1 dir2
                                compare Nef files common to directories
                                dir1 and dir2. Write output *.txt for each
                                file into the output directory specified below

        -o outDir, --outdir     Specify output directory for batch processing

        -s, --screen            Output batch processing to screen, default is to .txt files
                                may be used with -d

        -r, --replace           Replace existing .txt files. If false then files are
                                appended with '(n)' before the extension, where n is
                                the next available number

        -c, --create            Automatically create directories as required

        -I, --ignorecase        Ignore case when comparing items

        --same                  output similarities between Nef files
                                default is differences

        -a, --almostequal       Consider float/complex numbers to be equal if within the
                                relative tolerance
                                
        -p, --places            Specify the number of decimal places for the relative
                                tolerance

    --verify                Verify Nef files

                            Can be used with switches: -f, -d
                                                        
Searches through all objects: dataExtents, dataBlocks, saveFrames and Loops within the files.
Comparisons are made for all data structures that have the same name.
Differences for items within a column are listed in the form:
  dataExtent:dataBlock:saveFrame:Loop:  <Column>: columnName  <rowIndex>: row  -->  value1 != value2
  
dataExtents, dataBlocks, saveFrames, Loops, columns present only in one file are listed.
  dataExtent:dataBlock:saveFrame:Loop: contains --> parameter1
                                                    parameter2
                                                    ...
                                                    parameterN

Please see the contents of nef.py for functions available to python.
"""

if __name__ == '__main__':
    parser = defineArguments()
    commandLineArguments = parser.parse_args()

    processArguments(commandLineArguments)
