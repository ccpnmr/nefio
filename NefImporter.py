"""
NefImporter - a series of routines for reading a Nef file and examining the contents.

Module Contents
===============

NefImporter contains:

  importFile        read a Nef file into a dataStructure


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
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-07-07 16:32:41 +0100 (Fri, July 07, 2017) $"
__version__ = "$Revision: 3.0.b2 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: Ed Brooksbank $"
__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

from . import GenericStarParser
from . import StarIo

class NefDict(GenericStarParser.DataBlock):
  """Top level data block for accessing object tree"""
  # put functions in here to read the contents of the dict.
  # superclassed from DataBlock which is of type StarContainer

  def __init__(self, inDict):

    import inspect

    # keep a copy of the original dataExtent
    super(NefDict, self).__init__(name=inDict.name)
    self._nefDict = inDict

    property_names = [p for p in dir(GenericStarParser.DataBlock) if isinstance(getattr(GenericStarParser.DataBlock, p), property)]
    method_names = inspect.getmembers(inDict, predicate=inspect.ismethod)

    # add the method to point to our loaded dataExtent
    for met in method_names:
      setattr(self.__class__, met[0], met[1])

def importFile(fileName, mode='standard'):
  nefDataExtent = StarIo.parseNefFile(fileName=fileName, mode=mode)
  outDict = list(nefDataExtent.values())
  if len(outDict) > 1:
    print(
      'More than one datablock in a NEF file is not allowed.  Using the first and discarding the rest.')
  outDict = outDict[0]

  return NefDict(outDict)
