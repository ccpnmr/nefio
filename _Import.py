"""NEF file import

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date$"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon P Skinner, Geerten W Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/license"
              "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/license"
                " or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author$"
__date__ = "$Date$"
__version__ = "$Revision$"

#=========================================================================================
# Start of code
#=========================================================================================


# NBNB FILE OBSOLETE



import collections
import os
import sys

# from ccpn.util import Common as commonUtil
from . import _Export
from ccpn.util.Bmrb import bmrb

def bmrbEntry2Data1(entry:bmrb.entry)-> dict:
  """Convert Bmrb entry to nested OrderedDict data structure:
  - an entry is an OrderedDict of saveframe OrderedDicts
  - a saveframe OrderedDict contains tag:value and loopPrefix:loop
  - a loop is a list of namedtuples with row data named after column names"""


  result = {}
  #
  for frame in entry.frame_list:
    frameDict = dict((tt[0], None if tt[1] == '.' else tt[1])
                            for tt in frame.tags)
    result[frameDict['sf_framecode']] = frameDict
    for loop in frame.loops:
      # Struct = collections.namedtuple(underScore2camelCase(loopPrefix[1:]),
      #                                 [underScore2camelCase(x) for x in loop.columns])
      prefix = loop.category
      columns = loop.columns
      frameDict[prefix] = [dict(zip(columns, (None if x == '.' else x for x in row)))
                           for row in loop.data]
  #
  return result








def bmrbEntry2Data(entry:bmrb.entry)-> collections.OrderedDict:
  """Convert Bmrb entry to nested OrderedDict data structure:
  - an entry is an OrderedDict of saveframe OrderedDicts
  - a saveframe OrderedDict contains tag:value and loopPrefix:loop
  - a loop is a list of namedtuples with row data named after column names"""

  OrderedDict = collections.OrderedDict

  result = OrderedDict()
  #
  for frame in entry.frame_list:
    frameDict = OrderedDict((underScore2camelCase(tt[0]), None if tt[1] == '.' else tt[1])
                            for tt in frame.tags)
    result[frameDict['sfFramecode']] = frameDict
    for loop in frame.loops:
      # Struct = collections.namedtuple(underScore2camelCase(loopPrefix[1:]),
      #                                 [underScore2camelCase(x) for x in loop.columns])
      prefix = underScore2camelCase(loop.category[1:])
      columns = [underScore2camelCase(x) for x in loop.columns]
      frameDict[prefix] = [OrderedDict(zip(columns, (None if x == '.' else x for x in row)))
                           for row in loop.data]
  #
  return result


# def bmrbEntry2Data(entry:bmrb.entry)-> collections.OrderedDict:
#   """Convert Bmrb entry to nested OrderedDict data structure:
#   - an entry is an OrderedDict of saveframe OrderedDicts
#   - a saveframe OrderedDict contains tag:value and loopPrefix:loop
#   - a loop is a list of namedtuples with row data named after column names"""
#
#   OrderedDict = collections.OrderedDict
#
#   result = OrderedDict()
#   #
#   for frame in entry.frame_list:
#     frameDict = OrderedDict((underScore2camelCase(tt[0]), None if tt[1] == '.' else tt[1])
#                             for tt in frame.tags)
#     result[frameDict['sfFramecode']] = frameDict
#     for loop in frame.loops:
#       # Struct = collections.namedtuple(underScore2camelCase(loopPrefix[1:]),
#       #                                 [underScore2camelCase(x) for x in loop.columns])
#       prefix = underScore2camelCase(loop.category[1:])
#
#       columns = loop.columns
#       listTags = [x for x in columns if x.endswith('_1')]
#       tag2List = {}
#       listLengths = {}
#       for listTag in listTags:
#         tagRoot = listTag[:-1]
#         newTag = underScore2camelCase(listTag[:-2])
#
#         count = 0
#         while True:
#           count += 1
#           tag = tagRoot + str(count)
#           if tag in columns:
#             tag2List[tag] = (newTag, count - 1)
#           else:
#             listLengths[newTag] = count
#             break
#
#
#
#       columns = [underScore2camelCase(x) for x in loop.columns]
#       frameDict[prefix] = [OrderedDict(zip(columns, (None if x == '.' else x for x in row)))
#                            for row in loop.data]
#   #
#   return result

def underScore2camelCase(string):
  """Converts under_score separated names to camelCase"""

  parts = string.lower().split('_')
  for ii in range(1,len(parts)):
    ss = parts[ii]
    parts[ii] = ss[0].upper() + ss[1:]
  return ''.join(parts)

def importNefFile(project, nefPath):
  """Import nef file nefPath into (wrapper) project
  returning NmrCalcStore describing the imported data"""

  print("Reading BMRB entry %s" % nefPath)
  entry = bmrb.entry.fromFile(nefPath)
  entry.printTree()
  bmrbData = bmrbEntry2Data(entry)
  del entry

  nmrProject = project._wrappedData
  nmrCalcStore = (nmrProject.findFirstNmrCalcStore(name='nefIO', nmrProject=nmrProject) or
                  nmrProject.root.newNmrCalcStore(name='nefIO', nmrProject=nmrProject))
  nmrCalcRun = nmrCalcStore.newRun(details="NEF I/O instance", status='active')

  frameCount = {}
  for frameCode,frameDict in bmrbData.items():
    # Read frames in order
    category = frameDict.get("sfCategory")

    if category == 'nef_nmr_meta_data':
      if category in frameCount:
        raise ValueError("NEF file has more than %s %s frames"
                         % (frameCount[category], category))
      else:
        readNmrMetaData(frameDict, project, nmrCalcRun)
        frameCount[category] = 1

    elif category == 'nef_molecular_system':
      if category in frameCount:
        raise ValueError("NEF file has more than %s %s frames"
                         % (frameCount[category], category))
      else:
        readMolecularSystem(frameDict, project, nmrCalcRun)
        frameCount[category] = 1

    elif category == 'nef_chemical_shift_list':
      count = frameCount.get(category, 0)
      readChemicalShiftList(frameDict, project, nmrCalcRun)
      frameCount[category] = count + 11

    elif category == 'nef_distance_restraint_list':
      count = frameCount.get(category, 0)
      readDistanceRestraintList(frameDict, project, nmrCalcRun)
      frameCount[category] = count + 1

    elif category == 'nef_dihedral_restraint_list':
      count = frameCount.get(category, 0)
      readDihedralRestraintList(frameDict, project, nmrCalcRun)
      frameCount[category] = count + 11

    elif category == 'nef_rdc_restraint_list':
      count = frameCount.get(category, 0)
      readRdcRestraintList(frameDict, project, nmrCalcRun)
      frameCount[category] = count + 1

    elif category == 'nef_nmr_spectrumt':
      count = frameCount.get(category, 0)
      readNmrSpectrum(frameDict, project, nmrCalcRun)
      frameCount[category] = count + 1

    elif category == 'nef_peak_restraint_links':
      if category in frameCount:
        raise ValueError("NEF file has more than %s %s frames"
                         % (frameCount[category], category))
      else:
        readPeakRestraintLinks(frameDict, project, nmrCalcRun)
        frameCount[category] = 1

    else:
      readAdditionalFrame(frameDict, project, nmrCalcRun)


def readNmrMetaData(frameDict, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet"
         % frameDict.get("sf_category"))

def readMolecularSystem(frameDict, project, nmrCalcRun):
  """Read Molecular System from frameDict structure
  """
  # Other linkings, notably None, are treated as 'single' if outside a chain, 'middle within it

  convertLinking = {
    'start':'start',
    'end':'end',
    'single':'none',
  }

  ccpnMolSystem = project._wrappedData.molSystem
  sequence = []
  chainCodes = []
  for dd in frameDict.get('nefSequence'):
    chainCode = dd['chainCode']
    if chainCode not in chainCodes:
      chainCodes.append(chainCode)

    # NBNB TBD descriptor is not set
    sequence.append((
      chainCode,
      dd['sequenceCode'],
      dd['residueType'],
      convertLinking.get(dd['linking']),
      None
    ))
  chains = project._makeChains(sequence)
  for ii, chain in chains:
    nmrCalcRun.newMolResidueData(name=chain.code, code=chainCodes[ii], ioRole='output',
                                   molSystemCode=ccpnMolSystem.code,
                                   chainCodes=(chain.code,), )




def readChemicalShiftList(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readDistanceRestraintList(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readDihedralRestraintList(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readRdcRestraintList(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readNmrSpectrum(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readPeakRestraintLinks(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

def readAdditionalFrame(frame, project, nmrCalcRun):
  print ("WARNING: reading of %s saveframes not implemented yet" % frame.get("sf_category"))

if __name__ == '__main__':

  from ccpn.core._implementation import Io as ccpnIo

  if len(sys.argv) >= 2:

    from ccpnmodel.ccpncore.lib.Io import Api as apiIo
    from ccpn.core.Project import Project

    # set up input
    junk, readPath = sys.argv[:2]
    outputDir = os.getcwd()
    nefDir, nefFile = os.path.split(readPath)
    nefName, junk = os.path.splitext(nefFile)

    if len(sys.argv) == 2:
      pp = ccpnIo.newProject(nefName, path=outputDir)

    else:
      projectDir = sys.argv[3]
      # NBNB TODO change this - no direct apiIo.load
      ccpnProject = apiIo.loadProject(projectDir)
      nmrProj = ccpnProject.findFirstNmrProject()
      _Export.prepareNmrProject(nmrProj)

      pp = Project(nmrProj)

    newNmrCalcStore = importNefFile(pp, readPath)

  else:
    print ("Error. Parameters are: ccpnProjectDirectory outputDirectory [constraintStoreSerial] ")

  # Export to JSON
  # NBNB TBD continue here
  bmrbData = bmrbEntry2Data1(bmrb.entry.fromFile(sys.argv[1]))
