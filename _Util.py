"""Module Documentation here

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

from ccpn.core.lib import CcpnSorting
from ccpn.util.Bmrb import bmrb

# Make sure we get parse eror warnings:
bmrb.raise_parse_warnings = True


stdTagStructure = [
    {
        "category": "nef_nmr_meta_data",
        "loops": {
            "_nef_program_script": [
                "program_name",
                "script_name",
                "script",
                "ccpn_special_par"
            ],
            "_nef_related_entries": [
                "database_name",
                "database_accession_code"
            ],
            "_nef_run_history": [
                "run_ordinal",
                "program_name",
                "program_version",
                "script_name",
                "script"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "format_name",
            "format_version",
            "program_name",
            "program_version",
            "creation_date",
            "coordinate_file_name",
            "_nef_related_entries",
            "_nef_program_script",
            "_nef_run_history"
        ]
    },
    {
        "category": "nef_molecular_system",
        "loops": {
            "_nef_covalent_links": [
                "chain_code_1",
                "sequence_code_1",
                "residue_type_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_type_2",
                "atom_name_2"
            ],
            "_nef_sequence": [
                "chain_code",
                "sequence_code",
                "residue_type",
                "residue_variant",
                "linking",
                "cross_linking"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "_nef_sequence",
            "_nef_covalent_links"
        ]
    },
    {
        "category": "nef_chemical_shift_list",
        "loops": {
            "_nef_chemical_shift": [
                "chain_code",
                "sequence_code",
                "residue_type",
                "atom_name",
                "value",
                "value_uncertainty"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "atom_chemical_shift_units",
            "_nef_chemical_shift"
        ]
    },
    {
        "category": "nef_dihedral_restraint_list",
        "loops": {
            "_nef_dihedral_restraint": [
                "restraint_id",
                "restraint_combination_id",
                "chain_code_1",
                "sequence_code_1",
                "residue_type_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_type_2",
                "atom_name_2",
                "chain_code_3",
                "sequence_code_3",
                "residue_type_3",
                "atom_name_3",
                "chain_code_4",
                "sequence_code_4",
                "residue_type_4",
                "atom_name_4",
                "weight",
                "target_value",
                "target_value_uncertainty",
                "lower_linear_limit",
                "lower_limit",
                "upper_limit",
                "upper_linear_limit"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "potential_type",
            "_nef_dihedral_restraint"
        ]
    },
    {
        "category": "nef_distance_restraint_list",
        "loops": {
            "_nef_distance_restraint": [
                "restraint_id",
                "restraint_combination_id",
                "chain_code_1",
                "sequence_code_1",
                "residue_type_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_type_2",
                "atom_name_2",
                "weight",
                "target_value",
                "target_value_uncertainty",
                "lower_linear_limit",
                "lower_limit",
                "upper_limit",
                "upper_linear_limit"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "potential_type",
            "_nef_distance_restraint"
        ]
    },
    {
        "category": "nef_rdc_restraint_list",
        "loops": {
            "_nef_rdc_restraint": [
                "restraint_id",
                "restraint_combination_id",
                "chain_code_1",
                "sequence_code_1",
                "residue_type_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_type_2",
                "atom_name_2",
                "weight",
                "target_value",
                "target_value_uncertainty",
                "lower_linear_limit",
                "lower_limit",
                "upper_limit",
                "upper_linear_limit"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "potential_type",
            "tensor_magnitude",
            "tensor_rhombicity",
            "tensor_chain_code",
            "tensor_sequence_code",
            "tensor_residue_type",
            "_nef_rdc_restraint"
        ]
    },
    {
        "category": "nef_nmr_spectrum",
        "loops": {
            "_nef_peak": [
                "peak_id",
                "volume",
                "volume_uncertainty",
                "height",
                "height_uncertainty",
                "position_1",
                "position_uncertainty_1",
                "position_2",
                "position_uncertainty_2",
                "position_3",
                "position_uncertainty_3",
                "chain_code_1",
                "sequence_code_1",
                "residue_type_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_type_2",
                "atom_name_2",
                "chain_code_3",
                "sequence_code_3",
                "residue_type_3",
                "atom_name_3"
            ],
            "_nef_spectrum_dimension": [
                "dimension_id",
                "axis_unit",
                "axis_code",
                "spectrometer_frequency",
                "spectral_width",
                "value_first_point",
                "folding",
                "absolute_peak_positions",
                "is_acquisition"
            ],
            "_nef_spectrum_dimension_transfer": [
                "dimension_1",
                "dimension_2",
                "transfer_type",
                "is_indirect"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "num_dimensions",
            "chemical_shift_list",
            "experiment_classification",
            "experiment_type",
            "_nef_spectrum_dimension",
            "_nef_spectrum_dimension_transfer",
            "_nef_peak"
        ]
    },
    {
        "category": "nef_peak_restraint_links",
        "loops": {
            "_nef_peak_restraint_link": [
                "nmr_spectrum_id",
                "peak_id",
                "restraint_list_id",
                "restraint_id"
            ]
        },
        "tags": [
            "sf_category",
            "sf_framecode",
            "_nef_peak_restraint_link"
        ]
    }
]


# def removeCommon(bmrbData1, bmrbData2):
#   """Remove identical elements in two bmrbData blocks to simplify comparison"""
#   dataBlocks = [bmrbData1, bmrbData2]
#   for frameCode in dataBlocks[0]:
#     if frameCode in dataBlocks[1]:
#       # Identical frames
#       frames = [dataBlocks[0][frameCode], dataBlocks[1][frameCode]]
#
#       # Check identical tags
#       for i0 in (0,1):
#         i1 = 1 - i0
#         for tag in frames[i0]:
#           if tag in frames[i1]:
#             # identical tags
#             if frames[i0].get(tag) == frames[i1].get(tag):
#               for frame in frames:
#                 if tag in frame:
#                   # NBNB Done like this because we might have None v. Missing
#                   del frame[tag]
#           elif isinstance(frames[i0][tag], list):
#             frames[i0][tag] = 'LOOP ONLY HERE. %s rows' % len(frames[i0][tag])
#
#       for tag in frames[0]:
#         if tag in frames[1] and isinstance(frames[0][tag], list):
#           # Compare loop contents
#           loops = [x[tag] for x in frames]
#           if loops[0] and loops[1]:
#
#             # get non-empty columns
#             ziploops = [[zip(*(y for y in x))] for x in loops]
#             ll = [x.keys() for x in loops]
#             keyLists = [[],[]]
#             for ii in 0,1:
#               for jj,column in ll[ii]:
#                 if any(*ziploops[ii][jj]):
#                   keyLists[ii].append(column)
#             commonColumns = keyLists[0] + [x for x in keyLists[1] if x not in keyLists[0]]
#             del ziploops
#
#             newLoops = []
#             for ii in 0,1:
#               loop = loops[ii]
#               newLoops.append(sorted([tuple(dd.get(tag) for tag in commonColumns)
#                                       for dd in loop]))
#             loop0 = newLoops[0]
#             loop1 = newLoops[1]
#             for ii in range(len(loop0)-1, -1, -1):
#               for jj in range(len(loop1)-1, -1, -1):
#                 if loop0[ii] == loop1[jj]:
#                   del loop0[ii]
#                 del loop1[jj]
#
#
#
#
#
#
#
#
#     else:
#       dataBlocks[0][frameCode] = 'FRAME ONLY HERE'
#
#   for frameCode in dataBlocks[1]:
#     if frameCode not in dataBlocks[0]:
#         dataBlocks[0][frameCode] = 'FRAME ONLY HERE'


def _getTagStructure(nefFile):
  """Quick hack to get tag structure and order from a file"""
  OrderedDict = collections.OrderedDict
  entry = bmrb.entry.fromFile(nefFile)
  result = []
  for frame in entry.frame_list:
    tagDict = OrderedDict(frame.tags)
    category = tagDict['sf_category']
    if not category.startswith('nef'):
      continue
    for categoryDict in result:
      if categoryDict['category'] == category:
        break
    else:
      categoryDict  = {'category':category, 'tags':[], 'loops':{}}
      result.append(categoryDict)

    # get tags
    tags0 = list(tagDict.keys())
    ll = [tags0, categoryDict['tags']]
    if len(ll[0]) < len(ll[1]):
      ll.reverse()
    categoryDict['tags'] = tags = ll[0] + [x for x in ll[1] if x not in ll[0]]

    # get loops
    for loop in frame.loops:
      prefix = loop.category
      if prefix in tags:
        columns = loop.columns
        ll = [columns, categoryDict['loops'][prefix]]
        if len(ll[0]) < len(ll[1]):
          ll.reverse()
        categoryDict['loops'][prefix] = ll[0] + [x for x in ll[1] if x not in ll[0]]
      else:
        tags.append(prefix)
        categoryDict['loops'][prefix] = loop.columns

  for dd in result:

    tags = dd['tags']
    loops = dd['loops']
    dd['tags']= [x for x in tags if x not in loops] + [x for x in tags if x in loops]

  #
  return result



def bmrbEntry2Dict(entry:bmrb.entry)-> dict:
  """Convert Bmrb entry to nested OrderedDict data structure:
  - an entry is an OrderedDict of saveframe OrderedDicts
  - a saveframe OrderedDict contains tag:value and loopPrefix:loop
  - a loop is a list of namedtuples with row data named after column names"""

  OrderedDict = collections.OrderedDict
  result = OrderedDict()
  #
  for frame in entry.frame_list:
    # frameDict = OrderedDict((tt[0], None if tt[1] == '.' else tt[1]) for tt in frame.tags)
    frameDict = OrderedDict(frame.tags)
    result[frameDict['sf_framecode']] = frameDict
    for loop in frame.loops:
      prefix = loop.category
      columns = loop.columns
      frameDict[prefix] = [OrderedDict(zip(columns, row))
                           for row in loop.data]
  #
  return result

def regulariseNefFile(inPath, outPath=None):
  """Re-export file with regular frame, tag. and column order
  and saveframes sorted by category, frameCode

  Export file name has '_r' appended to main file name"""

  raise NotImplementedError(("Still to be upgraded"))

  nefDir, nefFile = os.path.split(inPath)
  nefName, nefExt = os.path.splitext(nefFile)
  if outPath is None:
    outPath = os.path.join(nefDir, '_r'.join((nefName, nefExt)))

  entryIn = bmrb.entry.fromFile(inPath)
  entryOut = regulariseEntry(entryIn)

  open(outPath, 'w').write(str(entryOut))

def regulariseEntry(entryIn):
  """Re-export file with regular frame, tag. and column order
  and saveframes sorted by category, frameCode

  Export file name has '_r' appended to main file name"""

  dictIn = bmrbEntry2Dict(entryIn)
  categories = set()

  result = bmrb.entry.fromScratch(entryIn.bmrb_id)

  for frameTags in stdTagStructure:
    category = frameTags['category']
    tags = frameTags['tags']
    loopTags = frameTags['loops']
    categories.add(category)

    for frameCode, frameData in sorted(dictIn.items()):

      if frameData['sf_category'] != category:
        continue

      saveframe = bmrb.saveframe.fromScratch(saveframe_name=frameCode,
                                             tag_prefix=category)
      result.addSaveframe(saveframe)

      # add std tags
      saveframe.addTags([(tag, frameData.get(tag, '.')) for tag in tags
                         if tag != 'sf_cat'
                                   'categories = set{}egory' and tag not in loopTags])

      # add additional tags
      saveframe.addTags([(tag, val) for tag,val in frameData.items()
                         if tag != 'sf_category' and tag not in tags
                         and isinstance(val,str)])

      # add std loops
      for prefix in tags:
        if prefix in loopTags:
          data0 = frameData.get(prefix)
          if data0:
            useTags = (loopTags[prefix] +
                       list(sorted(x for x in data0[0] if x not in loopTags[prefix])))
            # data = []
            # for dd in data0:
            #   row = [dd.get(x, '.') for x in useTags]
            #   # data.append([(int(ss) if isinstance(ss, str) and all(x.isdigit() for x in ss) else ss) for ss in row])
            #   data.append(row)
            #
            obj = bmrb.loop.fromScratch(category=prefix)
            saveframe.addLoop(obj)
            for tag in useTags:
              obj.addColumn(tag)
            # for row in sorted(data, key=integerStringSortKey):
            for row in sorted(([dd.get(x, '.') for x in useTags] for dd in data0),
                              key=CcpnSorting.universalSortKey):
              obj.addData(row)

      # Add additional loops
      for prefix, data0 in frameData.items():
        if prefix not in tags and not isinstance(data0, str):

          if data0:
            useTags = list(data0[0].keys())
            # data = []
            # for dd in data0:
            #   data.append([(int(ss) if isinstance(ss, str) and all(x.isdigit() for x in ss) else ss) for ss in dd.values()])
            #
            obj = bmrb.loop.fromScratch(category=prefix)
            saveframe.addLoop(obj)
            for tag in useTags:
              obj.addColumn(tag)
            for row in sorted(data0, key=CcpnSorting.universalSortKey):
              obj.addData(row)

    # write out additional frames
    for frameCode, frameData in sorted(dictIn.items()):

      category = frameData['sf_category']
      if category in categories:
        continue

      saveframe = bmrb.saveframe.fromScratch(saveframe_name=category,
                                             tag_prefix=category)

      # add additional tags
      saveframe.addTags([(tag, val) for tag,val in frameData.items()
                         if tag != 'sf_category'
                         and isinstance(val,str)])

      # Add additional loops
      for prefix, data0 in frameData.items():
        if isinstance(data0, list):

          if data0:
            useTags = list(data0[0].keys())
            # data = []
            # for dd in data0:
            #   data.append([(int(ss) if isinstance(ss, str) and all(x.isdigit() for x in ss) else ss) for ss in dd.values()])
            #
            obj = bmrb.loop.fromScratch(category=prefix)
            saveframe.addLoop(obj)
            for tag in useTags:
              obj.addColumn(tag)
            for row in sorted(data0, key=CcpnSorting.universalSortKey):
              try:
                obj.addData(row)
              except:
                print ("     %s" % row)
                raise
  return result

def compareNefFiles(nefPath1, nefPath2):
  entry1 = bmrb.entry.fromFile(nefPath1)
  entry2 = bmrb.entry.fromFile(nefPath2)
  print(compareEntries(entry1, entry2))

def compareEntries(entry1, entry2):

  data1 = bmrbEntry2Dict(regulariseEntry(entry1))
  data2 = bmrbEntry2Dict(regulariseEntry(entry2))
  data = [data1, data2]
  files = [x.bmrb_id for x in (entry1, entry2)]
  if files[0] == files[1]:
    files = [files[0]+'_1', files[0] + '_2']

  diffs = []
  oneline = '%s\t\t%s\t\t%s'

  val1, val2 = entry1.bmrb_id, entry2.bmrb_id
  if val1 != val2:
    diffs.append(oneline % ('bmrb_id', val1, val2))

  if diffs and diffs[-1] != '\n':
    diffs.append('\n')

  for ii in (0,1):
    jj = 1 - ii
    for name in data[ii]:
      if name not in data[jj]:
        diffs.append("Saveframe %s only in %s" %(name, files[ii]))
    if diffs and diffs[-1] != '\n':
      diffs.append('\n')

  for name in data[0]:
    if name in data[1]:
      compareSaveframes(name, diffs, [data[0][name], data[1][name]], files)


  if diffs:
    return "Files differ\n" + '\n'.join(diffs)
  else:
    return "Files are identical\n"

def compareSaveframes(prefix, diffs, data, files):

  oneline = "Difference: %s\t%s\t%s\t%s"

  for ii in (0,1):
    jj = 1 - ii
    for name in data[ii]:
      if name not in data[jj]:
        diffs.append("Saveframe,tag %s %s only in %s" % (prefix, name, files[ii]))

  if diffs and diffs[-1] != '\n':
    diffs.append('\n')

  for name in data[0]:
    if name in data[1]:
      val1, val2 = data[0][name], data[1][name]
      if isinstance(val1, list):
        # This is a loop
        compareLoops(prefix, name, diffs, [val1, val2], files)
      elif val1 != val2:
        diffs.append(oneline % (prefix, name, val1, val2))


def compareLoops(frameCode, prefix, diffs, data, files):

  maxDifferences = 10

  for ii in (0,1):
    jj = 1 - ii
    extraCols = [x for x in data[ii][0] if x not in  data[jj][0]]
    if extraCols:
      diffs.append("%s\t %s\t Extra columns %s in %s\n"
                   % (frameCode, prefix, extraCols, files[ii]))

  columns = [x for x in data[0][0] if x in data[1][0]]

  data = [[tuple(dd[x] for x in columns) for dd in data[ii]] for ii in (0,1)]

  if any(data):
    sets = [set(ll) for ll in data]
    extraRows = [[x for x in data[ii] if x not in sets[1-ii]] for ii in (0,1)]
    if any(extraRows):
      commonCount = len(data[0]) - len(extraRows[0])
      if not commonCount and len(data[0]) == len(data[1]):
        print("All %s rows differ" % commonCount)
      else:
        print("Loops differ: common rows: %s, %s only:  %s; %s only: %s" %
              (len(data[0]) - len(extraRows[0]), files[0], len(extraRows[0]),
               files[1], len(extraRows[1])))
      if max(len(x) for x in extraRows) <= maxDifferences:
        for ii, rows in enumerate(extraRows):
          ll = ["First %s extra rows for %s %s in %s" %
                (maxDifferences, frameCode, prefix, files[ii])]
          ll.append('\t'.join(columns))
          for row in rows[:maxDifferences]:
            ll.append('\t'.join(str(x) for x in row))
          ll.append('')
          diffs.append('\n'.join(ll))


if __name__ == '__main__':

  # Export to JSON
  # bmrbData = bmrbEntry2Data1(bmrb.entry.fromFile(sys.argv[1]))

  # data = _getTagStructure(sys.argv[1])
  # print(json.dumps(data, sort_keys=True, indent=4))

  # # regularise file
  inPath = sys.argv[1]
  if len(sys.argv) > 2:
    outPath = sys.argv[2]
  else:
    outPath = None
  if os.path.isdir(inPath):
    for fileName in os.listdir(inPath):
      print("Regularising %s" % fileName)
      files = list(os.path.join(x, fileName) for x in (inPath, outPath))
      try:
        regulariseNefFile(*files)
      except:
        print ("File %s raised an ERROR" % files)
  else:
      regulariseNefFile(inPath, outPath)


  # Compare files
  # inDir = sys.argv[1]
  # fileNames = os.listdir(inDir)
  #
  # for fileName in fileNames:
  #   tt = fileName.split('-', 1)
  #   if len(tt) == 2:
  #     originalName = tt[1]
  #     files = tuple(os.path.join(inDir, x) for x in (fileName, originalName))
  #     print("\n\n\n\nCOMPARING %s to %s :" % (fileName, originalName))
  #     try:
  #       compareNefFiles(*files)
  #     except:
  #       print("Comparison raised an ERROR")