"""NEF (Nmr Exchange Format) exporter code
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


import datetime
import os
import random
import sys
from typing import Sequence

from ccpn.core.lib import CcpnSorting
from ccpn.core.lib import RestraintLib
from ccpn.util import Path
from ccpn.core.lib import Pid
from ccpn.util.Bmrb import bmrb
from ccpnmodel.ccpncore.lib.spectrum import Spectrum as libSpectrum
from ccpnmodel.ccpncore.memops import Version

nefExtension = 'nef'
# Max value used for random integer. Set to be expressible as a signed 32-bit integer.
maxRandomInt = 2000000000

formatVersion = 0.9

# potential tags to use depending on potential type - underscore version
# NBNB TBD move to more appropriate location?
# NB target+value and target_value_error are not included as they are used each time
nefTagsByPotentialType = {
  'log-normal' : (),
  'parabolic' : (),
  'square-well-parabolic' : ( 'lower_limit', 'upper_limit'),
  'upper-bound-parabolic' : ('upper_limit',),
  'lower-bound-parabolic' : ('lower_limit',),
  'square-well-parabolic-linear' : ('lower_linear_limit', 'lower_limit','upper_limit',
                                    'upper_linear_limit'),
  'upper_bound-parabolic-linear' : ('upper_limit','upper_linear_limit'),
  'lower_bound-parabolic-linear' : ('lower_linear_limit', 'lower_limit'),
  'unknown' : ('lower_linear_limit', 'lower_limit','upper_limit', 'upper_linear_limit'),
}

nef2CcpnTags = {
  'lower_limit' : 'lowerLimit',
  'upper_limit' : 'upperLimit',
  'lower_linear_limit' : 'additionalLowerLimit',
  'upper_linear_limit' : 'additionalUpperLimit',
}


# NBNB TBD we might want an exportNmrCalc2Nef wrapper as well

def exportRestraintStore(dataSet, dataName=None, directory=None,
                         forceFirstShiftList=False):
  """Export dataSet and associated data to NEF

  Nb forceFirstShiftList is a hack for cases where shiftlists are not properly set"""

  # nmrProject
  project = dataSet._project

  restraintLists = dataSet.restraintLists

  # PeakLists
  # Now objects are not hashable we need to do all this to remove duplicates
  # As there will be few peaklists, sorting is probably not warranted
  peakLists = []
  for peakList in (z.peakList for x in restraintLists
                   for y in x.restraints
                   for z in y.peaks):
    if peakList not in peakLists:
      peakLists.append(peakList)
  peakLists.sort()

  # NBNB TBD temporary hack to get right result for CASD entries:
  if not peakLists:
    peakLists = [y for x in project.spectra for y in x.peakLists]

  shiftLists = project.chemicalShiftLists
  if len(shiftLists) == 1 or forceFirstShiftList or not peakLists:
    shiftList = shiftLists[0]
  else:
    shiftList = None

  if dataName:
    ll = dataName.rsplit('.',1)
    if len(ll) == 2 and ll[1] == nefExtension:
      dataName = ll[0]
  else:
    dataName = project.name.translate(Pid.remapSeparators)
    dataName = '%s-%s' % (dataName, dataSet.serial)

  entry = _makeStarEntry(project, dataName, project.chains, peakLists,
                        restraintLists, shiftList=shiftList)

  # Export object tree
  if not directory:
    # set directory
    directory = os.getcwd()
  if not os.path.exists(directory):
    os.makedirs(directory)

  # write file
  filePath = os.path.join(directory, dataName)
  filePath = '.'.join((filePath, nefExtension))
  open(filePath, 'w').write(entry.nefString())


def _makeStarEntry(project, dataName, chains=(), peakLists=(), restraintLists=(),
                  shiftList=None, programName:str='CCPN', programVersion:str=None,
                      coordinateFileName:str=None, relatedEntries:Sequence=(),
                      programScripts:Sequence=(), programHistory:Sequence=(), uuid=None):
  """Make Bmrb Sans entry for export.
   shift lists are taken from peaklists, otherwise from shiftList parameter.
   siftList is used for PeakLists that have none set."""

  peakLists = list(sorted(peakLists))
  restraintLists = list(sorted(restraintLists))
  chains = list(sorted(chains))

  # Set up parameters and sanity check

  # Shift list check
  aSet = set()
  if shiftList is not None:
    aSet.add(shiftList)
  for peakList in peakLists:
    aSet.add(peakList.chemicalShiftList)
  shiftLists = [x for x in project.chemicalShiftLists if x in aSet]
  if not shiftLists:
    ll = project.chemicalShiftLists
    if len(ll) == 1:
      shiftLists = ll
      shiftList = shiftLists[0]
  if not shiftLists:
    raise ValueError("Function must have either a shiftList or peakLists with shiftLists")
  elif len(shiftLists) == 1:
    shiftList = shiftLists[0]

  # Peak list check
  for peakList in peakLists:
    if peakList._project is not project:
      raise ValueError("PeakList %s does not match Project %s"
                       % (peakList, project))
    if peakList.chemicalShiftList is None:
      # NBNB This is a side effect - but shift lists should always be set anyway
      peakList.chemicalShiftList = shiftList

  # Chain check
  ll = [x.id for x in chains]
  if len(ll) != len(set(ll)):
    raise ValueError("Chains parameter contains duplicates: %s" % chains)
  for chain in chains:
    if chain._project is not project:
      raise ValueError("Chain %s does not match Project %s"
                       % (chain, project))

  # restraint list check
  aSet = set(x.dataSet.id for x in restraintLists)
  if len(aSet) == 1:
    if restraintLists[0].dataSet._project is not project:
      raise ValueError("Restraint lists do not match Project")
  else:
    raise ValueError("Restraint lists are not from a single DataSet")

  # Make BMRB object tree :

  # Make Entry
  entry = bmrb.entry.fromScratch(dataName)

  # MetaData saveframe
  entry.addSaveframe(_createMetaDataFrame(programName=programName, programVersion=programVersion,
                      coordinateFileName=coordinateFileName, relatedEntries=relatedEntries,
                      programScripts=programScripts, programHistory=programHistory, uuid=uuid))

  # Make molecular system
  entry.addSaveframe(_createMolecularSystemFrame(chains))

  # Make shift lists
  for shiftList in shiftLists:
    entry.addSaveframe(createShiftListFrame(shiftList))

  # Make restraint lists
  for restraintList in restraintLists:
    entry.addSaveframe(createRestraintListFrame(restraintList))

  # Make peak list and spectrum frames
  for peakList in peakLists:
    entry.addSaveframe(createPeakListFrame(peakList))

  # Make Peak-restraint links frame
  entry.addSaveframe(createPeakRestraintLinksFrame(restraintLists, peakLists))

  return entry


def _createMetaDataFrame(programName:str='CCPN', programVersion:str=None,
                      coordinateFileName:str=None, relatedEntries:Sequence=(),
                      programScripts:Sequence=(), programHistory:Sequence=(),
                      uuid=None):
  """Make Metadata singleton saveframe.
  NBNB currently all values are dummy - later the function will need parameters"""

  if programVersion is None:
    programVersion = str(Version.currentModelVersion)

  # Make meta_data
  saveframe = bmrb.saveframe.fromScratch(saveframe_name='nef_nmr_meta_data',
                                         tag_prefix='nef_nmr_meta_data')

  saveframe.addTags([
    ('sf_category','nef_nmr_meta_data'),
    ('sf_framecode','nef_nmr_meta_data'),
    ('format_name','Nmr_Exchange_Format'),
    ('format_version',formatVersion),
    ('program_name',programName),
    ('program_version',programVersion),
  ])
  timeStamp = datetime.datetime.today().isoformat()
  saveframe.addTag('creation_date', timeStamp)
  if uuid is None:
    # NBNB TBD do we want to have microseconds in the timestamp?
    uuid = '%s-%s-%s' % (programName, timeStamp, random.randint(0,maxRandomInt))
  saveframe.addTag('uuid', uuid)
  saveframe.addTag('coordinate_file_name', coordinateFileName)

  # related entries loop
  if relatedEntries:
    loop = bmrb.loop.fromScratch(category='nef_related_entries')
    saveframe.addLoop(loop)
    for tag in ('database_name', 'database_accession_code'):
      loop.addColumn(tag)
    for dataTuple in relatedEntries:
      loop.addData(dataTuple)

  # Program script loop
  if programScripts:
    loop = bmrb.loop.fromScratch(category='nef_program_script')
    saveframe.addLoop(loop)
    for tag in ('program_name', 'script_name', 'script', 'ccpn_special_par'):
      loop.addColumn(tag)
    for dataTuple in programScripts:
      loop.addData(dataTuple)

  # Program history loop
  if programHistory:
    loop = bmrb.loop.fromScratch(category='nef_run_history')
    saveframe.addLoop(loop)
    for tag in ('run_ordinal', 'program_name', 'program_version', 'script_name', 'script'):
      loop.addColumn(tag)
    for dataTuple in programHistory:
      loop.addData(dataTuple)
  #
  return saveframe


def _createMolecularSystemFrame(chains):
  """ Make molecular system frame"""

  project = chains[0]._project

  # Header block
  category = 'nef_molecular_system'
  saveframe = bmrb.saveframe.fromScratch(saveframe_name=category,
                                         tag_prefix=category)
  saveframe.addTags([
    ('sf_category',category),
    ('sf_framecode',category),
  ])

  # sequence loop
  loop = bmrb.loop.fromScratch(category='nef_sequence')
  saveframe.addLoop(loop)
  for tag in ('chain_code', 'sequence_code', 'residue_type', 'linking', 'residue_variant'):
    loop.addColumn(tag)


  chainCodes = [x.shortName for x in chains]
  if len(chainCodes) != len(set(chainCodes)):
      raise ValueError("Duplicate chains not allowed: %s" % (chains,))
  for chain in chains:
    chainCode = chain.shortName
    residues = chain.residues
    if not residues:
      continue

    apiMolResidues = [x._wrappedData.molResidue for x in residues]
    isCyclic = (apiMolResidues[0].previousMolResidue is apiMolResidues[-1]
                and apiMolResidues[-1].nextMolResidue is apiMolResidues[0])

    for ii, residue in enumerate(residues):
      sequenceCode = residue.sequenceCode
      residueType = residue.residueType
      apiMolResidue = apiMolResidues[ii]
      isLinearPolymer = apiMolResidue.chemComp.isLinearPolymer

      # Set linking
      linking = useLinking = residue.linking
      if not isLinearPolymer:
        useLinking = 'nonlinear'

      elif linking == 'none':
        useLinking = 'single'

      elif ii == 0:
        if isCyclic:
          useLinking = 'cyclic'
        elif linking!= 'start':
          useLinking = 'break'

      elif residue is residues[-1]:
        if isCyclic:
          useLinking = 'cyclic'
        elif linking!= 'end':
          useLinking = 'break'

      else:
        if  apiMolResidue.nextMolResidue is not apiMolResidues[ii+1] and linking != 'end':
          useLinking = 'break'
        if apiMolResidue.previousMolResidue is not apiMolResidues[ii-1] and linking != 'start':
          useLinking = 'break'

      # NBNB TBD check residue variants - this is according to new proposal 1/3/2015
      chemCompVar = residue._wrappedData.chemCompVar
      chemComp = chemCompVar.chemComp
      if chemCompVar.isDefaultVar:
        residueVariant=None
      else:
        defaultVar = None
        if isLinearPolymer:
          defaultVar = chemComp.findFirstChemCompVar(isDefaultVar=True, linking=linking)
        if defaultVar is None:
          defaultVar = chemComp.findFirstChemCompVar(isDefaultVar=True, linking='none')
        if defaultVar is None:
          defaultVar = chemComp.findFirstChemCompVar(isDefaultVar=False, linking=linking)
        if defaultVar is None:
          defaultVar = chemComp.sortedChemCompVars()[0]

        atoms = chemCompVar.findAllChemAtoms(className='ChemAtom')
        defAtoms = defaultVar.findAllChemAtoms(className='ChemAtom')
        addNames = [x.name.translate(Pid.unmapSeparators) for x in (atoms - defAtoms)]
        removeNames = [x.name.translate(Pid.unmapSeparators) for x in (defAtoms - atoms)]
        residueVariant = ','.join(['+%s' %x for x in sorted(addNames)] +
                                  ['-%s' %x for x in sorted(removeNames)])
        if not residueVariant:
          residueVariant = None

      loop.addData([chainCode, sequenceCode, residueType, linking, residueVariant])

  # covalent cross-links loop

  atomPairs = []
  for chain in chains:
    molecule = chain._apiChain.molecule
    for molResLink in molecule.findAllMolResLinks(isStdLinear=False):
      pair = []
      atomPairs.append(pair)
      for molResLinkEnd in molResLink.molResLinkEnds:
        apiResidue = chain.findFirstResidue(seqId=molResLinkEnd.molResidue.serial)
        atomName = molResLinkEnd.linkEnd.boundChemAtom.name
        apiAtom = apiResidue.findFirstAtom(name=atomName)
        pair.append(project._data2Obj[apiAtom]._id.split('.'))
      pair.sort(key=CcpnSorting.universalSortKey)

  for molSystemLink in project._apiNmrProject.molSystem.molSystemLinks:
    pair = []
    atomPairs.append(pair)
    for molSystemLinkEnd in molSystemLink.molSystemLinkEnds:
      atomName = molSystemLinkEnd.linkEnd.boundChemAtom.name
      apiAtom = molSystemLinkEnd.residue.findFirstAtom(name=atomName)
      pair.append(project._data2Obj[apiAtom]._id.split('.'))
    pair.sort(key=CcpnSorting.universalSortKey)

  # NB Loop may be empty. The entry.export_string function takes care of this
  loop = bmrb.loop.fromScratch(category='nef_covalent_links')
  saveframe.addLoop(loop)
  for tag in ('chain_code_1', 'sequence_code_1', 'residue_type_1', 'atom_name_1',
              'chain_code_2', 'sequence_code_2', 'residue_type_2', 'atom_name_2',):
    loop.addColumn(tag)


  for ll in sorted(atomPairs, key=CcpnSorting.universalSortKey):
    values = ll[0]._id.split('.') +  ll[1]._id.split('.')
    loop.addData(values)
  #
  return saveframe


def createShiftListFrame(shiftList):
  """make a saveframe for a shift list"""

  sf_category = 'nef_chemical_shift_list'
  framecode = '%s_%s' % (sf_category, shiftList.name.translate(Pid.unmapSeparators))
  saveframe = bmrb.saveframe.fromScratch(saveframe_name=framecode,
                                         tag_prefix='nef_chemical_shift_list')

  saveframe.addTags([
    ('sf_category',sf_category),
    ('sf_framecode',framecode),
  ])

  saveframe.addTag('atom_chemical_shift_units', shiftList.unit)

  # chemical shift loop
  loop = bmrb.loop.fromScratch(category='nef_chemical_shift')
  saveframe.addLoop(loop)
  for tag in ('chain_code', 'sequence_code', 'residue_type', 'atom_name',
              'value', 'value_uncertainty',):
    loop.addColumn(tag)

  for x in shiftList.chemicalShifts:
    loop.addData((Pid.splitId(x.nmrAtom._id) + (x.value, x.valueError)))

  # sortkey = shiftList._project.universalSortKey
  # for row in sorted(((Pid.splitId(x.nmrAtom._id) + (x.value, x.valueError))
  #                   for x in shiftList.chemicalShifts), key=sortkey):
  #   loop.addData(row)
  #
  return saveframe

def createRestraintListFrame(restraintList):
  """make a saveframe for a restraint list of whatever type."""
  restraintType = restraintList.restraintType
  potentialType = restraintList.potentialType
  restraintItemLength = restraintList.itemLength

  # if restraintType == 'HBond':
  #   restraintListTag = 'distance'
  #   if restraintList.origin is None:
  #     restraintList.origin = 'hbond'
  # else:
  restraintListTag = restraintType.lower()
  sf_category = 'nef_%s_restraint_list' % restraintListTag
  framecode = '%s_%s' % (sf_category, restraintList.name.translate(Pid.unmapSeparators))

  saveframe = bmrb.saveframe.fromScratch(saveframe_name=framecode, tag_prefix=sf_category)

  saveframe.addTags([
    ('sf_category',sf_category),
    ('sf_framecode',framecode),
    ('potential_type',restraintList.potentialType),
    ('origin',restraintList.origin)
  ])

  # Rdc-specific tags:
  if restraintType == 'Rdc':
    saveframe.addTags([
      ('tensor_magnitude',restraintList.tensorMagnitude),
      ('tensor_rhombicity',restraintList.tensorRhombicity),
      ('tensor_chain_code',restraintList.tensorChainCode),
      ('tensor_sequence_code',restraintList.tensorSequenceCode),
      ('tensor_residue_type',restraintList.tensorResidueType),
    ])

  # restraint loop
  loop = bmrb.loop.fromScratch(category=sf_category[:-5])
  saveframe.addLoop(loop)
  # ID tags:
  for tag in ('ordinal', 'restraint_id', 'restraint_combination_id',):
    loop.addColumn(tag)
  # Assignment tags:
  for ii in range(restraintItemLength):
    for tag in ('chain_code_', 'sequence_code_', 'residue_type_', 'atom_name_'):
      ss = str(ii + 1)
      loop.addColumn(tag + ss)
  # fixed tags:
  for tag in ('weight', 'target_value', 'target_value_uncertainty',):
    loop.addColumn(tag)
  # Rdc-specific tags:
  if restraintType == 'Rdc':
    for tag in ('scale', 'distance_dependent'):
      loop.addColumn(tag)
  # potential tags
  tags = nefTagsByPotentialType.get(potentialType)
  if tags is None:
    tags = nefTagsByPotentialType.get('unknown')
  for tag in tags:
    loop.addColumn(tag)
  # Dihedral-specific tags:
  if restraintType == 'Dihedral':
    loop.addColumn('name')

  # Add data
  ordinal = 0
  for restraint in restraintList.restraints:
    serial = restraint.serial

    for contribution in restraint.restraintContributions:
      row1 = [serial, contribution.combinationId]
      row3 = [contribution.weight, contribution.targetValue, contribution.error]
      if restraintType == 'Rdc':
        row3.extend((contribution.scale, contribution.isDistanceDependent))
      for ss in nefTagsByPotentialType[restraintList.potentialType]:
        tag = nef2CcpnTags.get(ss,ss)
        row3.append(getattr(contribution, tag))

      for restraintItem in contribution.restraintItems:
        row2 = []
        # Assignment tags:
        for atomId in restraintItem:
          row2.extend(Pid.splitId(atomId))
        ordinal += 1
        ll = [ordinal] + row1 + row2 + row3
        if restraintType == 'Dihedral':
          ll.append(RestraintLib.dihedralName(project, restraintItem))
        loop.addData(ll)
  #
  return saveframe

def createPeakListFrame(peakList):
  """make saveframe for peakList an containing spectrum"""

  # Set up variables
  spectrum = peakList.spectrum
  dimensionCount = spectrum.dimensionCount
  apiDataDims = spectrum._apiDataSource.sortedDataDims()
  if peakList.title is None:
    peakList.title = '%s-%s' % (spectrum.name.translate(Pid.unmapSeparators), peakList.serial)
  category = 'nef_nmr_spectrum'
  framecode = '%s_%s' % (category, peakList.title)
  saveframe = bmrb.saveframe.fromScratch(saveframe_name=framecode,
                                         tag_prefix=category)
  # Top tagged values
  obj = peakList.chemicalShiftList
  if obj is None:
    shiftListString = None
  else:
    shiftListString = 'nef_chemical_shift_list_%s' % obj.name.translate(Pid.unmapSeparators)
  saveframe.addTags([
    ('sf_category',category),
    ('sf_framecode',framecode),
    ('num_dimensions',dimensionCount),
    ('chemical_shift_list', shiftListString)
  ])

  # Experiment type
  refExperiment = spectrum._apiDataSource.experiment.refExperiment
  if refExperiment is None:
    name = synonym = None
  else:
    name = refExperiment.name
    synonym = refExperiment.synonym
  saveframe.addTags([
    ('experiment_classification',name),
    ('expriment_type',synonym),
  ])


  # experiment dimensions loop
  loop = bmrb.loop.fromScratch(category='nef_spectrum_dimension')
  saveframe.addLoop(loop)
  for tag in ('dimension_id', 'axis_unit', 'axis_code', 'spectrometer_frequency',
              'spectral_width', 'value_first_point', 'folding', 'absolute_peak_positions',
              'is_acquisition'):
    loop.addColumn(tag)

  rows = zip(
    range(1, dimensionCount + 1),
    spectrum.axisUnits,
    spectrum.isotopeCodes,
    spectrum.spectrometerFrequencies,
    spectrum.spectralWidths,
    tuple(x.primaryDataDimRef.pointToValue(1. - x.pointOffset) for x in apiDataDims),
    # NBNB TBD this can break for non-Freq dimensions
    spectrum.foldingModes,
    (True,) * dimensionCount,
    # NBNB TBD we have no way of representing data for absPeakPositions=False
    tuple(x.expDim.isAcquisition for x in apiDataDims),
  )
  for row in rows:
    loop.addData(row)

  # dimension connection loop
  loop = bmrb.loop.fromScratch(category='nef_spectrum_dimension_transfer')
  saveframe.addLoop(loop)
  for tag in ('dimension_1', 'dimension_2', 'transfer_type', 'is_indirect'):
    loop.addColumn(tag)

  ref2Dim = {}
  expDimRefs = []
  for dim,dataDim in enumerate(apiDataDims):
    # Find ExpDimRef to use. NB, could break for unusual cases
    expDimRef = dataDim.expDim.sortedExpDimRefs()[0]
    expDimRefs.append(expDimRef)
    ref2Dim[expDimRef] = dim + 1

  for xd1 in expDimRefs:
    for xd2 in expDimRefs:
      dim1 = ref2Dim[xd1]
      dim2 = ref2Dim[xd2]
      if dim1 < dim2:
        tt = libSpectrum._expDimRefTransferType(xd1, xd2)
        if tt is not None:
          loop.addData([dim1, dim2, tt[0], not(tt[1])])


  # main peak loop
  loop = bmrb.loop.fromScratch(category='nef_peak')
  saveframe.addLoop(loop)
  for tag in ('ordinal', 'peak_id', 'volume', 'volume_uncertainty', 'height', 'height_uncertainty'):
    loop.addColumn(tag)
  for ii in range(1, dimensionCount + 1):
    loop.addColumn('position_%s' % ii)
    loop.addColumn('position_uncertainty_%s' % ii)
  for ii in range(1, dimensionCount + 1):
    loop.addColumn('chain_code_%s' % ii)
    loop.addColumn('sequence_code_%s' % ii)
    loop.addColumn('residue_type_%s' % ii)
    loop.addColumn('atom_name_%s' % ii)

  ordinal = 0
  for peak in peakList.peaks:

    # serial and intensities
    firstpart = [peak.serial, peak.volume, peak.volumeError, peak.height, peak.heightError]

    # peak position
    for ll in zip(peak.position, peak.positionError):
      firstpart.extend(ll)

    # assignments
    assignedNmrAtoms = peak.assignedNmrAtoms
    if assignedNmrAtoms:
      for assignment in assignedNmrAtoms:
        ordinal += 1
        row = [ordinal] + firstpart
        for nmrAtom in assignment:
          if nmrAtom is None:
            row.extend((None, None, None, None))
          else:
            row.extend(Pid.splitId(nmrAtom._id))
        loop.addData(row)
    else:
      # Unassigned peak
      ordinal += 1
      loop.addData([ordinal] + firstpart + [None] * (dimensionCount * 4))
  #
  return saveframe


def createPeakRestraintLinksFrame(restraintLists, peakLists):
  """ Make peak-constraint links frame, with links tha match both constraintLists and peakLists"""

  # header block
  category = 'nef_peak_restraint_links'
  saveframe = bmrb.saveframe.fromScratch(saveframe_name=category,
                                         tag_prefix=category)
  saveframe.addTags([
    ('sf_category',category),
    ('sf_framecode',category),
  ])

  # peak restraint links loop
  loop = bmrb.loop.fromScratch(category='nef_peak_restraint_link')
  saveframe.addLoop(loop)
  for tag in ('nmr_spectrum_id', 'peak_id', 'restraint_list_id', 'restraint_id'):
    loop.addColumn(tag)

  data = []
  for restraintList in restraintLists:
    restraintType = restraintList.restraintType
    if restraintType == 'HBond':
      restraintListTag = 'distance'
    else:
      restraintListTag = restraintType.lower()
    restraint_list_id = ( 'nef_%s_restraint_list_%s' %
                          (restraintListTag, restraintList.name.translate(Pid.unmapSeparators)))

    for restraint in restraintList.restraints:
      restraint_id = restraint.serial
      for peak in restraint.peaks:
        if peak.peakList in peakLists:
          data.append(['nef_nmr_spectrum_' + peak.peakList.title, peak.serial,
                       restraint_list_id, restraint_id])
  for ll in sorted(data):
      loop.addData(ll)
  #
  return saveframe


# # No longer necessary - it is all dealt with elsewhere
# def prepareNmrProject(nmrProject):
#   """Fix non-CCPN-generated NmrProjects that violate normal assumptions"""
#
#   # If there is only one shiftList, make sure all experiments have it set
#   shiftLists = nmrProject.findAllMeasurementLists(className='ShiftList')
#   if len(shiftLists) == 1:
#     shiftList = shiftLists.pop()
#     for experiment in nmrProject.experiments:
#       experiment.shiftList = shiftList
#

  #
  # # fix restraint list names to be unique
  # constraintLists = [y for x in nmrProject.sortedNmrConstraintStores()
  #                      for y in x.sortedConstraintLists()]
  #
  # names = [x.name for x in constraintLists]
  # for ii,name in enumerate(names):
  #   obj = constraintLists[ii]
  #   if name:
  #     names[ii] =  '_'.join(obj.name.split())
  #   else:
  #     names[ii] = ("RestraintList:%s.%s,%s" %
  #                  (obj.nmrConstraintStore.serial, obj.className[:-14], obj.serial))
  #
  # for ii,name in enumerate(names):
  #   obj = constraintLists[ii]
  #   if names.count(name) > 1:
  #     name = ("RestraintList:%s.%s,%s" %
  #             (obj.nmrConstraintStore.serial, obj.className[:-14], obj.serial))
  #     names[ii] = name
  #   constraintLists[ii].name = name


if __name__ == '__main__':

  if len(sys.argv) >= 3:

    from ccpn.core._implementation import Io as ccpnIo

    # from ccpnmodel.ccpncore.lib.Io.Api import loadProject
    # from ccpn.core.Project import Project
    #
    # # set up input
    # junk, projectDir, outputDir = sys.argv[:3]
    # apiProject = loadProject(projectDir)
    #
    # # # NBNB TBD this must be integrated in upgrade once we are in a stable version
    # # from ccpnmodel.v_3_0_2.upgrade import correctFinalResult
    # # correctFinalResult(apiProject)
    #
    # nmrProject = apiProject.findFirstNmrProject()
    # prepareNmrProject(nmrProject)
    # pp = Project(nmrProject)

    # set up input
    junk, projectDir, outputDir = sys.argv[:3]
    project = ccpnIo.loadProject(projectDir)
    apiProject = project._wrappedData.root

    tt = Path.splitPath(projectDir)
    if not tt[1]:
      tt = Path.splitPath(tt[0])

    if len(sys.argv) >= 4:
      # set up input
      constraintStoreSerial = int(sys.argv[4])
      constraintStore = apiProject.findFirstNmrConstraintStore(serial=constraintStoreSerial)
    else:
      constraintStore = apiProject.findFirstNmrConstraintStore()

    exportRestraintStore(project._data2Obj[constraintStore], dataName=tt[1],
                         directory=outputDir, forceFirstShiftList=True)

  else:
    print ("Error. Parameters are: ccpnProjectDirectory outputDirectory [constraintStoreSerial] ")