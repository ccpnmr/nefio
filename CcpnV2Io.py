"""NEF I/O for CCPN V2 release, data model version 2.1.2

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================

__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan"
               "Simon P Skinner & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = (
"For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license"
"or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = ": Rasmus H Fogh $"
__dateModified__ = ": 2017-04-07 11:40:45 +0100 (Fri, April 07, 2017) $"
__version__ = ": 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = ": Rasmus H Fogh $"

__date__ = ": 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

# NB must be Python 2.7 and 3.x compatible


# Max value used for random integer. Set to be expressible as a signed 32-bit integer.
maxRandomInt =  2000000000

#  - saveframe category names in reading order
# The order is significant, because setting of crosslinks relies on the order frames are read
# Frames are read in correct order regardless of how they are in the file
saveFrameReadingOrder = [
  'nef_nmr_meta_data',
  'nef_molecular_system',
  'ccpn_sample',
  'ccpn_substance',
  'ccpn_assignment',
  'nef_chemical_shift_list',
  'ccpn_dataset',
  'nef_distance_restraint_list',
  'nef_dihedral_restraint_list',
  'nef_rdc_restraint_list',
  'nef_nmr_spectrum',
  'nef_peak_restraint_links',
  'ccpn_complex',
  'ccpn_spectrum_group',
  'ccpn_restraint_list',
  'ccpn_notes',
  'ccpn_additional_data'
]


def loadNefFile(memopsRoot, path):
  """Load NEF file at path into memopsRoot"""

  nefReader = CcpnNefReader()
  dataBlock = nefReader.getNefData(path)
  nefReader.importNewProject(memopsRoot, dataBlock)


class CcpnNefReader:
  # Importer functions - used for converting saveframes and loops
  importers = {}

  def __init__(self, application:str, specificationFile:str=None, mode:str='standard',
               testing:bool=False):

    self.application = application
    self.mode=mode
    self.saveFrameName = None
    self.warnings = []
    self.errors = []
    self.testing = testing

    # Map for resolving crosslinks in NEF file
    self.frameCode2Object = {}

    # Map for speeding up restraint reading
    self._dataSet2ItemMap = None
    self._nmrResidueMap = None

    self.defaultDataSetSerial = None
    self.defaultNmrChain = None
    self.mainDataSetSerial = None
    self.defaultChemicalShiftList = None


  def getNefData(self, path:str):
    """Get NEF data structure from file"""
    nmrDataExtent = StarIo.parseNefFile(path)
    dataBlocks = list(nmrDataExtent.values())
    dataBlock = dataBlocks[0]

    # Initialise afresh for every file read
    self._dataSet2ItemMap = {}
    self._nmrResidueMap = {}
    #
    return dataBlock

  def _getSaveFramesInOrder(self, dataBlock:StarIo.NmrDataBlock) -> OD:
    """Get saveframes in fixed reading order as Ordereddict(category:[saveframe,])"""
    result = OD(((x, []) for x in saveFrameReadingOrder))
    result['other'] = otherFrames = []
    for saveFrameName, saveFrame in dataBlock.items():
      sf_category = saveFrame.get('sf_category')
      ll = result.get(sf_category)
      if ll is None:
        ll = otherFrames
      ll.append(saveFrame)
    #
    return result

  def importNewProject(self, project:Project, dataBlock:StarIo.NmrDataBlock,
                       projectIsEmpty:bool=True):
    """Import entire project from dataBlock into empty Project"""

    t0 = time.time()


    # TODO Add error handling

    self.warnings = []
    self.project = project
    self.defaultChainCode = None

    saveframeOrderedDict = self._getSaveFramesInOrder(dataBlock)

    # Load metadata and molecular system first
    metaDataFrame = dataBlock['nef_nmr_meta_data']
    self.saveFrameName = 'nef_nmr_meta_data'
    self.load_nef_nmr_meta_data(project, metaDataFrame)
    del saveframeOrderedDict['nef_nmr_meta_data']

    saveFrame = dataBlock.get('nef_molecular_system')
    if saveFrame:
      self.saveFrameName = 'nef_molecular_system'
      self.load_nef_molecular_system(project, saveFrame)
    del saveframeOrderedDict['nef_molecular_system']

    # Load assignments, or preload from shiftlists
    # to make sure '@' and '#' identifiers match the right serials
    saveFrame = dataBlock.get('ccpn_assignment')
    if saveFrame:
      self.saveFrameName = 'ccpn_assignment'
      self.load_ccpn_assignment(project, saveFrame)
      del saveframeOrderedDict['ccpn_assignment']
    else:
      self.preloadAssignmentData(dataBlock)

    # t1 = time.time()
    # print ('@~@~ NEF load starting frames', t1-t0)

    for sf_category, saveFrames in saveframeOrderedDict.items():
      for saveFrame in saveFrames:
        saveFrameName = self.saveFrameName = saveFrame.name

        importer = self.importers.get(sf_category)
        if importer is None:
          print ("WARNING, unknown saveframe category", sf_category, saveFrameName)
        else:
          # NB - newObject may be project, for some saveframes.
          result = importer(self, project, saveFrame)
          if isinstance(result, AbstractWrapperObject):
            self.frameCode2Object[saveFrameName] = result
          # elif not isinstance(result, list):
          #   self.warning("Unexpected return %s while reading %s" %
          #                (result, saveFrameName))

          # Handle unmapped elements
          extraTags = [x for x in saveFrame
                       if x not in nef2CcpnMap[sf_category]
                       and x not in ('sf_category', 'sf_framecode')]
          if extraTags:
            print("WARNING - unused tags in saveframe %s: %s" % (saveFrameName, extraTags))
            # TODO put here function that stashes data in object, or something
            # ues newObject here
          # t2 = time.time()
          # print('@~@~ loaded', saveFrameName, t2-t1)
          # t1 = t2

    # Put metadata in main dataset
    self.updateMetaData(metaDataFrame)


    t2 = time.time()
    print('Loaded NEF file, time = ', t2-t0)

    for msg in self.warnings:
      print ('====> ', msg)
    self.project = None


  def load_nef_nmr_meta_data(self, project:Project, saveFrame:StarIo.NmrSaveFrame):
    """load nef_nmr_meta_data saveFrame"""

    # Other data are read in here at the end of the load
    self.mainDataSetSerial = saveFrame.get('ccpn_dataset_serial')

    return None

    # TODO - store data in this saveframe
    # for now we store none of this, as the storage slots are in DataSet, not Project
    # Maybe for another load function?
  #
  importers['nef_nmr_meta_data'] = load_nef_nmr_meta_data


  def load_nef_molecular_system(self, project:Project, saveFrame:StarIo.NmrSaveFrame):
    """load nef_molecular_system saveFrame"""

    mapping = nef2CcpnMap['nef_molecular_system']
    for tag, ccpnTag in mapping.items():
      if ccpnTag == _isALoop:
        loop = saveFrame.get(tag)
        if loop:
          importer = self.importers[tag]
          importer(self, project, loop)
    #
    return None
  #
  importers['nef_molecular_system'] = load_nef_molecular_system


  def load_nef_sequence(self, project:Project, loop:StarIo.NmrLoop):
    """Load nef_sequence loop"""

    result = []

    chainData = {}
    for row in loop.data:
      chainCode = row['chain_code']
      ll = chainData.get(chainCode)
      if ll is None:
        chainData[chainCode] = [row]
      else:
        ll.append(row)

    defaultChainCode = None
    if None in chainData:
      defaultChainCode = 'A'
      # Replace chainCode None with default chainCode
      # Selecting the first value that is not already taken.
      while defaultChainCode in chainData:
        defaultChainCode = commonUtil.incrementName(defaultChainCode)
      chainData[defaultChainCode] = chainData.pop(None)
    self.defaultChainCode = defaultChainCode


    sequence2Chain = {}
    tags =('residue_name', 'linking', 'residue_variant')
    for chainCode, rows in sorted(chainData.items()):
      compoundName = rows[0].get('ccpn_compound_name')
      role = rows[0].get('ccpn_chain_role')
      comment = rows[0].get('ccpn_chain_comment')
      for row in rows:
        if row.get('linking') == 'dummy':
          row['residue_name'] = 'dummy.' + row['residue_name']
      sequence = tuple(tuple(row.get(tag) for tag in tags) for row in rows)

      lastChain = sequence2Chain.get(sequence)
      if lastChain is None:
        newSubstance = project.fetchNefSubstance(sequence=rows, name=compoundName)
        newChain = newSubstance.createChain(shortName=chainCode, role=role,
                                            comment=comment)
        sequence2Chain[sequence] = newChain

        # Set variant codes:
        for ii, residue in enumerate(newChain.residues):
          variantCode = sequence[ii][2]

          if variantCode:

            atomNamesRemoved, atomNamesAdded = residue._wrappedData.getAtomNameDifferences()


            for code in variantCode.split(','):
              code = code.strip()  # Should not be necessary but costs nothing to catch those errors
              atom = residue.getAtom(code[1:])
              if code[0] == '-':
                if atom is None:
                  residue._project._logger.error(
                    "Incorrect variantCode %s: No atom named %s found in %s. Skipping ..."
                    % (variantCode, code, residue)
                  )
                else:
                  atom.delete()

              elif code[0] == '+':
                if atom is None:
                  residue.newAtom(name=code[1:])
                else:
                  residue._project._logger.error(
                    "Incorrect variantCode %s: Atom named %s already present in %s. Skipping ..."
                    % (variantCode, code, residue)
                  )

              else:
                residue._project._logger.error(
                  "Incorrect variantCode %s: must start with '+' or '-'. Skipping ..."
                  % variantCode
                )

      else:
        newChain = lastChain.clone(shortName=chainCode)
        newChain.role = role
        newChain.comment = comment

      for apiResidue in newChain._wrappedData.sortedResidues():
        # Necessary to guarantee against name clashes
        # Direct access to avoid unnecessary notifiers
        apiResidue.__dict__['seqInsertCode'] = '__@~@~__'
      for ii,apiResidue in enumerate(newChain._wrappedData.sortedResidues()):
        # NB we have to loop over API residues to be sure we get the residues
        # in creation order rather than sorted order
        residue = project._data2Obj[apiResidue]
        residue.rename(rows[ii].get('sequence_code'))
        residue._resetIds()

      # Necessary as notification is blanked here:
      newChain._resetIds()

      #
      result.append(newChain)

      # Add Residue comments
      for ii, apiResidue in enumerate(newChain.residues):
        comment = rows[ii].get('ccpn_comment')
        if comment:
          apiResidue.details = comment
    #
    return result
  #
  importers['nef_sequence'] = load_nef_sequence


  def load_nef_covalent_links(self, project:Project, loop:StarIo.NmrLoop):
    """Load nef_sequence loop"""

    result = []

    for row in loop.data:
      id1 = Pid.createId(*(row[x] for x in ('chain_code_1', 'sequence_code_1',
                                           'residue_name_1', 'atom_name_1', )))
      id2 = Pid.createId(*(row[x] for x in ('chain_code_2', 'sequence_code_2',
                                           'residue_name_2', 'atom_name_2', )))
      atom1 = project.getAtom(id1)
      atom2 = project.getAtom(id2)
      if atom1 is None:
        self.warning("Unknown atom %s for bond to %s. Skipping..." % (id1, id2))
      elif atom2 is None:
        self.warning("Unknown atom %s for bond to %s. Skipping..." % (id2, id1))
      else:
        result.append((atom1, atom2))
        atom1.addInterAtomBond(atom2)
    #
    return result
  #
  importers['nef_covalent_links'] = load_nef_covalent_links



  def preloadAssignmentData(self, dataBlock:StarIo.NmrDataBlock):
    """Set up NmrChains and NmrResidues with reserved names to ensure the serials are OK
    and create NmrResidues in connencted nmrChains in order

    NB later we can store serials in CCPN projects, but something is needed that works anyway

    NB, without CCPN-specific tags you can NOT guarantee that connected stretches are stable,
    and that serials are put back where they came from.
    This heuristic creates NmrResidues in connected stretches in the order they are found,
    but this will break if connected tretches appear in multiple shiftlists and some are partial."""


    project = self.project

    for saveFrameName, saveFrame in dataBlock.items():

      # get all NmrResidue data in chemicalshift lists
      assignmentData = {}
      assignmentData2 = {}
      if saveFrameName.startswith('nef_chemical_shift_list'):
        loop = saveFrame.get('nef_chemical_shift')
        if loop:
          for row in loop.data:
            chainCode = row['chain_code']
            nmrResidues = assignmentData.get(chainCode, OD())
            assignmentData[chainCode] = nmrResidues
            nmrResidues[(row['sequence_code'], row['residue_name'])] = None

      # Create objects with reserved names
      for chainCode in sorted(assignmentData):
        if chainCode[0] in '@#' and chainCode[1:].isdigit():
          # reserved name - make chain
          try:
            project.fetchNmrChain(chainCode)
          except ValueError:
            # Could not be done, probably because we have NmrChain '@1'. Leave for later
            pass

      for chainCode, nmrResidues in sorted(assignmentData.items()):

        # Create NmrChain
        try:
          nmrChain = project.fetchNmrChain(chainCode)
        except ValueError:
          nmrChain = project.fetchNmrChain('`%s`' % chainCode)

        if nmrChain.isConnected:
          # Save data for later processing
          assignmentData2[nmrChain] = nmrResidues
        else:
          # Create non-assigned NmrResidues to reserve the serials. The rest can wait
          for sequenceCode, residueType in list(nmrResidues.keys()):
            if sequenceCode[0] == '@' and sequenceCode[1:].isdigit():
              nmrChain.fetchNmrResidue(sequenceCode=sequenceCode, residueType=residueType)

      for nmrChain, nmrResidues in sorted(assignmentData2.items()):
        # Create NmrResidues in order, to preserve connection order
        for sequenceCode, residueType in list(nmrResidues.keys()):
          # This time we want all non-offset, regardless of type - as we must get them in order
          if (len(sequenceCode) < 2 or sequenceCode[-2] not in '+-'
              or not sequenceCode[-1].isdigit()):
            # I.e. for sequenceCodes that do not include an offset
            nmrChain.fetchNmrResidue(sequenceCode=sequenceCode, residueType=residueType)


  def load_nef_chemical_shift_list(self, project:Project, saveFrame:StarIo.NmrSaveFrame):
    """load nef_chemical_shift_list saveFrame"""

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    parameters['name'] = framecode[len(category) + 1:]

    # Make main object
    result = project.newChemicalShiftList(**parameters)

    if self.defaultChemicalShiftList is None:
      # ChemicalShiftList should default to the unique ChemicalShIftList in the file
      # A file with multiple ChemicalShiftLists MUST have explicit chemical shift lists
      # given for all spectra- but this is nto hte place for validity checking
      self.defaultChemicalShiftList = result

    if self.testing:
      # When testing you want the values to remain as read
      result.autoUpdate = False
      # NB The above is how it ought to work.
      # The below is how it is working as of July 2016
      result._wrappedData.topObject.shiftAveraging = False

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        importer(self, result, loop)
    #
    return result
  #

  importers['nef_chemical_shift_list'] = load_nef_chemical_shift_list

  def load_nef_chemical_shift(self, parent:ChemicalShiftList, loop:StarIo.NmrLoop):
    """load nef_chemical_shift loop"""

    result = []

    creatorFunc = parent.newChemicalShift

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    for row in loop.data:
      parameters = self._parametersFromLoopRow(row, map2)
      tt = tuple(row.get(tag) for tag in ('chain_code', 'sequence_code', 'residue_name',
                                          'atom_name'))
      element = row.get('element')
      isotope = row.get('isotope_number')
      if element:
        if isotope:
          isotopeCode = '%s%s' % (isotope, element.title())
        else:
          isotopeCode = Constants.DEFAULT_ISOTOPE_DICT.get(element.upper())
      elif isotope:
        element = commonUtil.name2ElementSymbol(tt[3])
        isotopeCode = '%s%s' % (isotope, element.title())
      else:
        isotopeCode = None
      try:
        nmrResidue = self.produceNmrResidue(*tt[:3])
        nmrAtom = self.produceNmrAtom(nmrResidue, tt[3], isotopeCode=isotopeCode)
        parameters['nmrAtom'] = nmrAtom
        result.append(creatorFunc(**parameters))

      except ValueError:
        self.warning("Cannot produce NmrAtom for assignment %s. Skipping ChemicalShift" % (tt,))
        # Should eventually be removed - raise while still testing
        raise
    #
    return result
  #
  importers['nef_chemical_shift'] = load_nef_chemical_shift

  def load_nef_restraint_list(self, project:Project, saveFrame:StarIo.NmrSaveFrame):
    """Serves to load nef_distance_restraint_list, nef_dihedral_restraint_list,
     nef_rdc_restraint_list and ccpn_restraint_list"""

    # Get ccpn-to-nef mapping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    if category == 'nef_distance_restraint_list':
      restraintType = 'Distance'
    elif category == 'nef_dihedral_restraint_list':
      restraintType = 'Dihedral'
    elif category == 'nef_rdc_restraint_list':
      restraintType = 'Rdc'
    else:
      restraintType = saveFrame.get('restraint_type')
      if not restraintType:
        self.warning("Missing restraint_type for saveFrame %s - value was %s" %
                     (framecode, restraintType))
        return
    parameters['restraintType'] = restraintType
    namePrefix = restraintType[:3].capitalize() + '-'

    # Get name from frameCode, add type disambiguation, and correct for ccpn dataSetSerial addition
    name = framecode[len(category) + 1:]
    dataSetSerial = saveFrame.get('ccpn_dataset_serial')
    if dataSetSerial is not None:
      ss = '`%s`' % dataSetSerial
      if name.startswith(ss):
        name = name[len(ss):]

    # Make main object
    dataSet = self.fetchDataSet(dataSetSerial)
    previous = dataSet.getRestraintList(name)
    if previous is not None:
      # NEF but NOT CCPN has separate namespaces for different restraint types
      # so we may get name clashes
      # We should preserve NEF names, but it cannot be helped.
      if not name.startswith(namePrefix):
        # Add prefix for disambiguation since NEF but NOT CCPN has separate namespaces
        # for different constraint types
        name = namePrefix + name
        while dataSet.getRestraintList(name) is not None:
          # This way we get a unique name even in the most bizarre cases
          name = '`%s`' % name

    parameters['name'] = name

    result = dataSet.newRestraintList(**parameters)

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        if loopName.endswith('_restraint'):
          # NBNB HACK: the restrain loop reader needs an itemLength.
          # There are no other loops currently, but if there ever is they will not need this
          itemLength = saveFrame.get('restraint_item_length')
          importer(self, result, loop, itemLength)
        else:
          importer(self, result, loop)
    #
    return result
  #
  importers['nef_distance_restraint_list'] = load_nef_restraint_list
  importers['nef_dihedral_restraint_list'] = load_nef_restraint_list
  importers['nef_rdc_restraint_list'] = load_nef_restraint_list
  importers['ccpn_restraint_list'] = load_nef_restraint_list

  def load_nef_restraint(self, restraintList:RestraintList, loop:StarIo.NmrLoop,
                         itemLength:int=None):
    """Serves to load nef_distance_restraint, nef_dihedral_restraint,
     nef_rdc_restraint and ccpn_restraint loops"""

    # NB Restraint.name - written out for dihedral restraints - is not read.
    # Which is probably OK, it is derived from the atoms.

    result = []

    string2ItemMap = self._dataSet2ItemMap[restraintList.dataSet]

    # set itemLength if not passed in:
    if not itemLength:
      itemLength = coreConstants.constraintListType2ItemLength.get(restraintList.restraintType)

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    contributionTags = sorted(map2.values())
    restraints = {}
    # assignTags = ('chain_code', 'sequence_code', 'residue_name', 'atom_name')

    max = itemLength + 1
    multipleAttributes = OD((
      ('chainCodes',tuple('chain_code_%s' % ii for ii in range(1, max))),
      ('sequenceCodes',tuple('sequence_code_%s' % ii for ii in range(1, max))),
      ('residueTypes',tuple('residue_name_%s' % ii for ii in range(1, max))),
      ('atomNames',tuple('atom_name_%s' % ii for ii in range(1, max))),
    ))

    parametersFromLoopRow = self._parametersFromLoopRow
    defaultChainCode = self.defaultChainCode
    for row in loop.data:

      # get or make restraint
      serial = row.get('restraint_id')
      restraint = restraints.get(serial)
      if restraint is None:
        valuesToContribution = {}
        dd = {'serial':serial}
        val = row.get('ccpn_vector_length')
        if val is not None:
          dd['vectorLength'] = val
        val = row.get('ccpn_figure_of_Merit')
        if val is not None:
          dd['figureOfMerit'] = val
        val = row.get('ccpn_comment')
        if val is not None:
          dd['comment'] = val
        restraint = restraintList.newRestraint(**dd)
        restraints[serial] = restraint
        result.append(restraint)

      # Get or make restraintContribution
      parameters = parametersFromLoopRow(row, map2)
      combinationId = parameters.get('combinationId')
      nonAssignmentValues = tuple(parameters.get(tag) for tag in contributionTags)
      if combinationId:
        # Items in a combination are ANDed, so each line has one contribution
        contribution = restraint.newRestraintContribution(**parameters)
      else:
        contribution = valuesToContribution.get(nonAssignmentValues)
        if contribution is None:
          contribution = restraint.newRestraintContribution(**parameters)
          valuesToContribution[nonAssignmentValues] = contribution

      # Add item
      # ll = [row._get(tag)[:itemLength] for tag in assignTags]
      ll = [list(row.get(x) for x in y) for y in multipleAttributes.values()]
      # Reset missing chain codes to default
      # ll[0] = [x or defaultChainCode for x in ll[0]]

      idStrings = []
      for item in zip(*ll):
        if defaultChainCode is not None and item[0] is None:
          # ChainCode missing - replace with default chain code
          item = (defaultChainCode,) + item[1:]
        idStrings.append(Pid.IDSEP.join(('' if x is None else str(x)) for x in item))
      try:
        contribution.addRestraintItem(idStrings, string2ItemMap)
      except ValueError:
        self.warning("Cannot Add restraintItem %s. Identical to previous. Skipping" % idStrings)

    #
    return result
  #
  importers['nef_distance_restraint'] = load_nef_restraint
  importers['nef_dihedral_restraint'] = load_nef_restraint
  importers['nef_rdc_restraint'] = load_nef_restraint
  importers['ccpn_restraint'] = load_nef_restraint

  def load_nef_nmr_spectrum(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    dimensionTransferTags = ('dimension_1', 'dimension_2', 'transfer_type', 'is_indirect')

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    # Get peakList parameters and make peakList
    peakListParameters, dummy = self._parametersFromSaveFrame(saveFrame, mapping)

    # Get spectrum parameters
    spectrumParameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping,
                                                                  ccpnPrefix='spectrum')

    # Get name from spectrum parameters, or from the frameCode
    spectrumName = framecode[len(category) + 1:]
    peakListSerial = peakListParameters.get('serial')
    if peakListSerial:
      ss = '`%s`' % peakListSerial
      # Remove peakList serial suffix (which was added for disambiguation)
      # So that multiple peakLists all go to one Spectrum
      if spectrumName.endswith(ss):
        spectrumName = spectrumName[:-len(ss)]

    spectrum = project.getSpectrum(spectrumName)
    if spectrum is None:
      # Spectrum does not already exist - create it.
      # NB For CCPN-exported projects spectra with multiple peakLists are handled this way

      frameCode = saveFrame.get('chemical_shift_list')
      if frameCode:
        spectrumParameters['chemicalShiftList'] = self.frameCode2Object[frameCode]
      else:
        # Defaults to first (there should be only one, but we want the read to work) ShiftList
        spectrumParameters['chemicalShiftList'] = self.defaultChemicalShiftList


      frameCode = saveFrame.get('ccpn_sample')
      if frameCode:
        spectrumParameters['sample'] = self.frameCode2Object[frameCode]

      # get per-dimension data - NB these are mandatory and cannot be worked around
      dimensionData = self.read_nef_spectrum_dimension(project,
                                                       saveFrame['nef_spectrum_dimension'])
      # read  dimension transfer data
      loopName = 'nef_spectrum_dimension_transfer'
      # Those are treated elsewhere
      loop = saveFrame.get(loopName)
      if loop:
        data = loop.data
        transferData = [
          SpectrumLib.MagnetisationTransferTuple(*(row.get(tag) for tag in dimensionTransferTags))
          for row in data
        ]
      else:
        transferData = []

      spectrum = createSpectrum(project, spectrumName, spectrumParameters, dimensionData,
                                transferData=transferData)

      # Set experiment transfers at the API level
      if transferData and not spectrum.magnetisationTransfers:
        spectrum._setMagnetisationTransfers(transferData)

      # Make data storage object
      filePath = saveFrame.get('ccpn_spectrum_file_path')
      if filePath:
        storageParameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping,
          ccpnPrefix='spectrum._wrappedData.dataStore'
        )
        storageParameters['numPoints'] = spectrum.pointCounts
        spectrum._wrappedData.addDataStore(filePath, **storageParameters)

      # Load CCPN dimensions before peaks
      loopName = 'ccpn_spectrum_dimension'
      # Those are treated elsewhere
      loop = saveFrame.get(loopName)
      if loop:
        self.load_ccpn_spectrum_dimension(spectrum, loop)

    # Make PeakLst
    peakList = spectrum.newPeakList(**peakListParameters)

    # Load peaks
    self.load_nef_peak(peakList, saveFrame.get('nef_peak'))

    # Load remaining loops, with spectrum as parent
    for loopName in loopNames:
      if loopName not in  ('nef_spectrum_dimension', 'ccpn_spectrum_dimension', 'nef_peak',
                           'nef_spectrum_dimension_transfer'):
        # Those are treated elsewhere
        loop = saveFrame.get(loopName)
        if loop:
          importer = self.importers[loopName]
          importer(self, spectrum, loop)
    #
    return peakList
  #
  importers['nef_nmr_spectrum'] = load_nef_nmr_spectrum


  def read_nef_spectrum_dimension_transfer(self, loop:StarIo.NmrLoop):

    transferTypes = ('onebond', 'Jcoupling', 'Jmultibond', 'relayed', 'through-space',
                     'relayed-alternate')

    result = []

    if loop:
      data = loop.data
      for row in data:
        ll = [row.get(tag) for tag in ('dimension_1', 'dimension_2', 'transfer_type',
                                       'is_indirect')]
        result.append(SpectrumLib.MagnetisationTransferTuple(*ll))
    #
    return result


  def load_nef_spectrum_dimension_transfer(self, spectrum:Spectrum, loop:StarIo.NmrLoop):

    transferTypes = ('onebond', 'Jcoupling', 'Jmultibond', 'relayed', 'through-space',
                     'relayed-alternate')

    result = []

    if loop:
      apiExperiment = spectrum._wrappedData.experiment

      data = loop.data
      # Remove invalid data rows
      for row in data:
        ll = [row.get(tag) for tag in ('dimension_1', 'dimension_2', 'transfer_type',
                                       'is_indirect')]
        if (apiExperiment.findFirstExpDim(dim=row['dimension_1']) is None or
            apiExperiment.findFirstExpDim(dim=row['dimension_2']) is None or
            row['transfer_type'] not in transferTypes):
          self.warning("Illegal values in nef_spectrum_dimension_transfer: %s"
                       % list(row.values()))
        else:
          result.append(SpectrumLib.MagnetisationTransferTuple(*ll))
    #
    return result

  def process_nef_spectrum_dimension_transfer(self, spectrum:Spectrum,
                                              dataLists:Sequence[Sequence]):
      # Store expTransfers in API as we can not be sure we will get a refExperiment

      apiExperiment = spectrum._wrappedData.experiment

      for ll in dataLists:
        expDimRefs = []
        for dim in ll[:2]:
          expDim = apiExperiment.findFirstExpDim(dim=dim)
          # After spectrum creation there will be one :
          expDimRefs.append(expDim.sortedExpDimRefs()[0])
        if apiExperiment.findFirstExpTransfer(expDimRefs=expDimRefs) is None:
          apiExperiment.newExpTransfer(expDimRefs=expDimRefs, transferType=ll[2],
                                       isDirect=not ll[3])
        else:
          self.warning("Duplicate nef_spectrum_dimension_transfer: %s" % (ll,))


  def load_ccpn_spectrum_dimension(self, spectrum:Spectrum, loop:StarIo.NmrLoop) -> dict:
    """Read ccpn_spectrum_dimension loop, set the relevant values,
    and return the spectrum and other parameters for further processing"""

    params = {}
    extras = {}
    nefTag2ccpnTag = nef2CcpnMap['ccpn_spectrum_dimension']

    if not loop.data:
      raise ValueError("ccpn_spectrum_dimension is empty")

    rows = sorted(loop.data, key=itemgetter('dimension_id'))

    # Get spectrum attributes
    for nefTag, ccpnTag in nefTag2ccpnTag.items():
      if nefTag in rows[0]:
        ll = list(row.get(nefTag) for row in rows)
        if any(x is not None for x in ll):
          if ccpnTag and not '.' in ccpnTag:
            params[ccpnTag] = ll
          else:
            extras[nefTag] = ll

    # Set main values
    for tag, value in params.items():
      if tag != 'referencePoints':
        setattr(spectrum, tag, value)

    referencePoints = params.get('referencePoints')
    points = []
    values = []
    if referencePoints is not None:
      spectrumReferences = spectrum.spectrumReferences
      for ii, spectrumReference in enumerate(spectrumReferences):
        if spectrumReference is None:
          points.append(None)
          values.append(None)
        else:
          point = referencePoints[ii]
          points.append(point)
          values.append(spectrumReference.pointToValue(point))
    spectrum.referencePoints = points
    spectrum.referenceValues = values




    # set storage attributes
    value = extras.get('dimension_is_complex')
    if value:
      spectrum._wrappedData.dataStore.isComplex = value
    value = extras.get('dimension_block_size')
    if value:
      spectrum._wrappedData.dataStore.blockSizes = value

    # set aliasingLimits
    defaultLimits = spectrum.dimensionCount * [None]
    lowerLimits = extras.get('lower_aliasing_limit', defaultLimits)
    higherLimits = extras.get('higher_aliasing_limit', defaultLimits)
    spectrum.aliasingLimits = list(zip(lowerLimits, higherLimits))
  #
  importers['ccpn_spectrum_dimension'] = load_ccpn_spectrum_dimension

  # def adjustAxisCodes(self, spectrum, dimensionData):
  #
  #   pass
  #   #print ('@~@~ CCPN data. Still TODO')
  #
  #   # Use data to rename axisCodes
  #   axisCodes = spectrum.axisCodes
  #   newCodes = list(axisCodes)
  #   atomTypes = [commonUtil.splitIntFromChars(x)[1] for x in spectrum.isotopeCodes]
  #   acquisitionAxisCode = spectrum.acquisitionAxisCode
  #   if acquisitionAxisCode is not None:
  #     acquisitionDim = axisCodes.index(acquisitionAxisCode) + 1
  #     if acquisitionAxisCode == atomTypes[acquisitionDim - 1]:
  #       # this axisCode needs improvement
  #       for pair in oneBondPairs:
  #         # First do acquisition dimension
  #         if acquisitionDim in pair:
  #           ll = pair.copy()
  #           ll.remove(acquisitionDim)
  #           otherDim = ll[0]
  #           otherCode = axisCodes[otherDim - 1]
  #           if otherCode == atomTypes[otherDim - 1]:


  def read_nef_spectrum_dimension(self, project:Project, loop:StarIo.NmrLoop):
    """Read nef_spectrum_dimension loop and convert data to a dictionary
    of ccpnTag:[per-dimension-value]"""

    # NB we are not using absolute_peak_positions - if false what could we do?

    result = {}
    nefTag2ccpnTag = nef2CcpnMap['nef_spectrum_dimension']

    rows = sorted(loop.data, key=itemgetter('dimension_id'))
    # rows = [(row['dimension_id'], row) for row in loop.data]
    # rows = [tt[1] for tt in sorted(rows)]

    if not rows:
      raise ValueError("nef_spectrum_dimension is missing or empty")

    # Get spectrum attributes
    for nefTag, ccpnTag in nefTag2ccpnTag.items():
      if nefTag in rows[0]:
        if ccpnTag:
          result[ccpnTag] = list(row.get(nefTag) for row in rows)
        else:
          # Unmapped attributes are passed with the original tag for later processing
          result[nefTag] = list(row.get(nefTag) for row in rows)
    #
    return result


  def load_ccpn_integral_list(self, spectrum:Spectrum,
                              loop:StarIo.NmrLoop) -> List[IntegralList]:

    result = []

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    creatorFunc = spectrum.newIntegralList
    for row in loop.data:
      parameters = self._parametersFromLoopRow(row, map2)
      integralList = creatorFunc(**parameters)
      integralList.resetSerial(row['serial'])
      # NB former call was BROKEN!
      # modelUtil.resetSerial(integralList, row['serial'], 'integralLists')
      result.append(integralList)
    #
    return result


  def load_ccpn_integral(self, spectrum:Spectrum,
                              loop:StarIo.NmrLoop) -> List[Integral]:

    result = []

    # Get name map for per-dimension attributes
    max = spectrum.dimensionCount + 1
    multipleAttributes = {
      'slopes':tuple('slopes_%s' % ii for ii in range(1, max)),
      'lowerLimits':tuple('lower_limits_%s' % ii for ii in range(1, max)),
      'upperLimits':tuple('upper_limits_%s' % ii for ii in range(1, max)),
    }

    serial2creatorFunc = dict((x.serial,x.newIntegral) for x in spectrum.integralLists)

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    creatorFunc = spectrum.newIntegralList
    for row in loop.data:
      parameters = self._parametersFromLoopRow(row, map2)
      integral = serial2creatorFunc[row['integral_list_serial']](**parameters)
      integral.slopes = tuple(row.get(x) for x in multipleAttributes['slopes'])
      lowerLimits = tuple(row.get(x) for x in multipleAttributes['lowerLimits'])
      upperLimits = tuple(row.get(x) for x in multipleAttributes['upperLimits'])
      # integral.slopes = row._get('slopes')
      # lowerLimits = row._get('lower_limits')
      # upperLimits = row._get('upper_limits')
      integral.limits = zip((lowerLimits, upperLimits))
      integral.resetSerial(row['integral_serial'])
      # NB former call was BROKEN!
      # modelUtil.resetSerial(integral, row['integral_serial'], 'integrals')
      result.append(integral)
    #
    return result

  def load_nef_peak(self, peakList:PeakList, loop:StarIo.NmrLoop) -> List[Peak]:
    """Serves to load nef_peak loop"""

    result = []

    dimensionCount = peakList.spectrum.dimensionCount

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])

    # Get name map for per-dimension attributes
    max = dimensionCount + 1
    multipleAttributes = {
      'position':tuple('position_%s' % ii for ii in range(1, max)),
      'positionError':tuple('position_uncertainty_%s' % ii for ii in range(1, max)),
      'chainCodes':tuple('chain_code_%s' % ii for ii in range(1, max)),
      'sequenceCodes':tuple('sequence_code_%s' % ii for ii in range(1, max)),
      'residueTypes':tuple('residue_name_%s' % ii for ii in range(1, max)),
      'atomNames':tuple('atom_name_%s' % ii for ii in range(1, max)),
    }

    peaks = {}
    assignedNmrAtoms = []

    for row in loop.data:

      parameters = self._parametersFromLoopRow(row, map2)

      # get or make peak
      serial = parameters['serial']
      peak = peaks.get(serial)
      # TODO check if peak parameters are the same for all rows, and do something about it
      # For now we simply use the first row that appears
      if peak is None:
        # start of a new peak

        # finalise last peak
        if result and assignedNmrAtoms:
          # There is a peak in result, and the peak has assignments to set
          result[-1].assignedNmrAtoms = assignedNmrAtoms
          assignedNmrAtoms.clear()

        # make new peak  multipleAttributes
        # parameters['position'] = row._get('position')[:dimensionCount]
        parameters['position'] = tuple(row.get(x) for x in multipleAttributes['position'])
        parameters['positionError'] = tuple(row.get(x) for x in multipleAttributes['positionError'])
        # parameters['positionError'] = row._get('position_uncertainty')[:dimensionCount]
        peak = peakList.newPeak(**parameters)
        peaks[serial] = peak
        result.append(peak)

      # Add assignment
      chainCodes = tuple(row.get(x) for x in multipleAttributes['chainCodes'])
      sequenceCodes = tuple(row.get(x) for x in multipleAttributes['sequenceCodes'])
      residueTypes = tuple(row.get(x) for x in multipleAttributes['residueTypes'])
      atomNames = tuple(row.get(x) for x in multipleAttributes['atomNames'])
      # chainCodes = row._get('chain_code')[:dimensionCount]
      # sequenceCodes = row._get('sequence_code')[:dimensionCount]
      # residueTypes = row._get('residue_name')[:dimensionCount]
      # atomNames = row._get('atom_name')[:dimensionCount]
      assignments = zip(chainCodes, sequenceCodes, residueTypes, atomNames)
      nmrAtoms = []
      foundAssignment = False
      for tt in assignments:
        if all(x is None for x in tt):
          # No assignment
          nmrAtoms.append(None)
        elif tt[1] and tt[3]:
          # Enough for an assignment - make it
          foundAssignment = True
          nmrResidue = self.produceNmrResidue(*tt[:3])
          nmrAtom = self.produceNmrAtom(nmrResidue, tt[3])
          nmrAtoms.append(nmrAtom)
        else:
          # partial and unusable assignment
          self.warning("Uninterpretable Peak assignment for peak %s: %s. Set to None"
                       % (peak.serial, tt))
          nmrAtoms.append(None)
      if foundAssignment:
        assignedNmrAtoms.append(nmrAtoms)

    # finalise last peak
    if result and assignedNmrAtoms:
      # There is a peak in result, and the peak has assignments to set
      result[-1].assignedNmrAtoms = assignedNmrAtoms
      assignedNmrAtoms.clear()
    #
    return result
  #
  importers['nef_peak'] = load_nef_peak


  def load_nef_peak_restraint_links(self, project:Project, saveFrame:StarIo.NmrSaveFrame):
    """load nef_peak_restraint_links saveFrame"""
    mapping = nef2CcpnMap['nef_peak_restraint_links']
    for tag, ccpnTag in mapping.items():
      if ccpnTag == _isALoop:
        loop = saveFrame.get(tag)
        if loop:
          importer = self.importers[tag]
          importer(self, project, loop)
    #
    return project
  #
  importers['nef_peak_restraint_links'] = load_nef_peak_restraint_links


  def load_nef_peak_restraint_link(self, project:Project, loop:StarIo.NmrLoop):
    """Load nef_peak_restraint_link loop"""

    links = {}

    # NBNB TODO. There was a very strange bug in this function
    # When I was using PeakList.getPeak(str(serial))
    # and RestraintList.getRestraint(str(serial), peaks and restraints were
    # sometimes missed even though the data were present.
    # Doing the test at tge API level (as now) fixed the problem
    # THIS SHOULD BE IMPOSSIBLE
    # At some point we ought to go back, reproduce the bug, and remove the reason for it.

    for row in loop.data:
      peakList = self.frameCode2Object.get(row.get('nmr_spectrum_id'))
      if peakList is None:
        self.warning(
          "No Spectrum saveframe found with framecode %s. Skipping peak_restraint_link"
          % row.get('nmr_spectrum_id')
        )
        continue
      restraintList = self.frameCode2Object.get(row.get('restraint_list_id'))
      if restraintList is None:
        self.warning(
          "No RestraintList saveframe found with framecode %s. Skipping peak_restraint_link"
          % row.get('restraint_list_id')
        )
        continue
      peak = peakList._wrappedData.findFirstPeak(serial=row.get('peak_id'))
      if peak is None:
        self.warning(
          "No peak %s found in %s Skipping peak_restraint_link"
          % (row.get('peak_id'), row.get('nmr_spectrum_id'))
        )
        continue
      restraint = restraintList._wrappedData.findFirstConstraint(serial=row.get('restraint_id'))
      if restraint is None:
        self.warning(
          "No restraint %s found in %s Skipping peak_restraint_link"
          % (row.get('restraint_id'), row.get('restraint_list_id'))
        )
        continue

      # Is all worked, now accumulate the lin
      ll = links.get(restraint, [])
      ll.append(peak)
      links[restraint] = ll

    # Set the actual links
    for restraint, peaks in links.items():
      restraint.peaks = peaks
    #
    return None
  #
  importers['nef_peak_restraint_link'] = load_nef_peak_restraint_link


  def load_ccpn_spectrum_group(self ,project:Project, saveFrame:StarIo.NmrSaveFrame):

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    # Make main object
    result = project.newSpectrumGroup(**parameters)

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        importer(self, result, loop)
    #
    return result
  #
  importers['ccpn_spectrum_group'] = load_ccpn_spectrum_group

  def load_ccpn_group_spectrum(self, parent:SpectrumGroup, loop:StarIo.NmrLoop):
    """load ccpn_group_spectrum loop"""

    spectra = []
    for row in loop.data:
      peakList = self.frameCode2Object.get(row.get('nmr_spectrum_id'))
      if peakList is None:
        self.warning(
          "No Spectrum saveframe found with framecode %s. Skipping Spectrum from SpectrumGroup"
          % row.get('nmr_spectrum_id')
        )
      else:
        spectra.append(peakList.spectrum)
    #
    parent.spectra = spectra
  #
  importers['ccpn_group_spectrum'] = load_ccpn_group_spectrum


  def load_ccpn_complex(self ,project:Project, saveFrame:StarIo.NmrSaveFrame):

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    # Make main object
    result = project.newComplex(**parameters)

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        importer(self, result, loop)
    #
    return result
  #
  importers['ccpn_complex'] = load_ccpn_complex

  def load_ccpn_complex_chain(self, parent:Complex, loop:StarIo.NmrLoop):
    """load ccpn_complex_chain loop"""

    chains = []
    for row in loop.data:
      chain = self.project.getChain(row.get('complex_chain_code'))
      if chain is None:
        self.warning(
          "No Chain found with code %s. Skipping Chain from Complex"
          % row.get('complex_chain_code')
        )
      else:
        chains.append(chain)
    #
    parent.chains = chains
  #
  importers['ccpn_complex_chain'] = load_ccpn_complex_chain


  def load_ccpn_sample(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    # NBNB TODO add crosslinks to spectrum (also for components)

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]

    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    # Make main object
    result = project.newSample(**parameters)

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        importer(self, result, loop)
    #
    return result
  #
  importers['ccpn_sample'] = load_ccpn_sample


  def load_ccpn_sample_component(self, parent:Sample, loop:StarIo.NmrLoop):
    """load ccpn_sample_component loop"""

    result = []

    creatorFunc = parent.newSampleComponent

    mapping = nef2CcpnMap[loop.name]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    for row in loop.data:
      parameters = self._parametersFromLoopRow(row, map2)
      result.append(creatorFunc(**parameters))
    #
    return result
  #
  importers['ccpn_sample_component'] = load_ccpn_sample_component


  def load_ccpn_substance(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]
    parameters, loopNames = self._parametersFromSaveFrame(saveFrame, mapping)

    name = parameters.pop('name')
    if 'labelling' in parameters:
      labelling = parameters.pop('labelling')
    else:
      labelling = None
    previous = [x for x in project.substances if x.name == name]
    sequence = saveFrame.get('sequence_string')
    if sequence and not previous:
      # We have a 'Molecule' type substance with a sequence and no previous occurrence
      # Create it as new polymer
      if ',' in sequence:
        sequence = list(sequence.split(','))
      params = {'molType':saveFrame.get('mol_type')}
      startNumber = saveFrame.get('start_number')
      if startNumber is not None:
        params['startNumber'] = startNumber
      isCyclic = saveFrame.get('is_cyclic')
      if isCyclic is not None:
        params['isCyclic'] = isCyclic
      #
      result = project.createPolymerSubstance(sequence, name, labelling, **params)

    else:
      # find or create substance
      # NB substance could legitimately be existing already, since substances are created
      # when a chain is created.
      result = project.fetchSubstance(name, labelling)
      if previous:
        # In case this is a new Substance, (known name, different labelling)
        # set the sequenceString, if any, to the same as previous
        sequenceString = previous[0].sequenceString
        if sequenceString is not None:
          result._wrappedData.seqString = sequenceString

    # Whether substance was pre-existing or not
    # overwrite the missing substance-specific parameters
    for tag, val in parameters.items():
      setattr(result, tag, val)

    # Load loops, with object as parent
    for loopName in loopNames:
      loop = saveFrame.get(loopName)
      if loop:
        importer = self.importers[loopName]
        importer(self, result, loop)
  #
  importers['ccpn_substance'] = load_ccpn_substance

  def load_ccpn_substance_synonym(self, parent:Substance, loop:StarIo.NmrLoop):
    """load ccpn_substance_synonym loop"""

    result = [row['synonym'] for row in loop.data]
    parent.synonyms = result
    #
    return result
  #
  importers['ccpn_substance_synonym'] = load_ccpn_substance_synonym

  def load_ccpn_assignment(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    # the saveframe contains nothing but three loops:
    nmrChainLoopName = 'nmr_chain'
    nmrResidueLoopName = 'nmr_residue'
    nmrAtomLoopName = 'nmr_atom'

    nmrChains = {}
    nmrResidues = {}

    # read nmr_chain loop
    mapping = nef2CcpnMap[nmrChainLoopName]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    creatorFunc = project.newNmrChain
    for row in saveFrame[nmrChainLoopName].data:
      parameters = self._parametersFromLoopRow(row, map2)
      nmrChain = creatorFunc(**parameters)
      nmrChain.resetSerial(row['serial'])
      # NB former call was BROKEN!
      # modelUtil.resetSerial(nmrChain, row['serial'], 'nmrChains')
      nmrChains[parameters['shortName']] = nmrChain

    # read nmr_residue loop
    mapping = nef2CcpnMap[nmrResidueLoopName]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    for row in saveFrame[nmrResidueLoopName].data:
      parameters = self._parametersFromLoopRow(row, map2)
      chainCode =  row['chain_code']
      nmrChain = nmrChains[chainCode]
      nmrResidue = nmrChain.newNmrResidue(**parameters)
      nmrResidue.resetSerial(row['serial'])
      # NB former call was BROKEN!
      # modelUtil.resetSerial(nmrResidue, row['serial'], 'nmrResidues')
      nmrResidues[(chainCode,parameters['sequenceCode'])] = nmrChain

    # read nmr_atom loop
    mapping = nef2CcpnMap[nmrAtomLoopName]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    for row in saveFrame[nmrAtomLoopName].data:
      parameters = self._parametersFromLoopRow(row, map2)
      chainCode =  row['chain_code']
      sequenceCode =  row['sequence_code']
      nmrResidue = nmrResidue[(chainCode, sequenceCode)]
      nmrAtom = nmrResidue.newNmrAtom(**parameters)
      nmrAtom.resetSerial(row['serial'])
      # NB former call was BROKEN!
      # modelUtil.resetSerial(nmrAtom, row['serial'], 'nmrAtoms')
  #
  importers['ccpn_assignment'] = load_ccpn_assignment


  def load_ccpn_notes(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    # ccpn_notes contains nothing except for the ccpn_note loop
    loopName = 'ccpn_notes'
    loop = saveFrame[loopName]
    creatorFunc = project.newNote

    result = []
    mapping = nef2CcpnMap[loopName]
    map2 = dict(item for item in mapping.items() if item[1] and '.' not in item[1])
    for row in loop.data:
      parameters = self._parametersFromLoopRow(row, map2)
      result.append(creatorFunc(**parameters))

      # load time stamps and serial = must bypass the API, as they are frozen
      apiNote = result._wrappedData
      created = row.get('created')
      if created:
        apiNote.__dict__['created'] = datetime.strptime(created, Constants.isoTimeFormat)
      lastModified = row.get('last_modified')
      if lastModified:
        apiNote.__dict__['lastModified'] = datetime.strptime(lastModified,
                                                             Constants.isoTimeFormat)
      serial = row.get('serial')
      if serial is not None:
        result.resetSerial(serial)
        # NB former call was BROKEN!
        # modelUtil.resetSerial(apiNote, serial, 'notes')

    #
    return result

  #
  importers['ccpn_notes'] = load_ccpn_notes


  def load_ccpn_additional_data(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    # ccpn_additional_data contains nothing except for the ccpn_internal_data loop
    loopName = 'ccpn_internal_data'
    loop = saveFrame[loopName]
    for row in loop.data:
      pid, data = row.values()
      project.getByPid(pid)._ccpnInternalData = jsonIo.loads(data)
  #
  importers['ccpn_additional_data'] = load_ccpn_additional_data


  def load_ccpn_dataset(self, project:Project, saveFrame:StarIo.NmrSaveFrame):

    print ("ccpn_dataset reading is not implemented yet")

    # Get ccpn-to-nef mappping for saveframe
    category = saveFrame['sf_category']
    framecode = saveFrame['sf_framecode']
    mapping = nef2CcpnMap[category]
  #
  importers['ccpn_dataset'] = load_ccpn_dataset


  def _parametersFromSaveFrame(self, saveFrame:StarIo.NmrSaveFrame, mapping:OD,
                               ccpnPrefix:str=None):
    """Extract {parameter:value} dictionary and list of loop names from saveFrame

    The mapping gives the map from NEF tags to ccpn tags.
    If the ccpn tag is of the form <ccpnPrefix.tag> it is ignored unless
    the first part of teh tag matches the passed-in ccpnPrefix
    (NB the ccpnPrefix may contain '.')."""

    # Get attributes that have a simple tag mapping, and make a separate loop list
    parameters = {}
    loopNames = []
    if ccpnPrefix is None:
      # Normal extraction from saveframe map
      for tag, ccpnTag in mapping.items():
        if ccpnTag == _isALoop:
          loopNames.append(tag)
        elif ccpnTag and '.' not in ccpnTag:
          val = saveFrame.get(tag)
          if val is not None:
            # necessary as tags like ccpn_serial should NOT be set if absent of None
            parameters[ccpnTag] = val
    else:
      # extracting tags of the form `ccpnPrefix`.tag
      for tag, ccpnTag in mapping.items():
        if ccpnTag == _isALoop:
          loopNames.append(tag)
        elif ccpnTag:
          parts = ccpnTag.rsplit('.', 1)
          if parts[0] == ccpnPrefix:
            val = saveFrame.get(tag)
            if val is not None:
              # necessary as tags like ccpn_serial should NOT be set if absent of None
              parameters[parts[1]] = val

    #
    return parameters, loopNames

  def warning(self, message:str):
    template = "WARNING in saveFrame%s\n%s"
    self.warnings.append(template % (self.saveFrameName, message))

  def _parametersFromLoopRow(self, row, mapping):
    parameters = {}
    for tag, ccpnTag in mapping.items():
      val = row.get(tag)
      if val is not None:
        parameters[ccpnTag] = val
    #
    return parameters

  def produceNmrChain(self, chainCode:str=None):
    """Get NmrChain, correcting for possible errors"""

    if chainCode is None:
        chainCode = coreConstants.defaultNmrChainCode
    newChainCode = chainCode
    while True:
      try:
        nmrChain = self.project.fetchNmrChain(newChainCode)
        return nmrChain
      except ValueError:
        newChainCode = '`%s`' % newChainCode
        self.warning("New NmrChain:%s name caused an error.  Renamed %s"
                     % (chainCode, newChainCode))

  def produceNmrResidue(self, chainCode:str=None, sequenceCode:str=None, residueType:str=None):
    """Get NmrResidue, correcting for possible errors"""

    inputTuple = (chainCode, sequenceCode, residueType)
    result = self._nmrResidueMap.get(inputTuple)
    if result is not None:
      return result

    if not sequenceCode:
      raise ValueError("Cannot produce NmrResidue for sequenceCode: %s" % repr(sequenceCode))

    if isinstance(sequenceCode, int):
      sequenceCode = str(sequenceCode)

    nmrChain = self.produceNmrChain(chainCode)

    rt = residueType or ''
    cc = nmrChain.shortName

    try:
      result = nmrChain.fetchNmrResidue(sequenceCode, residueType)
      # return result
    except ValueError:
      # This can happen legitimately, e.g. for offset residues
      # Further processing needed.

      # Parse values from sequenceCode
      seqNo, insertCode, offset = commonUtil.parseSequenceCode(sequenceCode)

      if offset is None:
        if residueType:
          # This could be a case where the NmrResidue had been created from an offset NmrResidue
          # with the residueType therefore missing.
          # If so, set the residueType
          # NBNB this also has the effect of having real entries with empty residueType
          # overridden by clashing entries with residueType set,
          # but that does make some sense and anyway cannot be helped.

          # NB this must be a pre-existing residue or it would have been created above.
          try:
            previous = nmrChain.fetchNmrResidue(sequenceCode)
            if not previous.residueType:
              result = previous
              result._wrappedData.residueType = residueType
              # return result
          except ValueError:
            # Deal with it below
            pass

      else:
        # Offset NmrResidue - get mainNmrResidue
        ss = '%+d' % offset
        try:
          # Fetch the mainNmrResidue. This will create it if it is missing (if possible)
          mainNmrResidue =  nmrChain.fetchNmrResidue(sequenceCode[:-len(ss)])
        except ValueError:
          # No go
          mainNmrResidue = None

        if mainNmrResidue is not None:
          try:
            result = nmrChain.fetchNmrResidue(sequenceCode, residueType)
            # return result
          except ValueError:
            # Handle lower down
            pass

    # If e get here, we could not make an NmrResidue that matched the input
    # Make with modified sequenceCode
    if result is None:
      newSequenceCode = sequenceCode
      while True:
        try:
          result = nmrChain.fetchNmrResidue(newSequenceCode, residueType)
          return result
        except ValueError:
          newSequenceCode = '`%s`' % newSequenceCode
          self.warning("New NmrResidue:%s.%s.%s name caused an error.  Renamed %s.%s.%s"
                       % (cc, sequenceCode, rt, cc, newSequenceCode, rt))
    #
    # Result cannot have been in map, so put it there
    self._nmrResidueMap[inputTuple] = result
    return result

  def produceNmrAtom(self, nmrResidue:NmrResidue, name:str, isotopeCode=None):
    """Get NmrAtom from NmrResidue and name, correcting for possible errors"""

    if not name:
      raise ValueError("Cannot produce NmrAtom for atom name: %s" % repr(name))

    if isotopeCode:
      prefix = isotopeCode.upper()
      if not name.startswith(prefix):
        prefix = commonUtil.isotopeCode2Nucleus(isotopeCode)
        if prefix is None:
          self.warning("Ignoring unsupported isotopeCode: %s for NmrAtom:%s.%s"
                       % (isotopeCode, nmrResidue._id, name))
        elif not name.startswith(prefix):
          newName = '%s@%s' % (prefix, name)
          self.warning("NmrAtom name %s does not match isotopeCode %s - renamed to %s"
                       % (isotopeCode, name, newName))
          name = newName

    newName = name
    while True:
      nmrAtom = nmrResidue.getNmrAtom(newName.translate(Pid.remapSeparators))
      if nmrAtom is None:
        try:
          nmrAtom = nmrResidue.newNmrAtom(newName, isotopeCode)
          return nmrAtom
        except ValueError:
          pass
      elif isotopeCode in (None, nmrAtom.isotopeCode):
        # We must ensure that the isotopeCodes match
        return nmrAtom
      else:
        # IsotopeCode mismatch. Try to
        if prefix == isotopeCode.upper():
          # Something wrong here.
          raise ValueError("Clash between NmrAtom %s (%s) and %s (%s) in %s"
                           % (nmrAtom.name, nmrAtom.isotopeCode, newName, isotopeCode, nmrResidue))
        else:
          newName = isotopeCode.upper() + newName[len(prefix):]
          continue

      # If we get here there was an error. Change name and try again
      tt = newName.rsplit('_', 1)
      if len(tt) == 2 and tt[1].isdigit():
        newName = '%s_%s' % (tt[0], int(tt[1]) + 1)
      else:
        newName += '_1'
      self.warning("New NmrAtom:%s.%s name caused an error.  Renamed %s.%s"
                   % (nmrResidue._id, name, nmrResidue._id, newName))

  def updateMetaData(self, metaDataFrame:StarIo.NmrSaveFrame):
    """Add meta information to main data set. Must be done at end of read"""

    # NBNB NOT WORKING YET!

    # dataSet = self.fetchDataSet(self.mainDataSetSerial)
    self.mainDataSetSerial = None

  def fetchDataSet(self, serial:int=None):
    """Fetch DataSet with given serial.
    If input is None, use self.defaultDataSetSerial
    If that too is None, create a new DataSet and use its serial as the default

    NB when reading, all DataSets with known serials should be instantiated BEFORE calling
    with input None"""

    if serial is None:
      serial = self.defaultDataSetSerial

    if serial is None:
      # default not set - create one
      dataSet = self.project.newDataSet()
      serial = dataSet.serial
      dataSet.title = 'Data_%s' % serial
      self.defaultDataSetSerial = serial
      self._dataSet2ItemMap[dataSet] = dataSet._getTempItemMap()

    else:
      # take or create dataSet matching serial
      dataSet = self.project.getDataSet(str(serial))
      if dataSet is None:
        dataSet = self.project.newDataSet()
        dataSet.resetSerial(serial)
        # NB former call was BROKEN!
        # modelUtil.resetSerial(dataSet._wrappedData, serial, 'nmrConstraintStores')
        dataSet._finaliseAction('rename')
        dataSet.title = 'Data_%s' % serial

        self._dataSet2ItemMap[dataSet] = dataSet._getTempItemMap()
    #
    return dataSet


def createSpectrum(project:Project, spectrumName:str, spectrumParameters:dict,
                   dimensionData:dict, transferData:Sequence[Tuple]=None):
  """Get or create spectrum using dictionaries of attributes, such as read in from NEF.

  :param spectrumParameters keyword-value dictionary of attribute to set on resulting spectrum

  :params Dictionary of keyword:list parameters, with per-dimension parameters.
  Either 'axisCodes' or 'isotopeCodes' must be present and fully populated.
  A number of other dimensionData are
  treated specially (see below)
  """

  spectrum = project.getSpectrum(spectrumName)
  if spectrum is None:
    # Spectrum did not already exist

    dimTags = list(dimensionData.keys())

    # First try to load it - we override the loaded attribute values below
    # but loading gives a more complete parameter set.
    spectrum = None
    filePath = spectrumParameters.get('filePath')
    if filePath and os.path.exists(filePath):
      try:
        dataType, subType, usePath = ioFormats.analyseUrl(path)
        if dataType == 'Spectrum':
          spectra = project.loadSpectrum(usePath, subType)
          if spectra:
            spectrum = spectra[0]
      except:
        # Deliberate - any error should be skipped
        pass
      if spectrum is None:
        project._logger.warning("Failed to load spectrum from spectrum path %s" % filePath)
      elif 'axisCodes' in dimensionData:
          # set axisCodes
          spectrum.axisCodes = dimensionData['axisCodes']

    acquisitionAxisIndex = None
    if 'is_acquisition' in dimensionData:
      dimTags.remove('is_acquisition')
      values = dimensionData['is_acquisition']
      if values.count(True) == 1:
        acquisitionAxisIndex = values.index(True)

    if spectrum is None:
      # Spectrum could not be loaded - now create a dummy spectrum

      if 'axisCodes' in dimTags:
        # We have the axisCodes, from ccpn
        dimTags.remove('axisCodes')
        axisCodes = dimensionData['axisCodes']

      else:
        if transferData is None:
          raise ValueError("Function needs either axisCodes or transferData")

        dimensionIds = dimensionData['dimension_id']

        # axisCodes were not set - produce a serviceable set
        axisCodes = makeNefAxisCodes(isotopeCodes=dimensionData['isotopeCodes'],
                                     dimensionIds=dimensionIds,
                                     acquisitionAxisIndex=acquisitionAxisIndex,
                                     transferData=transferData)


      # make new spectrum with default parameters
      spectrum = project.createDummySpectrum(axisCodes, spectrumName,
        chemicalShiftList=spectrumParameters.get('chemicalShiftList')
      )
      if acquisitionAxisIndex is not None:
        spectrum.acquisitionAxisCode = axisCodes[acquisitionAxisIndex]

      # Delete autocreated peaklist  and reset - we want any read-in peakList to be the first
      # If necessary an empty PeakList is added downstream
      spectrum.peakLists[0].delete()
      spectrum._wrappedData.__dict__['_serialDict']['peakLists'] = 0

    # (Re)set all spectrum attributes

    # First per-dimension ones
    dimTags.remove('dimension_id')
    if 'absolute_peak_positions' in dimensionData:
      # NB We are not using these. What could we do with them?
      dimTags.remove('absolute_peak_positions')
    if 'folding' in dimensionData:
      dimTags.remove('folding')
      values = [None if x == 'none' else x for x in dimensionData['folding']]
      spectrum.foldingModes = values
    if 'pointCounts' in dimensionData:
      dimTags.remove('pointCounts')
      spectrum.pointCounts = pointCounts = dimensionData['pointCounts']
      if 'totalPointCounts' in dimensionData:
        dimTags.remove('totalPointCounts')
        spectrum.totalPointCounts = dimensionData['totalPointCounts']
      else:
        spectrum.totalPointCounts = pointCounts
    # Needed below:
    if 'value_first_point' in dimensionData:
      dimTags.remove('value_first_point')
    if 'referencePoints' in dimensionData:
      dimTags.remove('referencePoints')
    # value_first_point = dimensionData.get('value_first_point')
    # if value_first_point is not None:
    #   dimensionData.pop('value_first_point')
    # referencePoints = dimensionData.get('referencePoints')
    # if referencePoints is not None:
    #   dimensionData.pop('referencePoints')

    # Remaining per-dimension values match the spectrum. Set them.
    # NB we use the old (default) values where the new value is None
    # - some attributes like spectralWidths do not accept None.
    if 'spectrometerFrequencies' in dimTags:
      # spectrometerFrequencies MUST be set before spectralWidths,
      # as the spectralWidths are otherwise modified
      dimTags.remove('spectrometerFrequencies')
      dimTags.insert(0, 'spectrometerFrequencies')
    for tag in dimTags:
      vals = dimensionData[tag]
      # Use old values where new ones are None
      oldVals = getattr(spectrum, tag)
      vals = [x if x is not None else oldVals[ii] for ii,x in enumerate(vals)]
      setattr(spectrum, tag, vals)

    # Set referencing.
    value_first_point = dimensionData.get('value_first_point')
    referencePoints = dimensionData.get('referencePoints')
    if value_first_point is None:
      # If reading NEF we must get value_first_point,
      #  but in other uses we might be getting referencePoints, referenceValues directly
      referenceValues = dimensionData.get('referenceValues')
      if referenceValues and referencePoints:
        spectrum.referencePoints = referencePoints
        spectrum.referenceValues = referenceValues
    else:
      if referencePoints is None:
        # not CCPN data
        referenceValues = spectrum.referenceValues
        for ii, val in enumerate(value_first_point):
          if val is None:
            value_first_point[ii] = referenceValues[ii]
        spectrum.referenceValues = value_first_point
        spectrum.referencePoints = [1] * len(referenceValues)
      else:
        points = list(spectrum.referencePoints)
        values = list(spectrum.referenceValues)
        sw = spectrum.spectralWidths
        pointCounts = spectrum.pointCounts
        for ii, refVal in enumerate(value_first_point):
          refPoint = referencePoints[ii]
          if refVal is not None and refPoint is not None:
            # if we are here refPoint should never be None, but OK, ...
            # Set reference to use refPoint
            points[ii] = refPoint
            refVal -= ((refPoint-1) * sw[ii] / pointCounts[ii])
            values[ii] = refVal
        spectrum.referencePoints = points
        spectrum.referenceValues = values

    # Then spectrum-level ones
    for tag, val in spectrumParameters.items():
      if tag != 'dimensionCount':
        # dimensionCount is handled already and not settable
        setattr(spectrum, tag, val)
    #
    return spectrum

  else:
    raise ValueError("Spectrum named %s already exists" % spectrumName)

def makeNefAxisCodes(isotopeCodes:Sequence[str], dimensionIds:List[int],
                     acquisitionAxisIndex:int, transferData:Sequence[Tuple]):

  nuclei = [commonUtil.splitIntFromChars(x)[1] for x in isotopeCodes]
  dimensionToNucleus = dict((zip(dimensionIds, nuclei)))
  dimensionToAxisCode = dimensionToNucleus.copy()

  oneBondConnections = {}
  for startNuc in 'FH':
    # look for onebond to F or H, the latter taking priority
    for dim1, dim2, transferType, isIndirect in transferData:
      if transferType == 'onebond':
        nuc1, nuc2 = dimensionToNucleus[dim1], dimensionToNucleus[dim2]
        if startNuc in (nuc1, nuc2):
          if startNuc == nuc1:
            oneBondConnections[dim1] = dim2
          else:
            oneBondConnections[dim2] = dim1
          dimensionToAxisCode[dim1] = nuc1 + nuc2.lower()
          dimensionToAxisCode[dim2] = nuc2 + nuc1.lower()

  resultMap = {}
  acquisitionAtEnd = False
  if acquisitionAxisIndex is not None:
    acquisitionAtEnd = acquisitionAxisIndex >= 0.5*len(isotopeCodes)
    # if acquisitionAtEnd:
    #   # reverse, because acquisition end of dimensions should
    #   # be the one WITHOUT number suffixes
    #   dimensionIds.reverse()

    # Put acquisition axis first, to make sure it gets the lowest number
    # even if it is not teh first to start with.
    acquisitionAxisId = dimensionIds[acquisitionAxisIndex]
    dimensionIds.remove(acquisitionAxisId)
    dimensionIds.insert(0, acquisitionAxisId)
  for dim in dimensionIds:
    axisCode = dimensionToAxisCode[dim]
    if axisCode in resultMap.values():
      ii = 0
      ss = axisCode
      while ss in resultMap.values():
        ii += 1
        ss = '%s%s' % (axisCode, ii)
      otherDim = oneBondConnections.get(dim)
      if otherDim is not None:
        # We are adding a suffix to e.g. Hc. Add the same suffix to equivalent Ch
        # NB this should only happen for certain 4D experiments.
        # NB not well tested, but better than leaving in a known error.
        ss = '%s%s' % (dimensionToAxisCode[otherDim], ii)
        if otherDim < dim:
          resultMap[otherDim] = ss
        dimensionToAxisCode[otherDim] = ss

    resultMap[dim] = axisCode
  dimensionIds.sort()
  # NBNB new attempt - may not work
  # if acquisitionAtEnd:
  #   # put result back in dimension order
  #   dimensionIds.reverse()
  result = list(resultMap[ii] for ii in dimensionIds)
  #
  return result