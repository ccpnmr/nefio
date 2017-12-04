"""
NefImporter - a series of routines for reading a Nef file and examining the contents.

Module Contents
===============

NefImporter contains:

  loadFile            read in the contents of a Nef file
  getCategories       return the current categories defined in the Nef structure
  getSaveFrameNames   return the names of the saveFrames with the file

  getChemicalShiftLists
  get<name>                   return the relevant structures of the Nef file
                              defined by the available categories

  getSaveFrame                return saveFrame of the given name
    sf.getTable               return table name form the saveFrame
    sf.hasTable               check if table exists
    sf.setTable               set the table



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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this is a fix to get the import to work when running the test case
import sys, importlib
from pathlib import Path

def import_parents(level=1):
  global __package__
  file = Path(__file__).resolve()
  parent, top = file.parent, file.parents[level]

  sys.path.append(str(top))
  try:
    sys.path.remove(str(parent))
  except ValueError:  # already removed
    pass

  __package__ = '.'.join(parent.parts[len(top.parts):])
  importlib.import_module(__package__)  # won't be needed after that

if __name__ == '__main__' and __package__ is None:
  import_parents()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MAJOR_VERSION = '0'
MINOR_VERSION = '8'
PATCH_LEVEL = '4'
__nef_version__ = '.'.join( (MAJOR_VERSION, MINOR_VERSION) )
# __version__ = '.'.join( (__nef_version__, PATCH_LEVEL) )


NEF_CATEGORIES = ['nef_nmr_meta_data',
                  'nef_molecular_system',
                  'nef_chemical_shift_list',
                  'nef_distance_restraint_list',
                  'nef_dihedral_restraint_list',
                  'nef_rdc_restraint_list',
                  'nef_nmr_spectrum',
                  'nef_peak_restraint_links']

NEF_REQUIRED_SAVEFRAME_BY_FRAMECODE = ['nef_nmr_meta_data',
                                       'nef_molecular_system']
NEF_REQUIRED_SAVEFRAME_BY_CATEGORY = ['nef_chemical_shift_list', ]

NEF_ALL_SAVEFRAME_REQUIRED_FIELDS = ['sf_category',
                                     'sf_framecode', ]

MD_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode',
                      'format_name',
                      'format_version',
                      'program_name',
                      'program_version',
                      'creation_date',
                      'uuid']
MD_OPTIONAL_FIELDS = ['coordinate_file_name', ]
MD_OPTIONAL_LOOPS = ['nef_related_entries',
                     'nef_program_script',
                     'nef_run_history']
MD_RE_REQUIRED_FIELDS = ['database_name',
                         'database_accession_code']
MD_PS_REQUIRED_FIELDS = ['program_name', ]
MD_RH_REQUIRED_FIELDS = ['run_ordinal',
                         'program_name']
MD_RH_OPTIONAL_FIELDS = ['program_version',
                         'script_name',
                         'script']

MS_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode']
MS_REQUIRED_LOOPS = ['nef_sequence']
MS_OPTIONAL_LOOPS = ['nef_covalent_links']
MS_NS_REQUIRED_FIELDS = ['chain_code',
                         'sequence_code',
                         'residue_type',
                         'linking',
                         'residue_variant']
MS_CL_REQUIRED_FIELDS = ['chain_code_1',
                         'sequence_code_1',
                         'residue_type_1',
                         'atom_name_1',
                         'chain_code_2',
                         'sequence_code_2',
                         'residue_type_2',
                         'atom_name_2']

CSL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'atom_chem_shift_units']
CSL_REQUIRED_LOOPS = ['nef_chemical_shift']
CSL_CS_REQUIRED_FIELDS = ['chain_code',
                          'sequence_code',
                          'residue_type',
                          'atom_name',
                          'value']
CSL_CS_OPTIONAL_FIELDS = ['value_uncertainty', ]

DRL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'potential_type']
DRL_REQUIRED_LOOPS = ['nef_distance_restraint']
DRL_OPTIONAL_FIELDS = ['restraint_origin', ]
DRL_DR_REQUIRED_FIELDS = ['ordinal',
                          'restraint_id',
                          'chain_code_1',
                          'sequence_code_1',
                          'residue_type_1',
                          'atom_name_1',
                          'chain_code_2',
                          'sequence_code_2',
                          'residue_type_2',
                          'atom_name_2',
                          'weight']
DRL_DR_OPTIONAL_FIELDS = ['restraint_combination_id',
                          'target_value',
                          'target_value_uncertainty',
                          'lower_linear_limit',
                          'lower_limit',
                          'upper_limit',
                          'upper_linear_limit']

DIHRL_REQUIRED_FIELDS = ['sf_category',
                         'sf_framecode',
                         'potential_type']
DIHRL_REQUIRED_LOOPS = ['nef_dihedral_restraint']
DIHRL_OPTIONAL_FIELDS = ['restraint_origin', ]
DIHRL_DIHR_REQUIRED_FIELDS = ['ordinal',
                              'restraint_id',
                              'restraint_combination_id',
                              'chain_code_1',
                              'sequence_code_1',
                              'residue_type_1',
                              'atom_name_1',
                              'chain_code_2',
                              'sequence_code_2',
                              'residue_type_2',
                              'atom_name_2',
                              'chain_code_3',
                              'sequence_code_3',
                              'residue_type_3',
                              'atom_name_3',
                              'chain_code_4',
                              'sequence_code_4',
                              'residue_type_4',
                              'atom_name_4',
                              'weight']
DIHRL_DIHR_OPTIONAL_FIELDS = ['target_value',
                              'target_value_uncertainty',
                              'lower_linear_limit',
                              'lower_limit',
                              'upper_limit',
                              'upper_linear_limit',
                              'name']

RRL_REQUIRED_FIELDS = ['sf_category',
                       'sf_framecode',
                       'potential_type']
RRL_REQUIRED_LOOPS = ['nef_rdc_restraint']
RRL_OPTIONAL_FIELDS = ['restraint_origin',
                       'tensor_magnitude',
                       'tensor_rhombicity',
                       'tensor_chain_code',
                       'tensor_sequence_code',
                       'tensor_residue_type', ]
RRL_RR_REQUIRED_FIELDS = ['ordinal',
                          'restraint_id',
                          'chain_code_1',
                          'sequence_code_1',
                          'residue_type_1',
                          'atom_name_1',
                          'chain_code_2',
                          'sequence_code_2',
                          'residue_type_2',
                          'atom_name_2',
                          'weight']
RRL_RR_OPTIONAL_FIELDS = ['restraint_combination_id',
                          'target_value',
                          'target_value_uncertainty',
                          'lower_linear_limit',
                          'lower_limit',
                          'upper_limit',
                          'upper_linear_limit',
                          'scale',
                          'distance_dependent', ]

PL_REQUIRED_FIELDS = ['sf_category',
                      'sf_framecode',
                      'num_dimensions',
                      'chemical_shift_list']
PL_REQUIRED_LOOPS = ['nef_spectrum_dimension',
                     'nef_spectrum_dimension_transfer',
                     'nef_peak']
PL_OPTIONAL_FIELDS = ['experiment_classification',
                      'experiment_type']
PL_SD_REQUIRED_FIELDS = ['dimension_id',
                         'axis_unit',
                         'axis_code']
PL_SD_OPTIONAL_FIELDS = ['spectrometer_frequency',
                         'spectral_width',
                         'value_first_point',
                         'folding',
                         'absolute_peak_positions',
                         'is_acquisition', ]
PL_SDT_REQUIRED_FIELDS = ['dimension_1',
                          'dimension_2',
                          'transfer_type']
PL_SDT_OPTIONAL_FIELDS = ['is_indirect', ]
PL_P_REQUIRED_FIELDS = ['ordinal',
                        'peak_id']
PL_P_REQUIRED_ALTERNATE_FIELDS = [['height', 'volume'], ]
PL_P_REQUIRED_FIELDS_PATTERN = ['position_{}',
                                'chain_code_{}',
                                'sequence_code_{}',
                                'residue_type_{}',
                                'atom_name_{}', ]
PL_P_OPTIONAL_ALTERNATE_FIELDS = {r'(height)':['{}_uncertainty', ],
                                  r'(volume)':['{}_uncertainty', ],
                                  r'position_([0-9]+)':['position_uncertainty_{}', ],
                                  }
PL_P_OPTIONAL_FIELDS_PATTERN = ['position_uncertainty_{}', ]

PRLS_REQUIRED_FIELDS = ['sf_category',
                        'sf_framecode']
PRLS_REQUIRED_LOOPS = ['nef_peak_restraint_link']
PRLS_PRL_REQUIRED_FIELDS = ['nmr_spectrum_id',
                            'peak_id',
                            'restraint_list_id',
                            'restraint_id']

NEF_RETURNALL = 'all'
NEF_RETURNNEF = 'nef_'
NEF_RETURNOTHER = 'other'
NEF_PREFIX = 'nef_'


from . import GenericStarParser
from . import StarIo


class NefDict(StarIo.NmrSaveFrame):
  # class to add functions to a saveFrame
  def __init__(self, inDict):
    super(NefDict, self).__init__(name=inDict.name)
    self._nefDict = inDict

  def getTable(self, name):
    # return table 'name' if exists else None
    pass

  def hasTable(self, name):
    # return True is the table exists in the saveFrame
    pass

  def setTable(self, name):
    # add the table 'name' to the saveFrame, or replace the existing
    # does this need to be here or in the main class?
    pass

class NefImporter():
  """Top level data block for accessing object tree"""
  # put functions in here to read the contents of the dict.
  # superclassed from DataBlock which is of type StarContainer
  def __init__(self, name=None, programName='Unknown', programVersion='Unknown', project=None, initialise=True):

    # import inspect

    # keep a copy of the original dataExtent
    self.name = name
    self._nefDict = StarIo.NmrDataBlock()     # empty block

    # these should all be wrapped
    # # property_names = [p for p in dir(GenericStarParser.DataBlock) if isinstance(getattr(GenericStarParser.DataBlock, p), property)]
    # methodNames = inspect.getmembers(inDict, predicate=inspect.ismethod)
    #
    # # add the method to point to our loaded dataExtent
    # for met in methodNames:
    #   setattr(self.__class__, met[0], met[1])

    # define a new Nef structure
    # if project is not None then build a new Nef structure from the Ccpn project

    if project:
      # new project here - program and version name from project
      pass

    else:
      # new, empty project, set names as above
      pass

    if initialise:
      self.initialise()

  def initialise(self):
    self._nefDict['nef_nmr_meta_data'] = StarIo.NmrDataBlock()
    self._nefDict['nef_nmr_meta_data'].update({k:'' for k in MD_REQUIRED_FIELDS})
    self._nefDict['nef_nmr_meta_data']['sf_category'] = 'nef_nmr_meta_data'
    self._nefDict['nef_nmr_meta_data']['sf_framecode'] = 'nef_nmr_meta_data'
    self._nefDict['nef_nmr_meta_data']['format_name'] = 'Nmr_Exchange_Format'
    self._nefDict['nef_nmr_meta_data']['format_version'] = __nef_version__
    # for l in Nef.MD_REQUIRED_LOOPS:
    #     self['nef_nmr_meta_data'][l] = []

    self._nefDict['nef_molecular_system'] = StarIo.NmrDataBlock()
    self._nefDict['nef_molecular_system'].update({k:'' for k in MS_REQUIRED_FIELDS})
    self._nefDict['nef_molecular_system']['sf_category'] = 'nef_molecular_system'
    self._nefDict['nef_molecular_system']['sf_framecode'] = 'nef_molecular_system'
    for l in MS_REQUIRED_LOOPS:
      self._nefDict['nef_molecular_system'][l] = []

      self.add_chemical_shift_list('nef_chemical_shift_list_1', 'ppm')

  def addSaveFrame(self, name, category, required_fields=None, required_loops=None):
    self._nefDict[name] = StarIo.NmrSaveFrame()
    if required_fields is not None:
      self._nefDict[name].update({k:'' for k in required_fields})
      self._nefDict[name]['sf_category'] = category
      self._nefDict[name]['sf_framecode'] = name
    if required_loops is not None:
      for l in required_loops:
        self._nefDict[name][l] = []
    return self._nefDict[name]

  def add_chemical_shift_list(self, name, cs_units='ppm'):
      category = 'nef_chemical_shift_list'
      self.addSaveFrame(name=name, category=category,
                         required_fields=CSL_REQUIRED_FIELDS,
                         required_loops=CSL_REQUIRED_LOOPS)
      self._nefDict[name]['atom_chem_shift_units'] = cs_units
      return self._nefDict[name]

  def add_distance_restraint_list(self, name, potential_type,
                                  restraint_origin=None):
      category = 'nef_distance_restraint_list'
      self.addSaveFrame(name=name, category=category,
                         required_fields=DRL_REQUIRED_FIELDS,
                         required_loops=DRL_REQUIRED_LOOPS)
      self._nefDict[name]['potential_type'] = potential_type
      if restraint_origin is not None:
        self._nefDict[name]['restraint_origin'] = restraint_origin

      return self._nefDict[name]

  def add_dihedral_restraint_list(self, name, potential_type,
                                  restraint_origin=None):
      category = 'nef_dihedral_restraint_list'
      self.addSaveFrame(name=name, category=category,
                         required_fields=DIHRL_REQUIRED_FIELDS,
                         required_loops=DIHRL_REQUIRED_LOOPS)
      self._nefDict[name]['potential_type'] = potential_type
      if restraint_origin is not None:
        self._nefDict[name]['restraint_origin'] = restraint_origin

      return self._nefDict[name]

  def add_rdc_restraint_list(self, name, potential_type,
                             restraint_origin=None, tensor_magnitude=None,
                             tensor_rhombicity=None, tensor_chain_code=None,
                             tensor_sequence_code=None, tensor_residue_type=None):
      category = 'nef_rdc_restraint_list'
      self.addSaveFrame(name=name, category=category,
                         required_fields=DIHRL_REQUIRED_FIELDS,
                         required_loops=RRL_REQUIRED_LOOPS)
      self._nefDict[name]['potential_type'] = potential_type
      if restraint_origin is not None:
        self._nefDict[name]['restraint_origin'] = restraint_origin
      if tensor_magnitude is not None:
        self._nefDict[name]['tensor_magnitude'] = tensor_magnitude
      if tensor_rhombicity is not None:
        self._nefDict[name]['tensor_rhombicity'] = tensor_rhombicity
      if tensor_chain_code is not None:
        self._nefDict[name]['tensor_chain_code'] = tensor_chain_code
      if tensor_sequence_code is not None:
        self._nefDict[name]['tensor_sequence_code'] = tensor_sequence_code
      if tensor_residue_type is not None:
        self._nefDict[name]['tensor_residue_type'] = tensor_residue_type

      return self._nefDict[name]

  def add_peak_list(self, name, num_dimensions, chemical_shift_list,
                    experiment_classification=None,
                    experiment_type=None):
      category = 'nef_nmr_spectrum'
      if chemical_shift_list in self:
          if self._nefDict[chemical_shift_list]['sf_category'] == 'nef_chemical_shift_list':
              self.addSaveFrame(name=name, category=category,
                                 required_fields=PL_REQUIRED_FIELDS,
                                 required_loops=PL_REQUIRED_LOOPS)
              self._nefDict[name]['num_dimensions'] = num_dimensions
              self._nefDict[name]['chemical_shift_list'] = chemical_shift_list
              if experiment_classification is not None:
                self._nefDict[name]['experiment_classification'] = experiment_classification
              if experiment_type is not None:
                self._nefDict[name]['experiment_type'] = experiment_type

              return self._nefDict[name]
          raise Exception('{} is not a nef_chemical_shift_list.'.format(chemical_shift_list))
      raise Exception('{} does not exist.'.format(chemical_shift_list))

  def add_linkage_table(self):
      name = category = 'nef_peak_restraint_links'
      return self.addSaveFrame(name=name, category=category,
                                required_fields=PRLS_REQUIRED_FIELDS,
                                required_loops=PL_REQUIRED_LOOPS)

  def toString(self):
    try:
      return self._nefDict.toString()
    except:
      return None

  # should this be static
  def fromString(self, text, strict=True):
    # set the Nef from the contents of the string, opposite of toString
    dataExtent = StarIo.parseNef(text=text)
    if dataExtent:
      dbs = [dataExtent[db] for db in dataExtent.keys()]
      if dbs:
        self._nefDict = dbs[0]
    else:
      self._nefDict = None

  def loadFile(self, fileName=None, mode='standard'):
    try:
      nefDataExtent = StarIo.parseNefFile(fileName=fileName, mode=mode)
      self._nefDict = list(nefDataExtent.values())
      if len(self._nefDict) > 1:
        print(
          'More than one datablock in a NEF file is not allowed.  Using the first and discarding the rest.')
      self._nefDict = self._nefDict[0]

      return True
    except:
      return False        # trap any errors and return False

  def saveFile(self, fileName=None, mode='standard'):
    try:
      with open(fileName, 'w') as op:
        op.write(self._nefDict.toString())

      return True
    except:
      return False        # trap any errors and return False

  def getCategories(self):
    # return a list of the categories available in a Nef file
    return tuple(NEF_CATEGORIES)

  def getSaveFrameNames(self, filter=NEF_RETURNALL):
    # return a list of the saveFrames in the file
    names = [self._nefDict[db].name for db in self._nefDict.keys()
             if isinstance(self._nefDict[db], StarIo.NmrSaveFrame)]

    if filter == NEF_RETURNNEF:
      names = [nm for nm in names if nm.startswith(NEF_PREFIX)]
    elif filter == NEF_RETURNOTHER:
      names = [nm for nm in names if not nm.startswith(NEF_PREFIX)]

    return tuple(names)

  def getSaveFrame(self, sfName):
    # return the saveFrame 'name'
    if sfName in self._nefDict:
      return NefDict(self._nefDict[sfName])

    return None

if __name__ == '__main__':
  from ccpn.ui.gui.widgets.Application import TestApplication

  app = TestApplication()

  test = NefImporter()
  test.loadFile('/Users/ejb66/PycharmProjects/Sec5Part3.nef')

  print (test.getCategories())
  names = test.getSaveFrameNames()
  print(names)
  names = test.getSaveFrameNames(filter=NEF_RETURNALL)
  print(names)
  names = test.getSaveFrameNames(filter=NEF_RETURNNEF)
  print (names)
  names = test.getSaveFrameNames(filter=NEF_RETURNOTHER)
  print (names)

  sf1 = test.getSaveFrame(names[0])
  sf2 = test.getSaveFrame('notFound')

  ts = test.toString()
  test.fromString(ts)