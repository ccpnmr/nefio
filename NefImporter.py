"""
NefImporter - a series of routines for reading a Nef file and examining the contents.

Module Contents
===============

NefImporter consists of two classes: NefImporter - a class for handling the top-level object, and
NefDict for handling individual saveFrames in the dictionary.

NefImporter contains:

  initialise          initialise a new dictionary
  loadFile            read in the contents of a .nef file
  saveFile            save the dictionary to a .nef file

  getCategories       return the current categories defined in the Nef structure
  getSaveFrameNames   return the names of the saveFrames with the file
  hasSaveFrame        return True if the saveFrame exists
  getSaveFrame        return saveFrame of the given name
  addSaveFrame        add a new saveFrame to the dictionary


  get<name>           return the relevant structures of the Nef file
                      defined by the available categories, where <name> can be:

                          NmrMetaData
                          MolecularSystems
                          ChemicalShiftLists
                          DistanceRestraintLists
                          DihedralRestraintLists
                          RdcRestraintLists
                          NmrSpectra
                          PeakRestraintLinks

                      e.g. getChemicalShiftLists

  add<name>           add a new saveFrame to the dictionary
                      defined by the available categories, where <name> can be:

                          ChemicalShiftList
                          DistanceRestraintList
                          DihedralRestraintList
                          RdcRestraintList
                          NmrSpectra
                          PeakLists
                          LinkageTables

                      e.g. addChemicalShiftList

  toString            convert Nef dictionary to a string that can be written to a file
  fromString          convert string to Nef dictionary

  getAttributeNames   get a list of the attributes attached to the dictionary
  getAttribute        return the value of the attribute
  hasAttribute        return True if the attribute Exists

  lastError           error code of the last operation
  lastErrorString     error string of the last operation


NefDict contains handling routines:

  getTableNames   return a list of the tables in the saveFrame
  getTable        return table from the saveFrame, it can be returned as an OrderedDict
                  or as a Pandas DataFrame
  hasTable        return true of the table exists
  setTable        set the table - currently not implemented

  getAttributeNames   get a list of the attributes attached to the saveFrame
  getAttribute        return the value of the attribute
  hasAttribute        return True if the attribute Exists

  lastError           error code of the last operation
  lastErrorString     error string of the last operation


Error handling
--------------

Errors can be handled in three different modes:

  'silent'            errors are handled internally and can be interrogated with saveFrame.lastError
                      with no logging to the stderr

  'standard'          errors are handled internally, error messages are logged to stderr.

  'strict'            errors message are logged to stderr and errors are raised to be trapped by
                      the calling functions

  error handling mode can be set at the instantiation of the object, e.g.

    newObject = NefImporter(errorLogging='standard')


Nef File Contents
-----------------

More details of the contents of Nef files can be found in GenericStarParser
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
from functools import wraps

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



import numpy as np
import pandas as pd
from . import GenericStarParser
from . import StarIo
from collections import OrderedDict


MAJOR_VERSION = '0'
MINOR_VERSION = '8'
PATCH_LEVEL = '4'
__nef_version__ = '.'.join( (MAJOR_VERSION, MINOR_VERSION) )
# __version__ = '.'.join( (__nef_version__, PATCH_LEVEL) )


NEF_CATEGORIES = [('nef_nmr_meta_data', 'get_nmr_meta_data'),
                  ('nef_molecular_system', 'get_molecular_systems'),
                  ('nef_chemical_shift_list', 'get_chemical_shift_lists'),
                  ('nef_distance_restraint_list', 'get_distance_restraint_lists'),
                  ('nef_dihedral_restraint_list', 'get_dihedral_restraint_lists'),
                  ('nef_rdc_restraint_list', 'get_rdc_restraint_lists'),
                  ('nef_nmr_spectrum', 'get_nmr_spectra'),
                  ('nef_peak_restraint_links', 'get_peak_restraint_links')]

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

NEF_STANDARD = 'standard'
NEF_SILENT = 'silent'
NEF_STRICT = 'strict'

NEFVALID = 0
NEFERROR_GENERICGETTABLEERROR = -1
NEFERROR_TABLEDOESNOTEXIST = -2
NEFERROR_SAVEFRAMEDOESNOTEXIST = -3
NEFERROR_ERRORLOADINGFILE = -4
NEFERROR_ERRORSAVINGFILE = -5
NEFERROR_BADTOSTRING = -6
NEFERROR_BADFROMSTRING = -7
NEFERROR_BADMULTICOLUMNVALUES = -8
NEFERROR_BADTABLENAMES = -9
NEFERROR_LISTTYPEERROR = -10
NEFERROR_BADLISTTYPE = -11
NEFERROR_BADFUNCTION = -12
NEFERROR_BADCATEGORIES = -13
NEFERROR_BADADDSAVEFRAME = -14
NEFERROR_READATTRIBUTENAMES = -15
NEFERROR_READATTRIBUTE = -16

class _errorLog():
  """
  A class to facilitate Logging of errors to stderr.

  example:
    errorLogging.logError('Error: %s' % errorMessage)

  functions available are:

    logError              write message to the current output
    func = logger         return the current logger
    logger = func         set the current logger
    loggingMode = mode    set the logging mode where mode is:
                          'standard', 'silent', 'strict'

          current modes are:
              'standard'      errors are written to the stderr, no errors are raised
              'silent'        no errors are raised, no output to the stderr
              'strict'        errors are logged to stderr and errors are raised
                              to be handled by the calling functions

    mode = logginMode     return the current mode.
  """
  _availableModes = (NEF_STANDARD, NEF_SILENT, NEF_STRICT)
  NEFERRORS = {NEFERROR_BADADDSAVEFRAME: 'bad add saveFrame'
              , NEFERROR_BADCATEGORIES: 'bad categories'
              , NEFERROR_BADFUNCTION: ''
              , NEFERROR_BADLISTTYPE: 'bad listType'
              , NEFERROR_BADTABLENAMES: 'bad table names'
              , NEFERROR_LISTTYPEERROR: 'list type error'
              , NEFERROR_BADMULTICOLUMNVALUES: 'bad multiColumnValues'
              , NEFERROR_BADFROMSTRING: 'bad convert from string'
              , NEFERROR_BADTOSTRING: 'bad convert to string'
              , NEFERROR_ERRORSAVINGFILE: 'error saving file'
              , NEFERROR_ERRORLOADINGFILE: 'error loading file'
              , NEFERROR_SAVEFRAMEDOESNOTEXIST: 'saveFrame does not exist'
              , NEFERROR_TABLEDOESNOTEXIST: 'table does not exist'
              , NEFERROR_GENERICGETTABLEERROR: 'table error'
              , NEFERROR_READATTRIBUTENAMES: 'error reading attribute names'
              , NEFERROR_READATTRIBUTE: 'error reading attribute'}

  def __init__(self, logOutput=sys.stderr.write, loggingMode=NEF_STANDARD, errorCode=NEFVALID):
    self._logOutput = logOutput
    self._loggingMode = loggingMode
    self._lastError = errorCode

  def __call__(self, func):
    # method that is called when used as a decorator for functions
    import traceback

    @wraps(func)
    def errortesting(obj, *args, **kwargs):
      try:
        obj._logError(errorCode=NEFVALID)
        return func(obj, *args, **kwargs)
      except Exception as es:
        _type, _value, _traceback = sys.exc_info()
        obj._logError(errorCode=NEFERROR_BADFUNCTION, errorString=str(obj)+str(_type)+str(es))
        return None
    return errortesting

  @property
  def logger(self):
    # return the current logger
    return self._logOutput

  @logger.setter
  def logger(self, func):
    # set the current logger
    self._logOutput = func

  @property
  def loggingMode(self):
    # return the current logger
    return self._loggingMode

  @loggingMode.setter
  def loggingMode(self, loggingMode):
    # set the current logger
    if loggingMode in self._availableModes:
      self._loggingMode = loggingMode

  @property
  def lastErrorString(self):
    # return the error code of the last action
    return self._lastErrorString

  @property
  def lastError(self):
    # return the error code of the last action
    return self._lastError

  def _clearError(self):
    self._lastError = NEFVALID
    self._lastErrorString = ''

  def _logError(self, errorCode=NEFVALID, errorString=''):
    # log errorCode to the current logger
    self._clearError()
    if errorCode != NEFVALID:
      self._lastError = errorCode
      if not errorString:
        errorString = self.NEFERRORS[errorCode]
      self._lastErrorString = 'runtimeError: '+str(errorCode)+' '+errorString+'\n'
      if self._loggingMode != NEF_SILENT:
        try:
          self._logOutput(self._lastErrorString)
        except Exception as es:
          raise es
      if self._loggingMode == NEF_STRICT:
        raise Exception(str(self._lastErrorString))
    pass

# class errorcheck():
#   """
#   Class to wrap functions with a standard error message
#   """
#   def __init__(self, errorCode=NEFERROR_BADFUNCTION):
#     self._errorCode = errorCode
#
#   def __call__(self, func):
#     @wraps(func)
#     def errortesting(obj, *args, **kwargs):
#       try:
#         obj.lastError = NEFVALID
#         return func(obj, *args, **kwargs)
#       except Exception as es:
#         obj.lastError = self._errorCode     #, func, str(es))
#         return None
#     return errortesting


class NefDict(StarIo.NmrSaveFrame):
  # class to add functions to a saveFrame
  def __init__(self, inFrame, _errorLogger=None):
    super(NefDict, self).__init__(name=inFrame.name)
    self._nefFrame = inFrame
    self._errorLogger = _errorLogger

  def _namedToOrderedDict(self, frame):
    # change a saveFrame into a normal OrderedDict
    newItem = OrderedDict()
    for ky in frame.keys():
      newItem[ky] = frame[ky]
    return newItem

  @_errorLog(errorCode=NEFERROR_BADTABLENAMES)
  def getTableNames(self):
    # return table 'name' if exists else None
    return tuple([db for db in self._nefFrame.keys()
              if isinstance(self._nefFrame[db], StarIo.NmrLoop)])

  @_errorLog(errorCode=NEFERROR_GENERICGETTABLEERROR)
  def getTable(self, name=None, asPandas=False):
    # return table 'name' if exists else None
    thisFrame = None
    if name:
      if name in self._nefFrame:
        thisFrame = self._nefFrame[name]
      else:
        # table not found
        self._logError(errorCode=NEFERROR_TABLEDOESNOTEXIST)
        return None
    else:
      tables = self.getTableNames()
      if tables:
        thisFrame = self._nefFrame[tables[0]]

    if asPandas:
      return self._convertToPandas(thisFrame)
    else:
      return [self._namedToOrderedDict(sf) for sf in thisFrame.data]

  @_errorLog(errorCode=NEFERROR_BADMULTICOLUMNVALUES)
  def multiColumnValues(self, column=None):
    return self._nefFrame.multiColumnValues(column=column)

  @_errorLog(errorCode=NEFERROR_TABLEDOESNOTEXIST)
  def hasTable(self, name):
    # return True is the table exists in the saveFrame
    return name in self._nefFrame

  def setTable(self, name):
    # add the table 'name' to the saveFrame, or replace the existing
    # does this need to be here or in the main class?
    print ('Not implemented yet.')

  def _convertToPandas(self, sf):
    try:
      df = pd.DataFrame(data=sf.data, columns=sf.columns)
      df.replace({'.': np.NAN, 'true': True, 'false': False}, inplace=True)
      return df
    except Exception as es:
      return None

  @_errorLog(errorCode=NEFERROR_READATTRIBUTENAMES)
  def getAttributeNames(self):
    return tuple([db for db in self._nefFrame.keys()
              if not isinstance(self._nefFrame[db], StarIo.NmrLoop)])

  @_errorLog(errorCode=NEFERROR_READATTRIBUTE)
  def getAttribute(self, name):
    return self._nefFrame[name]

  @_errorLog(errorCode=NEFERROR_READATTRIBUTE)
  def hasAttribute(self, name):
    return name in self._nefFrame

  @property
  def lastErrorString(self):
    # return the error code of the last action
    return self._errorLogger._lastErrorString

  @property
  def lastError(self):
    # return the error code of the last action
    return self._errorLogger.lastError

  def _logError(self, errorCode=NEFVALID, errorString=''):
    # return the error code of the last action
    self._errorLogger._logError(errorCode=errorCode, errorString=errorString)


class NefImporter():
  """Top level data block for accessing object tree"""
  # put functions in here to read the contents of the dict.
  # superclassed from DataBlock which is of type StarContainer
  def __init__(self, name=None
               , programName='Unknown'
               , programVersion='Unknown'
               , initialise=True
               , errorLogging=NEF_STANDARD):
    self.name = name
    self._nefDict = StarIo.NmrDataBlock()     # empty block
    self.programName = programName
    self.programVersion = programVersion
    self._errorLogger = _errorLog(loggingMode=errorLogging)

    if initialise:
      # initialise a basic object
      self.initialise()

    # short loop to add methods to the class based on the names in the NEF_CATEGORIES list
    # from functools import partial
    # for nefCategory in NEF_CATEGORIES:
    #   setattr(self.__class__, nefCategory[1], partial(self._getListType, nefCategory[0]))

  def _namedToNefDict(self, frame):
    # change a saveFrame into a normal OrderedDict
    newItem = NefDict(inFrame=frame, _errorLogger=self._errorLogger)
    for ky in frame.keys():
      newItem[ky] = frame[ky]
    return newItem

  @_errorLog(errorCode=NEFERROR_BADLISTTYPE)
  def _getListType(self, _listType):
    # return a list of '_listType' from the saveFrame, used with nefCategory
    if self._nefDict and isinstance(self._nefDict, OrderedDict):
      sfList = [self._nefDict[db] for db in self._nefDict.keys() if _listType in db]
      sfList = [self._namedToNefDict(sf) for sf in sfList]

      if len(sfList) > 1:
        return sfList
      elif sfList:
        return sfList[0]

    return None

  # routines to get the Nef specific data from the dictionary
  # this was done using the metafunction above but this doesn't give hints
  def getNmrMetaData(self):
    return self._getListType(NEF_CATEGORIES[0][0])

  def getMolecularSystems(self):
    return self._getListType(NEF_CATEGORIES[1][0])

  def getChemicalShiftLists(self):
    return self._getListType(NEF_CATEGORIES[2][0])

  def getDistanceRestraintLists(self):
    return self._getListType(NEF_CATEGORIES[3][0])

  def getDihedralRestraintLists(self):
    return self._getListType(NEF_CATEGORIES[4][0])

  def getRdcRestraintLists(self):
    return self._getListType(NEF_CATEGORIES[5][0])

  def getNmrSpectra(self):
    return self._getListType(NEF_CATEGORIES[6][0])

  def getPeakRestraintLinks(self):
    return self._getListType(NEF_CATEGORIES[7][0])

  def initialise(self):
    self._nefDict['nef_nmr_meta_data'] = StarIo.NmrDataBlock()
    self._nefDict['nef_nmr_meta_data'].update({k:'' for k in MD_REQUIRED_FIELDS})
    self._nefDict['nef_nmr_meta_data']['sf_category'] = 'nef_nmr_meta_data'
    self._nefDict['nef_nmr_meta_data']['sf_framecode'] = 'nef_nmr_meta_data'
    self._nefDict['nef_nmr_meta_data']['format_name'] = 'Nmr_Exchange_Format'
    self._nefDict['nef_nmr_meta_data']['format_version'] = __nef_version__
    self._nefDict['nef_nmr_meta_data']['program_name'] = self.programName
    self._nefDict['nef_nmr_meta_data']['program_version'] = self.programVersion
    # for l in Nef.MD_REQUIRED_LOOPS:
    #     self['nef_nmr_meta_data'][l] = []

    self._nefDict['nef_molecular_system'] = StarIo.NmrDataBlock()
    self._nefDict['nef_molecular_system'].update({k:'' for k in MS_REQUIRED_FIELDS})
    self._nefDict['nef_molecular_system']['sf_category'] = 'nef_molecular_system'
    self._nefDict['nef_molecular_system']['sf_framecode'] = 'nef_molecular_system'
    for l in MS_REQUIRED_LOOPS:
      self._nefDict['nef_molecular_system'][l] = []

      self.addChemicalShiftList('nef_chemical_shift_list_1', 'ppm')

    self._logError(errorCode=NEFVALID)

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
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

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addChemicalShiftList(self, name, cs_units='ppm'):
    category = 'nef_chemical_shift_list'
    self.addSaveFrame(name=name, category=category,
                       required_fields=CSL_REQUIRED_FIELDS,
                       required_loops=CSL_REQUIRED_LOOPS)
    self._nefDict[name]['atom_chem_shift_units'] = cs_units
    return self._nefDict[name]

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addDistanceRestraintList(self, name, potential_type,
                                  restraint_origin=None):
    category = 'nef_distance_restraint_list'
    self.addSaveFrame(name=name, category=category,
                       required_fields=DRL_REQUIRED_FIELDS,
                       required_loops=DRL_REQUIRED_LOOPS)
    self._nefDict[name]['potential_type'] = potential_type
    if restraint_origin is not None:
      self._nefDict[name]['restraint_origin'] = restraint_origin

    return self._nefDict[name]

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addDihedralRestraintList(self, name, potential_type,
                                  restraint_origin=None):
    category = 'nef_dihedral_restraint_list'
    self.addSaveFrame(name=name, category=category,
                       required_fields=DIHRL_REQUIRED_FIELDS,
                       required_loops=DIHRL_REQUIRED_LOOPS)
    self._nefDict[name]['potential_type'] = potential_type
    if restraint_origin is not None:
      self._nefDict[name]['restraint_origin'] = restraint_origin

    return self._nefDict[name]

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addRdcRestraintList(self, name, potential_type,
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

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addPeakList(self, name, num_dimensions, chemical_shift_list,
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

  @_errorLog(errorCode=NEFERROR_BADADDSAVEFRAME)
  def addLinkageTable(self):
    name = category = 'nef_peak_restraint_links'
    return self.addSaveFrame(name=name, category=category,
                              required_fields=PRLS_REQUIRED_FIELDS,
                              required_loops=PL_REQUIRED_LOOPS)

  @_errorLog(errorCode=NEFERROR_BADTOSTRING)
  def toString(self):
    return self._nefDict.toString()

  # should this be static
  @_errorLog(errorCode=NEFERROR_BADFROMSTRING)
  def fromString(self, text, strict=True):
    # set the Nef from the contents of the string, opposite of toString
    dataExtent = StarIo.parseNef(text=text)
    if dataExtent:
      dbs = [dataExtent[db] for db in dataExtent.keys()]
      if dbs:
        self._nefDict = dbs[0]
    else:
      self.lastError = NEFERROR_BADFROMSTRING
      self._nefDict = None

  @_errorLog(errorCode=NEFERROR_ERRORLOADINGFILE)
  def loadFile(self, fileName=None, mode='standard'):
    nefDataExtent = StarIo.parseNefFile(fileName=fileName, mode=mode)
    self._nefDict = list(nefDataExtent.values())
    if len(self._nefDict) > 1:
      print(
        'More than one datablock in a NEF file is not allowed.  Using the first and discarding the rest.')
    self._nefDict = self._nefDict[0]

    return True

  @_errorLog(errorCode=NEFERROR_ERRORSAVINGFILE)
  def saveFile(self, fileName=None, mode='standard'):
    with open(fileName, 'w') as op:
      op.write(self._nefDict.toString())

    return True

  @_errorLog(errorCode=NEFERROR_BADCATEGORIES)
  def getCategories(self):
    # return a list of the categories available in a Nef file
    return tuple([nm[0] for nm in NEF_CATEGORIES])

  @_errorLog(errorCode=NEFERROR_SAVEFRAMEDOESNOTEXIST)
  def getSaveFrameNames(self, filter=NEF_RETURNALL):
    # return a list of the saveFrames in the file
    names = [db for db in self._nefDict.keys()
             if isinstance(self._nefDict[db], StarIo.NmrSaveFrame)]

    if filter == NEF_RETURNNEF:
      names = [nm for nm in names if nm and nm.startswith(NEF_PREFIX)]
    elif filter == NEF_RETURNOTHER:
      names = [nm for nm in names if nm and not nm.startswith(NEF_PREFIX)]

    return tuple(names)

  @_errorLog(errorCode=NEFERROR_SAVEFRAMEDOESNOTEXIST)
  def hasSaveFrame(self, name):
    # return True if the saveFrame exists, else False
    return name in self._nefDict

  @_errorLog(errorCode=NEFERROR_SAVEFRAMEDOESNOTEXIST)
  def getSaveFrame(self, sfName):
    # return the saveFrame 'name'
    return NefDict(self._nefDict[sfName], _errorLogger=self._errorLogger)

  @_errorLog(errorCode=NEFERROR_READATTRIBUTENAMES)
  def getAttributeNames(self):
    names = [db for db in self._nefDict.keys()
             if not isinstance(self._nefDict[db], StarIo.NmrSaveFrame)]

    return tuple(names)

  @_errorLog(errorCode=NEFERROR_READATTRIBUTE)
  def getAttribute(self, name):
    return self._nefDict[name]

  @_errorLog(errorCode=NEFERROR_READATTRIBUTE)
  def hasAttribute(self, name):
    return name in self._nefDict

  @property
  def lastErrorString(self):
    # return the error code of the last action
    return self._errorLogger._lastErrorString

  @property
  def lastError(self):
    # return the error code of the last action
    return self._errorLogger.lastError

  def _logError(self, errorCode=NEFVALID, errorString=''):
    # return the error code of the last action
    self._errorLogger._logError(errorCode=errorCode, errorString=errorString)


if __name__ == '__main__':
  from ccpn.ui.gui.widgets.Application import TestApplication

  app = TestApplication()

  # TODO:ED write a test suite for this
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

  sf1 = None
  if names:
    sf1 = test.getSaveFrame(names[0])
  sf2 = test.getSaveFrame('notFound')

  print (test.hasSaveFrame('notFound'))
  if names:
    print (test.hasSaveFrame(names[0]))

  ts = test.toString()
  test.fromString(ts)

  if sf1 is not None:
    print (sf1.name)
    print (sf1.getTableNames())

    table = sf1.getTable()
    sf1.getTable('nmr_atom', asPandas=True)
    table = sf1.getTable('nmr_atom', asPandas=True)

    print (table)

    print (sf1.hasTable('notFound'))
    print (sf1.hasTable('nmr_residue'))

    print (sf1.getAttributeNames())
    print (sf1.getAttribute('notHere'))
    print (sf1.getAttribute('sf_category'))
    print (sf1.hasAttribute('sf_framecode'))
    print (sf1.hasAttribute('nothing'))

  print ('Testing getTable')
  try:
    print (test.getSaveFrame(sfName='ccpn_assignment').getTable())
  except Exception as es:
    print ('Error: %s' % str(es))

  try:
    print (test.getSaveFrame(sfName='ccpn_assignment').getTable(name='nmr_residue', asPandas=True))
  except Exception as es:
    print ('Error: %s' % str(es))

  try:
    print (test.getSaveFrame(sfName='ccpn_assignment').getTable(name='notFound', asPandas=True))
  except Exception as es:
    print ('Error: %s' % str(es))

  print ('Testing saveFile')
  print ('SAVE ', test.saveFile('/Users/ejb66/PycharmProjects/Sec5Part3testing.nef'))
  print (test.lastError)

  # test meta creation of category names
  print (test.getNmrMetaData())
  print (test.getMolecularSystems())
  print (test.getDihedralRestraintLists())
  print (test.getDistanceRestraintLists())
  print (test.getChemicalShiftLists())
  print (test.getNmrSpectra())
  print (test.getPeakRestraintLinks())
  print (test.getRdcRestraintLists())

  # sp = test.getNmrSpectra()
  # if sp:
  #   for spp in sp:
  #     print (spp.getTableNames())
  #     print (spp.getTable('nef_spectrum_dimension', asPandas=True))
  #
  # print (test.getMolecularSystems().getTable('nef_sequence'))
  # print (test.getChemicalShiftLists().getTable(asPandas=True))
  #
  # print (test._nefDict.tagPrefix)

  print (test.getAttributeNames())