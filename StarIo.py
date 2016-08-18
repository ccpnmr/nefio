"""I/O for NEF and NmrStar formats

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date$"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon Skinner, Geerten Vuister"
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

# NB Assumes that file was parsed with lowercaseTags = True

# Constraints on class hierarchy:
# - all plain tags in a saveframe start with the '<sf_category>.'
# - ALl loop column names start with '<loopcategory>.'
# - loopcategories share a namespace with tags within a saveframe
# - DataBlocks can contain only saveframes.
# - tag prefixes ('_', 'save_', 'data_' are stripped


# TODO NBNB some longer variable names;

import keyword
import os
import typing

from . import GenericStarParser

# # Translations from internal values to framecode-compatible strings
# replacementCharacter = '`'
# _translationDict = {
#   '#':'?', '"':'?', "'":'?'
# }
# for indx in range(128):
#   char = chr(indx)
#   if char.isspace:
#     # replace space characters with '_'
#     _translationDict[char] = '_'
#   elif not char.isprintable():
#     _translationDict[char] = '?'
# for index in range(128,256):
#   _translationDict[chr(indx)] =
# strToFramecode = str.maketrans(_translationDict)

# Make to-=ascii mapping at byte level'
# Space and unprintable characters map to '_'
# Characters outside the ascii range and #'" map to '?'
#
# NBNB this could be done better (more readable) with custom mappings
#
_replaceSpaceWithByte = ord('_')
_replaceWithByte = ord('?')
_b2 = bytearray(range(256))
for ii in range(128):
  _char = chr(ii)
  if _char.isspace():
    _b2[ii] = _replaceSpaceWithByte
  elif not _char.isprintable():
    _b2[ii] = _replaceSpaceWithByte
for _char in ('#', '"', "'"):
  _b2[ord(_char)] = _replaceWithByte
for ii in range(128,256):
  _b2[ii] = _replaceWithByte
latin_1_to_framecode_translator = bytes.maketrans(bytearray(range(256)), _b2)


def string2FramecodeString(text):

  # Get 8-bit bytes, replacing values >255 with ord('?')
  bts = text.encode('latin_1', 'replace')

  # Translate bytes, using preset translator
  bts2 = bts.translate(latin_1_to_framecode_translator)

  # convert back to string - all bytes are ascii, so no errors can occur
  return str(bts2.decode('ascii'))


def parseNmrStar(text:str, mode='standard'):
  """load NMRSTAR file"""
  dataExtent = GenericStarParser.parse(text, mode)
  converter = _StarDataConverter(dataExtent)
  converter.preValidate()
  result = converter.convert()
  #
  return result

def parseNef(text:str, mode='standard'):
  """load NEF from string"""

  dataExtent = GenericStarParser.parse(text, mode)
  converter = _StarDataConverter(dataExtent, fileType='nef')
  converter.preValidate()
  result = converter.convert()
  #
  return result

def parseNmrStarFile(fileName:str, mode:str='standard', wrapInDataBlock:bool=False):
  """parse NMRSTAR from file"""
  text = open(fileName).read()
  if wrapInDataBlock and 'save_' in text and not 'data_' in text:
    text = "data_dummy \n\n" + text
  dataExtent = GenericStarParser.parse(text, mode)
  converter = _StarDataConverter(dataExtent)
  converter.preValidate()
  result = converter.convert()
  #
  return result

def parseNefFile(fileName:str, mode:str='standard', wrapInDataBlock:bool=False):
  """parse NEF from file"""
  text = open(fileName).read()
  if wrapInDataBlock and 'save_' in text and not 'data_' in text:
    text = "data_dummy \n\n" + text
  dataExtent = GenericStarParser.parse(text, mode)
  converter = _StarDataConverter(dataExtent, fileType='nef')
  converter.preValidate()
  result = converter.convert()
  #
  return result

class StarValidationError(ValueError):
  pass

class NmrDataExtent(GenericStarParser.DataExtent):
  """Top level container (OrderedDict) for NMRSTAR/NEF object tree"""
  pass

# We insert these afterwards as we want the functions at the top of the file
# but can only annotate after DataExtent is created
parseNef.__annotations__['return'] = NmrDataExtent
parseNefFile.__annotations__['return'] = NmrDataExtent
parseNmrStar.__annotations__['return'] = NmrDataExtent
parseNmrStarFile.__annotations__['return'] = NmrDataExtent

class NmrLoop(GenericStarParser.Loop):
  """Loop for NMRSTAR/NEF object tree

  The contents, self.data is a list of OrderedDicts matching the column names.
  rows can be modified or deleted from data, but adding new rows directly is likely to
  break - use the newRow function."""

  @property
  def category(self) -> str:
    """Loop category tag - synonym for name (unlike the case of SaveFrame)"""
    return self.name

  @property
  def tagPrefix(self) -> str:
    """Prefix to use before item tags on output"""
    return '_%s.' % self.name

class NmrSaveFrame(GenericStarParser.SaveFrame):
  """SaveFrame (OrderedDict)for NMRSTAR/NEF object tree"""

  def __init__(self, name=None, category=None):
    super(NmrSaveFrame, self).__init__(name=name)
    self.category = category

  @property
  def tagPrefix(self) -> str:
    """Prefix to use before item tags on output"""
    return '_%s.' % self.category

  def newLoop(self, name, columns) -> NmrLoop:
    """Make new NmrLoop and add it to the NmrSaveFrame"""
    loop = NmrLoop(name, columns)
    self.addItem(name, loop)
    return loop

class NmrDataBlock(GenericStarParser.DataBlock):
  """DataBlock (OrderedDict)for NMRSTAR/NEF object tree"""

  def newSaveFrame(self, name:str, category:str) -> NmrSaveFrame:
    """Make new NmrSaveFrame and add it to the DataBlock"""
    name = string2FramecodeString(name)
    saveFrame = NmrSaveFrame(name, category=category)
    self.addItem(name, saveFrame)
    saveFrame.addItem('sf_category', category)
    saveFrame.addItem('sf_framecode', name)
    return saveFrame

  def addSaveFrame(self, saveFrame:NmrSaveFrame):
    """Add existing NmrSaveFrame to the DataBlock"""
    self.addItem(saveFrame['sf_framecode'], saveFrame)

class NmrLoopRow(GenericStarParser.LoopRow):
  pass


class _StarDataConverter:
  """Converter from output of a GeneralStarParser to a NEF or NMRSTAR nested data structure

  NB Function assumes valid data as output from GeneralStarParser with lowerCaseTags settings
  and does not double check validity."""


  validFileTypes = ('nef', 'star')

  def __init__(self, dataExtent: GenericStarParser.DataExtent, fileType='star',
               specification=None, convertColumnNames=True):

    # Set option settings
    if specification is None:
      self.specification = None
    else:
      raise NotImplementedError("_StarDataConverter specification input not yet implemented")
    fileType = fileType and fileType.lower()
    if fileType not in self.validFileTypes:
      raise StarValidationError("fileType %s must be one of %s" % (fileType, self.validFileTypes))
    self.fileType = fileType

    self.convertColumnNames = convertColumnNames

    self.dataExtent = dataExtent

    # Stack of objects parsed, to give context for error messages
    self.stack = []

  def preValidate(self):
    self.stack = []

    try:
      for dataBlock in self.dataExtent.values():
        self.preValidateDataBlock(dataBlock)
    except StarValidationError:
      raise
    except:
      print(self._errorMessage('System error:'))
      raise

  def convert(self):

    nmrDataExtent = NmrDataExtent(name=self.dataExtent.name)

    self.stack = []

    try:
      for dataBlock in self.dataExtent.values():
        newDataBlock =  self.convertDataBlock(dataBlock)
        nmrDataExtent.addItem(newDataBlock.name, newDataBlock)
    except StarValidationError:
      raise
    except:
      print(self._errorMessage('System error:'))
      raise
    #
    return nmrDataExtent

  def preValidateDataBlock(self, dataBlock):

    self.stack.append(dataBlock)

    name = dataBlock.name
    if name != 'global_' and not name.startswith('data_'):
      self.raiseValidationError("DataBlock name  must be 'global_' or start with 'data_'")

    for tag, saveFrame in dataBlock.items():

      if isinstance(saveFrame, GenericStarParser.SaveFrame):
        self.preValidateSaveFrame(saveFrame)
      else:
        self.raiseValidationError("%s file DataBlock contains non-saveframe element %s:%s"
                         % (self.fileType, tag, saveFrame))

    self.stack.pop()

  def convertDataBlock(self, dataBlock):

    self.stack.append(dataBlock)

    # get NmrDataBlock name
    name = dataBlock.name
    if name.startswith('data_'):
      name = name[5:] or '__MissingDataBlockName'
    elif name == 'global_':
      name = 'global'

    # Make NmrDataBlock and connect it
    nmrDataBlock = NmrDataBlock(name=name)

    for saveFrame in dataBlock.values():
      nmrSaveFrame = self.convertSaveFrame(saveFrame)
      nmrDataBlock.addItem(nmrSaveFrame.name, nmrSaveFrame)
    #
    self.stack.pop()
    return nmrDataBlock


  def preValidateSaveFrame(self, saveFrame):

    self.stack.append(saveFrame)
    commonPrefix = os.path.commonprefix([tt[0] for tt in saveFrame.items()
                                         if isinstance(tt[1], str)])
    tt = commonPrefix.split('.', 1)
    if len(tt) == 2:
      prefix = tt[0] + '.'
    else:
      self.raiseValidationError(
        "Saveframe tags do not start with a common dot-separated prefix: %s"
        % [tt[0] for tt in saveFrame.items() if isinstance(tt[1], str)]
      )

    sf_category = saveFrame.get(prefix + 'sf_category')
    if sf_category is None:
      self.raiseValidationError("SaveFrame lacks .sf_category item")
    sf_framecode = saveFrame.get( prefix + 'sf_framecode')
    if sf_framecode is None:
      self.raiseValidationError("SaveFrame lacks .sf_framecode item")

    sf_lowername = saveFrame.name   # NB tags are lower-cased from the parser
    if sf_lowername.startswith('save_'):
      sf_lowername = sf_lowername[5:]
    if sf_lowername != sf_framecode.lower():
      self.raiseValidationError("Saveframe.name %s does not match sf_framecode %s"
                       % (sf_lowername, sf_framecode))

    if self.fileType == 'nef':
      if not sf_framecode.startswith(sf_category):
        self.raiseValidationError("NEF file sf_framecode %s does not start with the sf_category %s" %
                         (sf_framecode, sf_category))
      if prefix[1:-1] != sf_category:
        self.raiseValidationError("NEF file sf_category %s does not match tag prefix %s" %
                         (sf_category, prefix))
    else:
      # NBNB TBD We do not check or store the tag prefix
      pass

    for tag, value in saveFrame.items():
      self.stack.append(tag)
      if isinstance(value, GenericStarParser.Loop):
        if tag == value.name:
          self.preValidateLoop(value)
      elif not isinstance(value, str):
        self.raiseValidationError("Saveframe contains item value of wrong type: %s"
                                  % value)
      self.stack.pop()

    self.stack.pop()

  def convertSaveFrame(self, saveFrame):

    self.stack.append(saveFrame)

    #Get common dot-separated prefix from non-loop items
    commonPrefix = os.path.commonprefix([tt[0] for tt in saveFrame.items()
                                         if isinstance(tt[1], str)])
    tt = commonPrefix.split('.', 1)
    if len(tt) == 2:
      prefix = tt[0] + '.'
    else:
      self.raiseValidationError(
        "Saveframe tags do not start with a common dot-separated prefix: %s"
        % [tt[0] for tt in saveFrame.items() if isinstance(tt[1], str)]
      )

    # get category and framecode
    # The prevalidation has already established that there is exactly one tag for each
    tags = [x for x in saveFrame if x.endswith('.sf_framecode')]
    sf_framecode = saveFrame[tags[0]]
    tags = [x for x in saveFrame if x.endswith('.sf_category')]
    sf_category = saveFrame[tags[0]]

    newSaveFrame = NmrSaveFrame(name=sf_framecode, category=sf_category)

    lowerCaseCategory = newSaveFrame.category.lower()
    for tag, value in saveFrame.items():

      self.stack.append(tag)

      if isinstance(value, str):
        value = self.convertValue(value, category=lowerCaseCategory, tag=tag)
        objname = tag[len(prefix):]
        newSaveFrame.addItem(objname, value)

      elif isinstance(value, GenericStarParser.Loop):
        if tag == value.columns[0]:
          # Only add loop on first appearance
          nmrLoop = self.convertLoop(value)
          newSaveFrame.addItem(nmrLoop.name, nmrLoop)
      self.stack.pop()
    #
    self.stack.pop()
    return newSaveFrame

  def preValidateLoop(self, loop):

    self.stack.append(loop)

    columns = loop.columns
    commonPrefix = os.path.commonprefix(columns)
    if len(commonPrefix.split('.', 1)) != 2:
      self.raiseValidationError(
        "Column names of %s do not start with a common dot-separated prefix: %s" % (loop, columns)
      )

    self.stack.pop()

  def convertLoop(self, loop):

    self.stack.append(loop)

    oldColumns = loop.columns
    commonPrefix = os.path.commonprefix(oldColumns)
    tt = commonPrefix.split('.', 1)
    if len(tt) == 2:
      category = tt[0]
      lenPrefix = len(category) + 1
      if category[0] == '_':
        category = category[1:]
    else:
      self.raiseValidationError(
        "Column names of %s do not start with a common dot-separated prefix: %s" % (loop,
                                                                                    oldColumns)
      )

    columns = []
    for ss in oldColumns:
      tag = ss[lenPrefix:]

      # Check for valid field names
      if tag.startswith('_') or not tag.isidentifier():
        if self.convertColumnNames:
          # make the name valid
          tag = ''.join(x if x.isalnum() else '_' for x in tag)
          while tag.startswith('_'):
            tag = tag[1:]
        else:
          raise ValueError("Invalid column name: %s" % ss)

      if keyword.iskeyword(tag):
        raise ValueError("column name (as modified) clashes with Python keyword: %s" % ss)

      columns.append(tag)

    newLoop = NmrLoop(category, columns)
    ff = self.convertValue #convertValue(value, category=lowerCaseCategory, tag=tag)
    for row in loop.data:
      values = [ff(x, category, columns[ii]) for ii,x in enumerate(row.values()) ]
      newLoop.newRow(values)

    #
    self.stack.pop()
    return newLoop

  def convertValue(self, value, category=None, tag=None):
    """NB category and tag are needed if  we want to use self.specification"""

    if self.specification:
      # Add specification-dependent processing here
      #
      return value

    if isinstance(value, GenericStarParser.UnquotedValue):
      # Convert special values
      if value == GenericStarParser.NULLSTRING :
        # null  value
        value = None
      elif value == GenericStarParser.TRUESTRING:
        # Boolean True
        value = True
      elif value == GenericStarParser.FALSESTRING:
        # Boolean False
        value = False
      elif value == GenericStarParser.UNKNOWNSTRING:
        value = None
      elif value[0] == '$':
        # SaveFrame reference
        value = value[1:]
      else:
        if not tag[-5:] in ('_code', '_name'):
          # HACK - tags ending in '_code' or '_name' are assumed to be string type
          # This takes care of e.g. 'sequence_code'
          # that often might evaluate to a number otherwise
          try:
            value = int(value)
          except ValueError:
            try:
              value = float(value)
            except ValueError:
              pass
    #
    return value

  def _errorMessage(self, msg):
    """Make standard error message"""
    template = "Error in context: %s\n%s"
    ll = [(x if isinstance(x, str) else x.name) for x in self.stack]
    return template % (ll, msg)

  def raiseValidationError(self, msg):
    raise StarValidationError(self._errorMessage(msg))


def splitNefSequence(rows:typing.Sequence[dict]) -> typing.List[typing.List[dict]]:
  """Split a sequence of nef_sequence dicts assumed to belong to the same chain
  into a list of lists of sequentially linked stretches following the NEF rules

  Note that missing linkings are treated as 'middle' and missing start/end tags
  are ignored, with the first/last residue treated, effectively, as linking 'break'

  Only unknown linking values and incorrect pairs of 'cyclic' tags raise an error"""

  result = []
  stretch = []
  inCyclic = False
  for row in rows:
    linking = row.get('linking')

    if inCyclic and linking not in ('middle', 'cyclic', None):
      raise ValueError(
          "Sequence contains 'cyclic' residue(s) that do not form a closed, cyclic molecule"
        )

    if linking == 'cyclic':
      if inCyclic:
        # End of cycle
        inCyclic = False
        stretch.append(row)
        result.append(stretch)
        stretch = []
      else:
        #start of cycle
        inCyclic = True
        if stretch:
          result.append(stretch)
        stretch = [row]

    elif linking in ('single', 'nonlinear', 'dummy'):
      # Always isolated. And last stretch, add new one, and prepare for the next one
      if stretch:
        result.append(stretch)
      result.append([row])
      stretch = []

    elif linking == 'start':
      # Start new stretch
      if stretch:
        result.append(stretch)
      stretch = [row]

    elif linking == 'end':
      # End stretch
      stretch.append(row)
      result.append(stretch)
      stretch = []

    elif linking in ('middle', None):
      # Continuation (we treat None as 'middle' as the most pragmatic approach
      # Validation of the NEF standard must be done elsewhere
      stretch.append(row)

    elif linking == 'break':

      # TODO NBNB This follows NEF spec as of July 2016 - which is rather confused
      # Propose change so that 'break' signals a chain break AFTER that residue,
      # and is ONLY used if there is a break between two 'middle' residues.

      if stretch:
        # inside stretch - end it
        stretch.append(row)
        result.append(stretch)
        stretch = []
      else:
        # Start of stretch - put row on
        stretch.append(row)


    else:
      raise ValueError("Illegal value of nef_sequence.linking: %s" % linking)

  if stretch:
    # Add final stretch if still open
    result.append(stretch)

  if inCyclic:
      raise ValueError(
          "Sequence contains 'cyclic' residue that is not terminated by matching 'cyclic residue"
        )
  #
  return result
