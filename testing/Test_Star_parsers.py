"""Module Documentation here

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2017"
__credits__ = ("Wayne Boucher, Ed Brooksbank, Rasmus H Fogh, Luca Mureddu, Timothy J Ragan"
               "Simon P Skinner & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpnLicense for licence text")
__reference__ = ("For publications, please use reference from http://www.ccpn.ac.uk/v3-software/downloads/license"
               "or ccpnmodel.ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2017-04-07 11:40:47 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"

__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import os
import time

from .. import GenericStarParser, StarIo
from ccpn.util import Path

TEST_FILE_PATH = os.path.join(Path.getTopDirectory(), 'internal', 'data', 'starExamples')

def _loadGeneralFile(path:str):
  usePath = path if path.startswith('/') else os.path.join(TEST_FILE_PATH, path)
  t0 = time.time()
  entry = GenericStarParser.parseFile(usePath) # 'lenient')
  print ("Parsing time %s for %s" % (time.time() - t0, path))
  return entry

def _loadNmrStarFile(path:str):
  usePath = path if path.startswith('/') else os.path.join(TEST_FILE_PATH, path)
  t0 = time.time()
  entry = StarIo.parseNmrStarFile(usePath) # 'lenient')
  print ("Parsing time %s for %s" % (time.time() - t0, path))
  return entry

def _loadNefFile(path:str):
  usePath = path if path.startswith('/') else os.path.join(TEST_FILE_PATH, path)
  t0 = time.time()
  entry = StarIo.parseNefFile(usePath) # 'lenient')
  print ("Parsing time %s for %s" % (time.time() - t0, path))
  return entry

def _printContents(object, indent=0):
  indentStep = 4
  indent += indentStep
  print (' '*indent, object, 'Tags: %s' %  len(object))
  for tag,obj in object.items():
    if isinstance(obj, GenericStarParser.Loop):
      if tag == obj.columns[0]:
        print (' '*(indent + indentStep), obj, 'Columns: %s' %  len(obj.columns))
    elif not isinstance(obj, str):
      _printContents(obj, indent)

def test_nmrstar_4267():
  print('\n\n', '# nmrstar_4267', '#'*60, '\n')
  _loadGeneralFile('4267_example.str')
  _loadNmrStarFile('4267_example.str')

def test_nef_commented_example():
  _loadGeneralFile('Commented_Example.nef')
  _loadNefFile('Commented_Example.nef')

def test_nef_2l9r_Paris_155():
  print('\n\n', '# Paris_155_nef', '#'*60, '\n')
  _loadGeneralFile('CCPN_2l9r_Paris_155.nef')
  _loadNefFile('CCPN_2l9r_Paris_155.nef')

def test_nef_1lci_Piscataway_179():
  _loadGeneralFile('CCPN_2lci_Piscataway_179.nef')
  _loadNefFile('CCPN_2lci_Piscataway_179.nef')

def test_nef_H1GI():
  _loadGeneralFile('CCPN_H1GI.nef')
  _loadNefFile('CCPN_H1GI.nef')

def test_mmcif_1bgl_1bgm():
  _loadGeneralFile('1bgl_1bgm.cif')

def test_dic_mmcif_nef():
  print('\n\n', '# nef_dic', '#'*60, '\n')
  # _printContents(_loadGeneralFile('mmcif_nef.dic'))
  _loadGeneralFile('mmcif_nef.dic')

def test_dic_mmcif_nmr_star():
  _loadGeneralFile('mmcif_nmr-star.dic')

def test_dic_mmcif_std():
  print('\n\n', '# mmcif_dic', '#'*60, '\n')
  # _printContents(_loadGeneralFile('mmcif_std.dic'))
  _loadGeneralFile('mmcif_std.dic')

def test_dic_mmcif_pdbx_v40():
  _loadGeneralFile('mmcif_pdbx_v40.dic')
