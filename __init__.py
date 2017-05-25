"""Star, Cif, and NEF I/O.

Files are converted to and from a tree of nested objects (essentially enhanced, ordered dictionaries)

The GenericStarParser.py will work for any valid Star file; see this file
for the precise reading and writing behaviour. Parsing comes in several modes, depending on how
strictly the STAR standard is enforced.

NmrStar and NEF files have a more restricted syntax, with all items and loops contained within
saveframes, and with saveframe and loop tags beginning with a prefix matching the saveframe and loop
category. The StarIo.py code is used for these files; it converts the output of the GenericStarParser
into a simpler object tree that makes use of these restrictions, strips mandatory tag prefixes,
and converts strings to numerical values with some limited heuristics.

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
__dateModified__ = "$dateModified: 2017-04-07 11:40:46 +0100 (Fri, April 07, 2017) $"
__version__ = "$Revision: 3.0.b1 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: CCPN $"

__date__ = "$Date: 2017-04-07 10:28:41 +0000 (Fri, April 07, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================



