"""
Module Documentation here
"""

from __future__ import unicode_literals, print_function, absolute_import, division


#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (http://www.ccpn.ac.uk) 2014 - 2019"
__credits__ = ("Ed Brooksbank, Luca Mureddu, Timothy J Ragan & Geerten W Vuister")
__licence__ = ("CCPN licence. See http://www.ccpn.ac.uk/v3-software/downloads/license")
__reference__ = ("Skinner, S.P., Fogh, R.H., Boucher, W., Ragan, T.J., Mureddu, L.G., & Vuister, G.W.",
                 "CcpNmr AnalysisAssign: a flexible platform for integrated NMR analysis",
                 "J.Biomol.Nmr (2016), 66, 111-124, http://doi.org/10.1007/s10858-016-0060-y")
#=========================================================================================
# Last code modification
#=========================================================================================
__modifiedBy__ = "$modifiedBy: Ed Brooksbank $"
__dateModified__ = "$dateModified: 2019-12-16 10:40:41 +0000 (Mon, December 16, 2019) $"
__version__ = "$Revision: 3.0.0 $"
#=========================================================================================
# Created
#=========================================================================================
__author__ = "$Author: TJ Ragan $"
__date__ = "$Date: 2017-03-23 16:50:22 +0000 (Thu, March 23, 2017) $"
#=========================================================================================
# Start of code
#=========================================================================================

import re

NMR_EXCHANGE_FORMAT = 'nmr_exchange_format'
FRAME_PREFIX = 'nef_saveframe_'
FRAME_SEARCH = r'{}(\w+)'.format(FRAME_PREFIX)
SF_CATEGORY = 'sf_category'
SF_FRAMECODE = 'sf_framecode'
NAME = 'name'
NEF_ITEM = 'nef_item'
IS_MANDATORY = 'is_mandatory'
IS_KEY = 'is_key'
LOOP_CATEGORY = 'loop_category'
CATEGORY = 'category'
SAVEFRAME = 'saveframe'
NEF_LOOP = 'nef_loop'
METADATA_KEY = 'METADATA'
METADATA_DICTKEY = 'nef_nmr_meta_data'
SPECIFICATION_KEY = 'nef_specification'
FORMAT_NAME = 'format_name'
FORMAT_VERSION = 'format_version'
CREATION_DATE = 'creation_date'
UUID = 'uuid'
VERSION = 'version'
CCPN_PREFIX = 'ccpn_'


class Validator(object):

    def __init__(self, nef=None, validateNefDict=None):
        self.nef = nef
        self.validateNefDict = validateNefDict
        self._validation_errors = None

    def isValid(self, nef=None, validNef=None):
        """Return whether the Nef file is valid
        """
        if nef is None:
            nef = self.nef
        if validNef is None:
            validNef = self.validateNefDict

        self._validateAll(nef, validNef)

        v = list(self._validation_errors.values())
        return not any(v)

    @property
    def validationErrors(self):
        """Return the dict of validation errors
        """
        if not self._validation_errors:
            self.isValid(self.nef, self.validateNefDict)
        return self._validation_errors

    def _validateAll(self, nef=None, validNef=None):
        """Validate a nef file (nef) against a nef dictionary (validNef)
        """
        if nef is None:
            nef = self.nef
        if validNef is None:
            validNef = self.validateNefDict
        if not nef:
            raise RuntimeError('Error: nef not defined')
        if not validNef:
            raise RuntimeError('Error: validateNefDict not defined')

        self._validation_errors = dict()
        self._validation_errors[SAVEFRAME] = []

        # validate meta_data
        e = self._validation_errors[SAVEFRAME]
        e += self._validate_nmr_meta_data(nef, validNef)

        # go through all the saveframes in the Nef object
        for sf_name, saveframe in nef.items():

            if SF_FRAMECODE not in saveframe or saveframe.name != saveframe[SF_FRAMECODE]:
                e = self._validation_errors[SAVEFRAME]
                e += ["Saveframe.name for sf_framecode '{}' is not defined correctly.".format(saveframe[SF_FRAMECODE]), ]
                break

            # check against the validation dictionary
            for vName, validFrame in validNef.items():

                # get the actual name from the end the the name - may need to be more complex later
                checkName = re.findall(FRAME_SEARCH, validFrame.name)

                if checkName and SF_CATEGORY in saveframe and saveframe[SF_CATEGORY] == checkName[0]:

                    ERROR_KEY = checkName[0]
                    e = self._validation_errors[ERROR_KEY] = []

                    # check items against loop_category = None, i.e., this saveframe
                    mandatoryFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is True and nm[LOOP_CATEGORY] is None]
                    optionalFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is False and nm[LOOP_CATEGORY] is None]
                    loopNames = [nm['category'] for nm in validFrame[NEF_LOOP].data]

                    # check for missing words/framecode is not correct/category is mismatched/bad fields (keys)
                    e += self.__sf_framecode_name_mismatch(saveframe, sf_name)
                    e += self.__dict_missing_keys(saveframe, mandatoryFields, label=sf_name)
                    e += self.__sf_category_name_mismatch(saveframe, checkName[0])
                    e += self.__dict_nonallowed_keys(saveframe, mandatoryFields + optionalFields + loopNames, label=sf_name)

                    # iterate through loops
                    loops = [kk for kk in saveframe.keys() if kk in loopNames]
                    for loop in loops:

                        # get the keys that belong to this loop
                        mandatoryLoopFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is True and nm[LOOP_CATEGORY] == loop]
                        optionalLoopFields = [nm[NAME] for nm in validFrame[NEF_ITEM].data if nm[IS_MANDATORY] is False and nm[LOOP_CATEGORY] == loop]

                        if saveframe[loop] and saveframe[loop].data:
                            # check for missing words/bad fields (keys)/malformed loops
                            e += self.__dict_missing_keys(saveframe[loop].data[0], mandatoryLoopFields, label='{}:{}'.format(sf_name, loop))
                            e += self.__dict_nonallowed_keys(saveframe[loop].data[0], mandatoryLoopFields + optionalLoopFields, label='{}:{}'.format(sf_name, loop))
                            e += self.__loop_entries_inconsistent_keys(saveframe[loop].data, label='{}:{}'.format(sf_name, loop))
                        else:

                            # this error is a catch-all as loadFile should test the integrity of the nef file before validation
                            e += ["Error reading loop '{}'.".format(loop), ]

                    break

                elif SF_CATEGORY in saveframe and saveframe[SF_CATEGORY].startswith(CCPN_PREFIX):

                    # skip ccpn_ specific categories for now.
                    break

            else:
                e = self._validation_errors[SAVEFRAME]
                e += ["No sf_category '{}' found (possibly bad name defined).".format(saveframe[SF_CATEGORY]),]

        return self._validation_errors

    def _validate_nmr_meta_data(self, nef, validNef):
        """Check that the information in the meta_data is correct for this version
        """
        DICT_KEY = METADATA_DICTKEY
        VALID_KEY = SPECIFICATION_KEY

        if DICT_KEY not in nef:
            return ['No {} saveframe.'.format(DICT_KEY)]
        else:

            # TODO:ED format name should be defined in the mmcif_nef.dic
            format_name = NMR_EXCHANGE_FORMAT

            if VALID_KEY not in validNef and VERSION not in validNef[VALID_KEY]:
                return ['Version not found']
            format_version = float(validNef[VALID_KEY][VERSION])

            md = nef[DICT_KEY]
            e = []

            if FORMAT_NAME in md:
                if md[FORMAT_NAME] != format_name:
                    e.append("{} must be '{}'.".format(format_name, NMR_EXCHANGE_FORMAT))

            if FORMAT_VERSION in md:
                mdVersion = float(md[FORMAT_VERSION])
                if mdVersion < format_version:
                    e.append('This reader (version {}) does not support {}.'.format(format_version, mdVersion))

            if CREATION_DATE in md:
                # TODO: ED How to validate the creation date?
                pass

            if UUID in md:
                # TODO: ED How to validate the uuid?
                pass

            return e

    # def isValid(self, nef=None):
    #     if nef is None:
    #         nef = self.nef
    #     self.validation_errors = dict()
    #
    #     # self.validation_errors.update(self._validate_datablock(nef))
    #     self.validation_errors.update(self._validate_saveframe_fields(nef))
    #     self.validation_errors.update(self._validate_required_saveframes(nef))
    #     self.validation_errors.update(self._validate_metadata(nef))
    #     self.validation_errors.update(self._validate_molecular_system(nef))
    #     self.validation_errors.update(self._validate_chemical_shift_lists(nef))
    #     self.validation_errors.update(self._validate_distance_restraint_lists(nef))
    #     self.validation_errors.update(self._validate_dihedral_restraint_lists(nef))
    #     self.validation_errors.update(self._validate_rdc_restraint_lists(nef))
    #     self.validation_errors.update(self._validate_peak_lists(nef))
    #     self.validation_errors.update(self._validate_linkage_table(nef))
    #
    #     v = list(self.validation_errors.values())
    #     return not any(v)

    # def _validate_datablock(self, nef=None):
    #     # not required as the validator is subclassed from NefImporter
    #     if nef is None:
    #         nef = self.nef
    #
    #     if not hasattr(nef, 'datablock'):
    #         return {'DATABLOCK': 'No data block specified'}
    #     return {'DATABLOCK': None}
    #
    # def _validate_saveframe_fields(self, nef=None):
    #     ERROR_KEY = 'SAVEFRAMES'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     for sf_name, saveframe in nef.items():
    #         e += self.__dict_missing_keys(saveframe,
    #                                       Nef.NEF_ALL_SAVEFRAME_REQUIRED_FIELDS,
    #                                       label=sf_name)
    #         e += self.__sf_framecode_name_mismatch(saveframe, sf_name)
    #     return errors
    #
    # def _validate_required_saveframes(self, nef=None):
    #     ERROR_KEY = 'REQUIRED_SAVEFRAMES'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     e += self.__dict_missing_keys(nef, Nef.NEF_REQUIRED_SAVEFRAME_BY_FRAMECODE)
    #     e += self.__dict_missing_value_with_key(nef, Nef.NEF_REQUIRED_SAVEFRAME_BY_CATEGORY)
    #     return errors
    #
    # def _validate_metadata(self, nef=None):
    #     ERROR_KEY = 'METADATA'
    #     DICT_KEY = 'nef_nmr_meta_data'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     if DICT_KEY not in nef:
    #         return {ERROR_KEY: 'No {} saveframe.'.format(DICT_KEY)}
    #     else:
    #         md = nef[DICT_KEY]
    #         e += self.__dict_missing_keys(md, Nef.MD_REQUIRED_FIELDS)
    #         e += self.__sf_framecode_name_mismatch(md, DICT_KEY)
    #         e += self.__sf_category_name_mismatch(md, DICT_KEY)
    #         e += self.__dict_nonallowed_keys(md, (Nef.MD_REQUIRED_FIELDS +
    #                                               Nef.MD_OPTIONAL_FIELDS +
    #                                               Nef.MD_OPTIONAL_LOOPS),
    #                                          label=DICT_KEY)
    #
    #         if 'format_name' in md:
    #             if md['format_name'] != 'Nmr_Exchange_Format':
    #                 e.append("format_name must be 'Nmr_Exchange_Format'.")
    #         if 'format_version' in md:
    #             major_version = str(md['format_version']).split('.')[0]
    #             if major_version != Nef.__nef_version__.split('.')[0]:
    #                 e.append('This reader does not support format version {}.'.format(major_version))
    #         if 'creation_date' in md:
    #             pass  # TODO: How to validate the creation date?
    #         if 'uuid' in md:
    #             pass  # TODO: How to validate the uuid?
    #
    #         if 'nef_related_entries' in md:
    #             for i, entry in enumerate(md['nef_related_entries']):
    #                 label = '{}:nef_related_entries entry {}'.format(DICT_KEY, i + 1)
    #                 e += self.__dict_missing_keys(entry, Nef.MD_RE_REQUIRED_FIELDS, label=label)
    #                 e += self.__dict_nonallowed_keys(entry, Nef.MD_RE_REQUIRED_FIELDS, label=label)
    #
    #         if 'nef_program_script' in md:
    #             for i, entry in enumerate(md['nef_program_script'].data):
    #                 label = '{}:nef_program_script entry {}'.format(DICT_KEY, i + 1)
    #                 e += self.__dict_missing_keys(entry, Nef.MD_PS_REQUIRED_FIELDS, label=label)
    #                 # Note: Because program specific parameters are allowed, there are not restrictions
    #                 # on what fields can be in this loop
    #             e += self.__loop_entries_inconsistent_keys(md['nef_program_script'].data,
    #                                                        label='{}:nef_program_script'.format(DICT_KEY))
    #
    #         if 'nef_run_history' in md:
    #             for i, entry in enumerate(md['nef_run_history']):
    #                 label = '{}:nef_run_history entry {}'.format(DICT_KEY, i + 1)
    #                 e += self.__dict_missing_keys(entry, Nef.MD_RH_REQUIRED_FIELDS, label=label)
    #                 e += self.__dict_nonallowed_keys(entry,
    #                                                  (Nef.MD_RH_REQUIRED_FIELDS +
    #                                                   Nef.MD_RH_OPTIONAL_FIELDS),
    #                                                  label=label)
    #             e += self.__loop_entries_inconsistent_keys(md['nef_run_history'],
    #                                                        label='{}:nef_run_history'.format(DICT_KEY))
    #
    #     return errors
    #
    # def _validate_molecular_system(self, nef=None):
    #     ERROR_KEY = 'MOLECULAR_SYSTEM'
    #     DICT_KEY = 'nef_molecular_system'
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     if 'nef_molecular_system' not in nef:
    #         return {ERROR_KEY: 'No {} saveframe.'.format(DICT_KEY)}
    #     else:
    #         ms = nef[DICT_KEY]
    #         e += self.__dict_missing_keys(ms, Nef.MS_REQUIRED_FIELDS + Nef.MS_REQUIRED_LOOPS)
    #         e += self.__dict_nonallowed_keys(ms, (Nef.MS_REQUIRED_FIELDS +
    #                                               Nef.MS_REQUIRED_LOOPS +
    #                                               Nef.MS_OPTIONAL_LOOPS),
    #                                          label=DICT_KEY)
    #         e += self.__sf_framecode_name_mismatch(ms, DICT_KEY)
    #         e += self.__sf_category_name_mismatch(ms, DICT_KEY)
    #
    #         if 'nef_sequence' in ms:
    #             if len(ms['nef_sequence'].data) == 0:
    #                 e.append('Empty nef_sequence.')
    #             else:
    #                 for i, entry in enumerate(ms['nef_sequence'].data):
    #                     label = '{}:nef_sequence entry {}'.format(DICT_KEY, i + 1)
    #                     e += self.__dict_missing_keys(entry, Nef.MS_NS_REQUIRED_FIELDS, label=label)
    #                     e += self.__dict_nonallowed_keys(entry, Nef.MS_NS_REQUIRED_FIELDS, label=label)
    #
    #         if 'nef_covalent_links' in ms:
    #             for i, entry in enumerate(ms['nef_covalent_links'].data):
    #                 label = '{}:nef_covalent_links entry {}'.format(DICT_KEY, i + 1)
    #                 e += self.__dict_missing_keys(entry, Nef.MS_CL_REQUIRED_FIELDS, label=label)
    #                 e += self.__dict_nonallowed_keys(entry, Nef.MS_CL_REQUIRED_FIELDS, label=label)
    #     return errors
    #
    # def _validate_chemical_shift_lists(self, nef=None):
    #     ERROR_KEY = 'CHEMICAL_SHIFT_LISTS'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     found_csl = False
    #     for saveframe_name, saveframe in nef.items():
    #         if 'sf_category' in saveframe:
    #             if saveframe['sf_category'] == 'nef_chemical_shift_list':
    #                 found_csl = True
    #                 e += self.__dict_missing_keys(saveframe,
    #                                               (Nef.CSL_REQUIRED_FIELDS +
    #                                                Nef.CSL_REQUIRED_LOOPS),
    #                                               label=saveframe_name)
    #                 e += self.__dict_nonallowed_keys(saveframe, (Nef.CSL_REQUIRED_FIELDS +
    #                                                              Nef.CSL_REQUIRED_LOOPS),
    #                                                  label=saveframe_name)
    #
    #                 if 'nef_chemical_shift' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_chemical_shift'].data):
    #                         label = '{}:nef_chemical_shift entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.CSL_CS_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.CSL_CS_REQUIRED_FIELDS +
    #                                                           Nef.CSL_CS_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_chemical_shift'].data,
    #                                                                label='{}:nef_chemical_shift'.format(saveframe_name))
    #     if not found_csl:
    #         e.append('No nef_chemical_shift_list saveframes found.')
    #     return errors
    #
    # def _validate_distance_restraint_lists(self, nef=None):
    #     ERROR_KEY = 'DISTANCE_RESTRAINT_LISTS'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     for saveframe_name, saveframe in nef.items():
    #         if 'sf_category' in saveframe:
    #             if saveframe['sf_category'] == 'nef_distance_restraint_list':
    #                 e += self.__dict_missing_keys(saveframe,
    #                                               (Nef.DRL_REQUIRED_FIELDS +
    #                                                Nef.DRL_REQUIRED_LOOPS),
    #                                               label=saveframe_name)
    #                 e += self.__dict_nonallowed_keys(saveframe, (Nef.DRL_REQUIRED_FIELDS +
    #                                                              Nef.DRL_REQUIRED_LOOPS +
    #                                                              Nef.DRL_OPTIONAL_FIELDS),
    #                                                  label=saveframe_name)
    #
    #                 if 'nef_distance_restraint' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_distance_restraint'].data):
    #                         label = '{}:nef_distance_restraint entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.DRL_DR_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.DRL_DR_REQUIRED_FIELDS +
    #                                                           Nef.DRL_DR_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_distance_restraint'].data,
    #                                                                label='{}:nef_distance_restraint'.format(saveframe_name))
    #     return errors
    #
    # def _validate_dihedral_restraint_lists(self, nef=None):
    #     ERROR_KEY = 'DIHEDRAL_RESTRAINT_LISTS'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     for saveframe_name, saveframe in nef.items():
    #         if 'sf_category' in saveframe:
    #             if saveframe['sf_category'] == 'nef_dihedral_restraint_list':
    #                 e += self.__dict_missing_keys(saveframe,
    #                                               (Nef.DIHRL_REQUIRED_FIELDS +
    #                                                Nef.DIHRL_REQUIRED_LOOPS),
    #                                               label=saveframe_name)
    #                 e += self.__dict_nonallowed_keys(saveframe, (Nef.DIHRL_REQUIRED_FIELDS +
    #                                                              Nef.DIHRL_REQUIRED_LOOPS +
    #                                                              Nef.DIHRL_OPTIONAL_FIELDS),
    #                                                  label=saveframe_name)
    #
    #                 if 'nef_dihedral_restraint' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_dihedral_restraint'].data):
    #                         label = '{}:nef_dihedral_restraint entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.DIHRL_DIHR_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.DIHRL_DIHR_REQUIRED_FIELDS +
    #                                                           Nef.DIHRL_DIHR_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_dihedral_restraint'].data,
    #                                                                label='{}:nef_dihedral_restraint'.format(saveframe_name))
    #     return errors
    #
    # def _validate_rdc_restraint_lists(self, nef=None):
    #     ERROR_KEY = 'RDC_RESTRAINT_LISTS'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     for saveframe_name, saveframe in nef.items():
    #         if 'sf_category' in saveframe:
    #             if saveframe['sf_category'] == 'nef_rdc_restraint_list':
    #                 e += self.__dict_missing_keys(saveframe,
    #                                               (Nef.RRL_REQUIRED_FIELDS +
    #                                                Nef.RRL_REQUIRED_LOOPS),
    #                                               label=saveframe_name)
    #                 e += self.__dict_nonallowed_keys(saveframe, (Nef.RRL_REQUIRED_FIELDS +
    #                                                              Nef.RRL_REQUIRED_LOOPS +
    #                                                              Nef.RRL_OPTIONAL_FIELDS),
    #                                                  label=saveframe_name)
    #
    #                 if 'nef_rdc_restraint' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_rdc_restraint'].data):
    #                         label = '{}:nef_rdc_restraint_list entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.RRL_RR_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.RRL_RR_REQUIRED_FIELDS +
    #                                                           Nef.RRL_RR_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_rdc_restraint'].data,
    #                                                                label='{}:nef_rdc_restraint'.format(saveframe_name))
    #     return errors
    #
    # def _validate_peak_lists(self, nef=None):
    #     ERROR_KEY = 'PEAK_LISTS'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     for saveframe_name, saveframe in nef.items():
    #         if 'sf_category' in saveframe:
    #             if saveframe['sf_category'] == 'nef_nmr_spectrum':
    #                 e += self.__dict_missing_keys(saveframe,
    #                                               (Nef.PL_REQUIRED_FIELDS +
    #                                                Nef.PL_REQUIRED_LOOPS),
    #                                               label=saveframe_name)
    #                 e += self.__dict_nonallowed_keys(saveframe, (Nef.PL_REQUIRED_FIELDS +
    #                                                              Nef.PL_REQUIRED_LOOPS +
    #                                                              Nef.PL_OPTIONAL_FIELDS),
    #                                                  label=saveframe_name)
    #                 if 'chemical_shift_list' in saveframe:
    #                     csl = saveframe['chemical_shift_list']
    #                     if csl not in nef:
    #                         e.append('{}: missing chemical_shift_list {}.'
    #                                  .format(saveframe['sf_framecode'], csl))
    #                 if 'nef_spectrum_dimension' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_spectrum_dimension'].data):
    #                         label = '{}:nef_rdc_restraint_list entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.PL_SD_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.PL_SD_REQUIRED_FIELDS +
    #                                                           Nef.PL_SD_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_spectrum_dimension'].data,
    #                                                                label='{}:nef_spectrum_dimension'.format(saveframe_name))
    #                 if 'nef_spectrum_dimension_transfer' in saveframe:
    #                     for i, entry in enumerate(saveframe['nef_spectrum_dimension_transfer'].data):
    #                         label = '{}:nef_rdc_restraint_list entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, Nef.PL_SDT_REQUIRED_FIELDS, label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.PL_SDT_REQUIRED_FIELDS +
    #                                                           Nef.PL_SDT_OPTIONAL_FIELDS),
    #                                                          label=label)
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_spectrum_dimension_transfer'].data,
    #                                                                label='{}:nef_spectrum_dimension_transfer'.format(saveframe_name))
    #                 if 'nef_peak' in saveframe:
    #                     peak_dimensions = len(saveframe['nef_spectrum_dimension'].data)
    #                     req_fields = []
    #                     opt_fields = []
    #                     for i in range(peak_dimensions):
    #                         for f in Nef.PL_P_REQUIRED_FIELDS_PATTERN:
    #                             req_fields.append(f.format(i + 1))
    #
    #                     if len(saveframe['nef_peak'].data) > 0:
    #                         for i in Nef.PL_P_REQUIRED_ALTERNATE_FIELDS:
    #                             found_alternate = False
    #                             for j in i:
    #                                 if j in saveframe['nef_peak'].data[0].keys():
    #                                     found_alternate = True
    #                                     req_fields.append(j)
    #                                     # for k in Nef.PL_P_OPTIONAL_ALTERNATE_FIELDS:
    #                                     #     req_fields.append(k.format(j))
    #                             if not found_alternate:
    #                                 e.append('test: found_alternate')
    #
    #                     for req_field in req_fields:
    #                         for optional_re, optional_re_val in Nef.PL_P_OPTIONAL_ALTERNATE_FIELDS.items():
    #                             match = re.search(optional_re, req_field)
    #                             if match:
    #                                 for orv in optional_re_val:
    #                                     opt_fields.append(orv.format(match.groups([0])[0]))
    #                     opt_fields = list(set(opt_fields))
    #                     for i, entry in enumerate(saveframe['nef_peak'].data):
    #                         label = '{}:nef_peak entry {}'.format(saveframe_name, i + 1)
    #                         e += self.__dict_missing_keys(entry, (Nef.PL_P_REQUIRED_FIELDS +
    #                                                               req_fields),
    #                                                       label=label)
    #                         e += self.__dict_nonallowed_keys(entry,
    #                                                          (Nef.PL_P_REQUIRED_FIELDS +
    #                                                           req_fields +
    #                                                           opt_fields),
    #                                                          label=label)
    #
    #                     e += self.__loop_entries_inconsistent_keys(saveframe['nef_peak'].data,
    #                                                                label='{}:nef_peak'.format(saveframe_name))
    #     return errors
    #
    # def _validate_linkage_table(self, nef=None):
    #     ERROR_KEY = 'LINKAGE_TABLES'
    #     DICT_KEY = 'nef_peak_restraint_links'
    #
    #     if nef is None:
    #         nef = self.nef
    #     errors = {ERROR_KEY: []}
    #     e = errors[ERROR_KEY]
    #
    #     if DICT_KEY in nef:
    #         prls = nef[DICT_KEY]
    #         e += self.__dict_missing_keys(prls, Nef.PRLS_REQUIRED_FIELDS)
    #         e += self.__sf_framecode_name_mismatch(prls, DICT_KEY)
    #         e += self.__sf_category_name_mismatch(prls, DICT_KEY)
    #         e += self.__dict_nonallowed_keys(prls, (Nef.PRLS_REQUIRED_FIELDS +
    #                                                 Nef.PRLS_REQUIRED_LOOPS),
    #                                          label=DICT_KEY)
    #         if 'nef_peak_restraint_link' in prls:
    #             for i, entry in enumerate(prls['nef_peak_restraint_link'].data):
    #                 label = '{}:nef_peak_restraint_link entry {}'.format(DICT_KEY, i + 1)
    #                 e += self.__dict_missing_keys(entry, Nef.PRLS_PRL_REQUIRED_FIELDS, label=label)
    #                 e += self.__dict_nonallowed_keys(entry, Nef.PRLS_PRL_REQUIRED_FIELDS, label=label)
    #
    #     return errors

    def __dict_missing_keys(self, dct, required_keys, label=None):
        if label is None:
            return ['Missing {} label.'.format(key) for key in required_keys if key not in dct]
        return ['{}: missing {} label.'.format(label, key) for key in required_keys if key not in dct]

    def __dict_missing_value_with_key(self, dct, keys):
        errors = []
        for key in keys:
            found_key = False
            for k, v in dct.items():
                if ('sf_category' in v) and (v['sf_category'] == key):
                    found_key = True
            if not found_key:
                errors.append('No saveframes with sf_category: {}.'.format(key))
        return errors

    def __sf_framecode_name_mismatch(self, dct, sf_framecode):
        if 'sf_framecode' in dct:
            if dct['sf_framecode'] != sf_framecode:
                return ["sf_framecode {} must match key {}.".format(dct['sf_framecode'], sf_framecode)]
        return []

    def __sf_category_name_mismatch(self, dct, sf_category):
        if 'sf_category' in dct:
            if dct['sf_category'] != sf_category:
                return ["sf_category {} must be {}.".format(dct['sf_category'], sf_category)]
        # else:
        #     return ["No sf_category.",]
        return []

    def __loop_entries_inconsistent_keys(self, loop, label):
        errors = []
        if len(loop) > 0:
            fields = list(loop[0].keys())
            fields_count = len(fields)
            finished = False
            while not finished:
                errors = []
                finished = True
                for i, entry in enumerate(loop):
                    for field in entry:
                        if field not in fields:
                            fields.append(field)
                    if len(fields) > fields_count:
                        fields_count = len(fields)
                        finished = False
                        break
                    errors += self.__dict_missing_keys(entry, fields, label=label + ' item {}'
                                                       .format(i))
        return errors

    def __dict_nonallowed_keys(self, dct, allowed_keys, label=None):
        if label is None:
            return ["Field '{}' not allowed.".format(key)
                    for key in dct.keys() if key not in allowed_keys and not key.startswith(CCPN_PREFIX)]
        return ["Field '{}' not allowed in {}.".format(key, label)
                for key in dct.keys() if key not in allowed_keys and not key.startswith(CCPN_PREFIX)]
