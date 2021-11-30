#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ICEscreen copyright Université de Lorraine - INRAE
# This file is part of ICEscreen.
# ICEscreen is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License Affero
# as published by the Free Software Foundation version 3 of the License.
# ICEscreen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License Affero for more details.
# You should have received a copy of the GNU General Public License Affero
# along with ICEscreen. If not, see <https://www.gnu.org/licenses/>.

"""
@author: tlacroix, jlao
"""
# import specific class OO for this script
import EMStructure


# Once an ICE or IME structure have been established, assign its type
def assignPutativeTypeStructure(ICEsIMEsStructureSent, maxNumberCDSForFilterIMESize):

    if ICEsIMEsStructureSent.delMerging_idxListUpstreamStructure >= 0:
        strTypeToSet = "deleted after merging event"
        ICEsIMEsStructureSent.catStructConjModule = strTypeToSet
    else:

        greenLightDistanceIME = ICEsIMEsStructureSent.isFilterIMESizeOk(maxNumberCDSForFilterIMESize)
        # totalNumberSPSureAndNotSure = len(ICEsIMEsStructureSent.listOrderedSPs) + len(ICEsIMEsStructureSent.setIntegraseLocusTagsToManuallyCheck) + len(ICEsIMEsStructureSent.setSPConjModuleLocusTagsToManuallyCheck)

        putativePartialICEVirB4 = False
        if len(ICEsIMEsStructureSent.setSPInConflict) >= 1:
            for currSPInConflict in ICEsIMEsStructureSent.setSPInConflict:
                if currSPInConflict.SPType == "VirB4":
                    strListInternalIdentifier = EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(currSPInConflict.setICEsIMEsStructureInConflict)
                    listInternalIdentifier = strListInternalIdentifier.split(", ")
                    if listInternalIdentifier[0] == str(ICEsIMEsStructureSent.internalIdentifier):
                        putativePartialICEVirB4 = True  # seemingly attribute virB4 in conflict to at least 1 structure

        # deal with main type
        strTypeToSet = ""
        if len(ICEsIMEsStructureSent.listRelaxase) == 1 and len(ICEsIMEsStructureSent.listCouplingProtein) == 1 and len(ICEsIMEsStructureSent.listVirB4) == 1:
            strTypeToSet = "ICE"  # canonical conjugaison module
        elif len(ICEsIMEsStructureSent.listRelaxase) >= 1 and len(ICEsIMEsStructureSent.listCouplingProtein) >= 1 and len(ICEsIMEsStructureSent.listVirB4) >= 1:
            strTypeToSet = "ICE"  # mulitple adjacent SPs conjugaison module
        elif len(ICEsIMEsStructureSent.listVirB4) >= 1 or putativePartialICEVirB4:
            strTypeToSet = "partial ICE VirB4"
        elif len(ICEsIMEsStructureSent.listVirB4) == 0 and len(ICEsIMEsStructureSent.listOrderedSPs) >= 2 and greenLightDistanceIME:
            strTypeToSet = "IME"
        elif len(ICEsIMEsStructureSent.listRelaxase) >= 1 or len(ICEsIMEsStructureSent.listCouplingProtein) >= 1:
            strTypeToSet = "partial conj module"
        else:
            strTypeToSet = "unsure"

        if strTypeToSet:
            ICEsIMEsStructureSent.catStructConjModule = strTypeToSet
        else:
            raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized category conj module for structure = {}".format(
                    str(ICEsIMEsStructureSent.internalIdentifier)))

        # deal with integrase
        strTypeToSet = ""
        if len(ICEsIMEsStructureSent.listIntegraseUpstream) == 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0:
            if ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType == "Tyrosine integrase":
                strTypeToSet = "one integrase - Tyr"
            elif ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType == "Serine integrase":
                strTypeToSet = "one integrase - Ser"
            elif ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType == "DDE transposase":
                strTypeToSet = "one integrase - DDE"
            else:
                raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized listIntegraseUpstream[0].SPType {} for locus tag = {}".format(
                        ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType,
                        ICEsIMEsStructureSent.listIntegraseUpstream[0].locusTag))
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 1:
            if ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "Tyrosine integrase":
                strTypeToSet = "one integrase - Tyr"
            elif ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "Serine integrase":
                strTypeToSet = "one integrase - Ser"
            elif ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "DDE transposase":
                strTypeToSet = "one integrase - DDE"
            else:
                raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized listIntegraseDownstream[0].SPType {} for locus tag = {}".format(
                        ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType,
                        ICEsIMEsStructureSent.listIntegraseDownstream[0].locusTag))
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) > 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0:
            setIntegraseType = set()
            for currIntegrase in ICEsIMEsStructureSent.listIntegraseUpstream:
                if currIntegrase.SPType == "Tyrosine integrase":
                    setIntegraseType.add("Tyr")
                elif currIntegrase.SPType == "Serine integrase":
                    setIntegraseType.add("Ser")
                elif currIntegrase.SPType == "DDE transposase":
                    setIntegraseType.add("DDE")
                else:
                    raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized currIntegrase.SPType {} for locus tag = {}".format(
                            currIntegrase.SPType, currIntegrase.locusTag))
            if len(setIntegraseType) == 1:
                strTypeToSet = "mulitple adjacent integrases - " + setIntegraseType.pop()
            elif len(setIntegraseType) == 0:
                raise RuntimeError("Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseUpstream) > 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0: len (setIntegraseType) == 0 for EM structure = {}".format(
                        str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "mulitple adjacent integrases - Mixed"
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureSent.listIntegraseDownstream) > 1:
            setIntegraseType = set()
            for currIntegrase in ICEsIMEsStructureSent.listIntegraseDownstream:
                if currIntegrase.SPType == "Tyrosine integrase":
                    setIntegraseType.add("Tyr")
                elif currIntegrase.SPType == "Serine integrase":
                    setIntegraseType.add("Ser")
                elif currIntegrase.SPType == "DDE transposase":
                    setIntegraseType.add("DDE")
                else:
                    raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized currIntegrase.SPType {} for locus tag = {}".format(
                            currIntegrase.SPType, currIntegrase.locusTag))
            if len(setIntegraseType) == 1:
                strTypeToSet = "mulitple adjacent integrases - " + setIntegraseType.pop()
            elif len(setIntegraseType) == 0:
                raise RuntimeError("Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseDownstream) > 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0: len (setIntegraseType) == 0 for EM structure = {}".format(
                        str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "mulitple adjacent integrases - Mixed"
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) >= 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) >= 1:
            setIntegraseType = set()
            for currIntegrase in ICEsIMEsStructureSent.listIntegraseUpstream:
                if currIntegrase.SPType == "Tyrosine integrase":
                    setIntegraseType.add("Tyr")
                elif currIntegrase.SPType == "Serine integrase":
                    setIntegraseType.add("Ser")
                elif currIntegrase.SPType == "DDE transposase":
                    setIntegraseType.add("DDE")
                else:
                    raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized currIntegrase.SPType {} for locus tag = {}".format(
                            currIntegrase.SPType, currIntegrase.locusTag))
            for currIntegrase in ICEsIMEsStructureSent.listIntegraseDownstream:
                if currIntegrase.SPType == "Tyrosine integrase":
                    setIntegraseType.add("Tyr")
                elif currIntegrase.SPType == "Serine integrase":
                    setIntegraseType.add("Ser")
                elif currIntegrase.SPType == "DDE transposase":
                    setIntegraseType.add("DDE")
                else:
                    raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized currIntegrase.SPType {} for locus tag = {}".format(
                            currIntegrase.SPType, currIntegrase.locusTag))
            if len(setIntegraseType) == 1:
                strTypeToSet = "both upstream or downstream integrases - " + setIntegraseType.pop()
            elif len(setIntegraseType) == 0:
                raise RuntimeError(
                        "Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseUpstream) >= 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) >= 1: len (setIntegraseType) == 0 for EM structure = {}".format(
                                str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "both upstream or downstream integrases - Mixed"
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0:
            strTypeToSet = "no integrase"
        else:
            raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized category for integrase for EM structure = {}".format(
                    str(ICEsIMEsStructureSent.internalIdentifier)))
        ICEsIMEsStructureSent.categoryRegardingIntegrase = strTypeToSet

        # fill in catStructWholeElem
        strTypeToSet = ""
        (strTypeToSet,
         is01_ICETyr,
         is01_ICESer,
         is01_ICEDDE,
         is01_ConjugationModulesNoIntegrase,
         is01_PartialICETyr,
         is01_PartialICESer,
         is01_PartialICEDDE,
         is01_PartialICENoIntegrase,
         is01_IMETyr,
         is01_IMESer,
         is01_IMEDDE,
         is01_MobilizableElementsNoIntegrase,
         is01_PartialConjModuleTyr,
         is01_PartialConjModuleSer,
         is01_PartialConjModuleDDE,
         is01_PartialConjModuleNoIntegrase,
         is01_Unsure) = getcatStructWholeElem(ICEsIMEsStructureSent, maxNumberCDSForFilterIMESize)
        ICEsIMEsStructureSent.catStructWholeElem = strTypeToSet


def getcatStructWholeElem(ICEIMEStructureSent, maxNumberCDSForFilterIMESize):
    strTypeToSet = ""
    is01_ICETyr = 0
    is01_ICESer = 0
    is01_ICEDDE = 0
    is01_ConjugationModulesNoIntegrase = 0
    is01_PartialICETyr = 0
    is01_PartialICESer = 0
    is01_PartialICEDDE = 0
    is01_PartialICENoIntegrase = 0
    is01_IMETyr = 0
    is01_IMESer = 0
    is01_IMEDDE = 0
    is01_MobilizableElementsNoIntegrase = 0
    is01_PartialConjModuleTyr = 0
    is01_PartialConjModuleSer = 0
    is01_PartialConjModuleDDE = 0
    is01_PartialConjModuleNoIntegrase = 0
    is01_Unsure = 0

    if ICEIMEStructureSent.catStructConjModule == "ICE":
        if ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Tyr":
            is01_ICETyr = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Ser":
            is01_ICESer = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - DDE":
            is01_ICEDDE = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Tyr":
            is01_ICETyr = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Ser":
            is01_ICESer = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - DDE":
            is01_ICEDDE = 1
            strTypeToSet = "Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Tyr":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases - Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Ser":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases - Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - DDE":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases - DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Mixed":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases - Mixed found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "no integrase":
            is01_ConjugationModulesNoIntegrase = 1
            strTypeToSet = "Conjugation module (R+C+V)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))
    elif ICEIMEStructureSent.catStructConjModule == "partial ICE VirB4":
        if ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Tyr":
            is01_PartialICETyr = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Ser":
            is01_PartialICESer = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - DDE":
            is01_PartialICEDDE = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Tyr":
            is01_PartialICETyr = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Ser":
            is01_PartialICESer = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - DDE":
            is01_PartialICEDDE = 1
            strTypeToSet = "Partial ICE (at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Tyr":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases - Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Ser":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases - Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - DDE":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases - DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Mixed":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases - Mixed found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "no integrase":
            is01_PartialICENoIntegrase = 1
            strTypeToSet = "Partial ICE (at least V)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))
    elif ICEIMEStructureSent.catStructConjModule == "IME":
        if ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Tyr":
            is01_IMETyr = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Ser":
            is01_IMESer = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - DDE":
            is01_IMEDDE = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Tyr":
            is01_IMETyr = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Ser":
            is01_IMESer = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - DDE":
            is01_IMEDDE = 1
            strTypeToSet = "Complete IME (R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Tyr":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases - Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Ser":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases - Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - DDE":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases - DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Mixed":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases - Mixed found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "no integrase":
            is01_MobilizableElementsNoIntegrase = 1
            strTypeToSet = "Mobilizable element (R+C with size <= " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))
    elif ICEIMEStructureSent.catStructConjModule == "partial conj module":
        if ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Tyr":
            is01_PartialConjModuleTyr = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - Ser":
            is01_PartialConjModuleSer = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "one integrase - DDE":
            is01_PartialConjModuleDDE = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Tyr":
            is01_PartialConjModuleTyr = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - Ser":
            is01_PartialConjModuleSer = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "mulitple adjacent integrases - DDE":
            is01_PartialConjModuleDDE = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Tyr":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases - Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Ser":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases - Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - DDE":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases - DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "both upstream or downstream integrases - Mixed":
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases - Mixed found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase == "no integrase":
            is01_PartialConjModuleNoIntegrase = 1
            strTypeToSet = "Other partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))
    elif ICEIMEStructureSent.catStructConjModule == "unsure":
        is01_Unsure = 1
        strTypeToSet = "Unsure, to check"
    else:
        raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.catStructConjModule = {} for ICE IME {}".format(
                ICEIMEStructureSent.catStructConjModule,
                ICEIMEStructureSent.internalIdentifier))

    return (strTypeToSet,
            is01_ICETyr,
            is01_ICESer,
            is01_ICEDDE,
            is01_ConjugationModulesNoIntegrase,
            is01_PartialICETyr,
            is01_PartialICESer,
            is01_PartialICEDDE,
            is01_PartialICENoIntegrase,
            is01_IMETyr,
            is01_IMESer,
            is01_IMEDDE,
            is01_MobilizableElementsNoIntegrase,
            is01_PartialConjModuleTyr,
            is01_PartialConjModuleSer,
            is01_PartialConjModuleDDE,
            is01_PartialConjModuleNoIntegrase,
            is01_Unsure)


def printOverallStatsToSummaryFile(
        totalNumberSP,
        totalNumberIntegrase,
        totalNumberUnassignedIntegrase,
        totalNumberRelaxase,
        totalNumberUnassignedRelaxase,
        totalNumberCoupling,
        totalNumberUnassignedCoupling,
        totalNumberVirb4,
        totalNumberUnassignedVirb4,
        listOfListAllICEsIMEsStructure,
        maxNumberCDSForSplitSPsByColocalizion,
        maxNumberCDSForFilterIMESize,
        summaryFile):

    # print("\n** Script version: \"{}\"".format(scriptVersion), file=summaryFile)
    # print("\n** Parameters:", file=summaryFile)
    # print("Maximum CDS numbers in initial segments (maxNumberCDSForSplitSPsByColocalizion): \"{}\"".format(maxNumberCDSForSplitSPsByColocalizion), file=summaryFile)
    # print("Maximum CDS numbers in an IME element (maxNumberCDSForFilterIMESize): \"{}\"".format(maxNumberCDSForFilterIMESize), file=summaryFile)

    totalNumberOfSegments = 0
    totalNumberOfSegmentsWithOneElement = 0
    totalNumberOfSegmentsWithSeveralElements = 0
    totalNumberOfSegmentsWithNestedElements = 0
    totalNumberOfSegmentsWithNoElement = 0
    countICETyrTotal = 0
    countICESerTotal = 0
    countICEDDETotal = 0
    countConjugationModulesNoIntegraseTotal = 0
    countPartialICETyrTotal = 0
    countPartialICESerTotal = 0
    countPartialICEDDETotal = 0
    countPartialICENoIntegraseTotal = 0
    countIMETyrTotal = 0
    countIMESerTotal = 0
    countIMEDDETotal = 0
    countMobilizableElementsNoIntegraseTotal = 0
    countPartialConjModuleTyrTotal = 0
    countPartialConjModuleSerTotal = 0
    countPartialConjModuleDDETotal = 0
    countPartialConjModuleNoIntegraseTotal = 0
    countUnsureTotal = 0
    countIntegrasesNotAttributedTyrTotal = 0
    countIntegrasesNotAttributedSerTotal = 0
    countIntegrasesNotAttributedDDETotal = 0
    countEMHostOnlyTotal = 0
    countEMGuestOnlyTotal = 0
    countEMHostAndGuestTotal = 0

    for idx, currListAllICEsIMEsStructure in enumerate(listOfListAllICEsIMEsStructure):

        countICETyrSegment = 0
        countICESerSegment = 0
        countICEDDESegment = 0
        countConjugationModulesNoIntegraseSegment = 0
        countPartialICETyrSegment = 0
        countPartialICESerSegment = 0
        countPartialICEDDESegment = 0
        countPartialICENoIntegraseSegment = 0
        countIMETyrSegment = 0
        countIMESerSegment = 0
        countIMEDDESegment = 0
        countMobilizableElementsNoIntegraseSegment = 0
        countPartialConjModuleTyrSegment = 0
        countPartialConjModuleSerSegment = 0
        countPartialConjModuleDDESegment = 0
        countPartialConjModuleNoIntegraseSegment = 0
        countUnsureSegment = 0
        countIntegrasesNotAttributedTyrSegment = 0
        countIntegrasesNotAttributedSerSegment = 0
        countIntegrasesNotAttributedDDESegment = 0
        countEMHostOnlySegment = 0
        countEMGuestOnlySegment = 0
        countEMHostAndGuestSegment = 0

        for currICEIMEStructure in currListAllICEsIMEsStructure:
            # print stat for each genomic segment and total

            if len(currICEIMEStructure.setHostNestedICEsIMEsStructure) >= 1 and len(currICEIMEStructure.setGuestsNestedICEsIMEsStructure) >= 1:
                countEMHostAndGuestSegment += 1
            elif len(currICEIMEStructure.setHostNestedICEsIMEsStructure) >= 1:
                countEMGuestOnlySegment += 1
            elif len(currICEIMEStructure.setGuestsNestedICEsIMEsStructure) >= 1:
                countEMHostOnlySegment += 1

            (strTypeToSet,
             is01_ICETyr,
             is01_ICESer,
             is01_ICEDDE,
             is01_ConjugationModulesNoIntegrase,
             is01_PartialICETyr,
             is01_PartialICESer,
             is01_PartialICEDDE,
             is01_PartialICENoIntegrase,
             is01_IMETyr,
             is01_IMESer,
             is01_IMEDDE,
             is01_MobilizableElementsNoIntegrase,
             is01_PartialConjModuleTyr,
             is01_PartialConjModuleSer,
             is01_PartialConjModuleDDE,
             is01_PartialConjModuleNoIntegrase,
             is01_Unsure) = getcatStructWholeElem(currICEIMEStructure, maxNumberCDSForFilterIMESize)
            countICETyrSegment += is01_ICETyr
            countICESerSegment += is01_ICESer
            countICEDDESegment += is01_ICEDDE
            countConjugationModulesNoIntegraseSegment += is01_ConjugationModulesNoIntegrase
            countPartialICETyrSegment += is01_PartialICETyr
            countPartialICESerSegment += is01_PartialICESer
            countPartialICEDDESegment += is01_PartialICEDDE
            countPartialICENoIntegraseSegment += is01_PartialICENoIntegrase
            countIMETyrSegment += is01_IMETyr
            countIMESerSegment += is01_IMESer
            countIMEDDESegment += is01_IMEDDE
            countMobilizableElementsNoIntegraseSegment += is01_MobilizableElementsNoIntegrase
            countPartialConjModuleTyrSegment += is01_PartialConjModuleTyr
            countPartialConjModuleSerSegment += is01_PartialConjModuleSer
            countPartialConjModuleDDESegment += is01_PartialConjModuleDDE
            countPartialConjModuleNoIntegraseSegment += is01_PartialConjModuleNoIntegrase
            countUnsureSegment += is01_Unsure

        totalNumberOfSegments += 1
        totalElementStructsInSegment = countICETyrSegment + countICESerSegment\
            + countICEDDESegment + countConjugationModulesNoIntegraseSegment\
            + countPartialICETyrSegment + countPartialICESerSegment\
            + countPartialICEDDESegment + countPartialICENoIntegraseSegment\
            + countIMETyrSegment + countIMESerSegment + countIMEDDESegment\
            + countMobilizableElementsNoIntegraseSegment\
            + countPartialConjModuleTyrSegment + countPartialConjModuleSerSegment\
            + countPartialConjModuleDDESegment\
            + countPartialConjModuleNoIntegraseSegment
        if (totalElementStructsInSegment) == 0:
            totalNumberOfSegmentsWithNoElement += 1
        elif (totalElementStructsInSegment) == 1:
            totalNumberOfSegmentsWithOneElement += 1
        else:
            totalNumberOfSegmentsWithSeveralElements += 1
        if countEMHostOnlySegment > 0:
            totalNumberOfSegmentsWithNestedElements += 1

        countICETyrTotal += countICETyrSegment
        countICESerTotal += countICESerSegment
        countICEDDETotal += countICEDDESegment
        countConjugationModulesNoIntegraseTotal += countConjugationModulesNoIntegraseSegment
        countPartialICETyrTotal += countPartialICETyrSegment
        countPartialICESerTotal += countPartialICESerSegment
        countPartialICEDDETotal += countPartialICEDDESegment
        countPartialICENoIntegraseTotal += countPartialICENoIntegraseSegment
        countIMETyrTotal += countIMETyrSegment
        countIMESerTotal += countIMESerSegment
        countIMEDDETotal += countIMEDDESegment
        countMobilizableElementsNoIntegraseTotal += countMobilizableElementsNoIntegraseSegment
        countPartialConjModuleTyrTotal += countPartialConjModuleTyrSegment
        countPartialConjModuleSerTotal += countPartialConjModuleSerSegment
        countPartialConjModuleDDETotal += countPartialConjModuleDDESegment
        countPartialConjModuleNoIntegraseTotal += countPartialConjModuleNoIntegraseSegment
        countUnsureTotal += countUnsureSegment
        countIntegrasesNotAttributedTyrTotal += countIntegrasesNotAttributedTyrSegment
        countIntegrasesNotAttributedSerTotal += countIntegrasesNotAttributedSerSegment
        countIntegrasesNotAttributedDDETotal += countIntegrasesNotAttributedDDESegment
        countEMHostOnlyTotal += countEMHostOnlySegment
        countEMGuestOnlyTotal += countEMGuestOnlySegment
        countEMHostAndGuestTotal += countEMHostAndGuestSegment

    # Elements
    totalNumberOfCompleteICE = countICETyrTotal + countICESerTotal + countICEDDETotal
    totalNumberOfConjugationModules = countConjugationModulesNoIntegraseTotal
    totalNumberOfPartialICE = countPartialICETyrTotal + countPartialICESerTotal + countPartialICEDDETotal + countPartialICENoIntegraseTotal
    totalNumberOfCompleteIME = countIMETyrTotal + countIMESerTotal + countIMEDDETotal
    totalNumberOfMobilizableElements = countMobilizableElementsNoIntegraseTotal
    totalNumberOfOtherPartialElements = countPartialConjModuleTyrTotal + countPartialConjModuleSerTotal + countPartialConjModuleDDETotal + countPartialConjModuleNoIntegraseTotal
    totalNumberOfNestedElementsInSegments = countEMHostOnlyTotal + countEMGuestOnlyTotal + countEMHostAndGuestTotal
    totalNumberOfHostElements = countEMHostOnlyTotal
    totalNumberOfGuestElements = countEMGuestOnlyTotal
    totalNumberOfHostAndGuestElements = countEMHostAndGuestTotal

    # print("\n** Overall statistics:", file=summaryFile)
    strTotalToPrint = (
        "##################### Mobile elements detection statistics #####################" +
        "\n# Number of complete elements" +
        "\nComplete ICE (4 types of SP): " + str(totalNumberOfCompleteICE) +
        "\nComplete IME (R+I, R+C+I): " + str(totalNumberOfCompleteIME) +
        "\n\n# Number of complete modules" +
        "\nConjugation module (R+C+V): " + str(totalNumberOfConjugationModules) +
        "\nMobilizable element (R+C with size <= " + str(maxNumberCDSForFilterIMESize) + " CDS): " + str(totalNumberOfMobilizableElements) +
        "\n\n# Number of partial elements" +
        "\nPartial ICE (at least V): " + str(totalNumberOfPartialICE) +
        "\nOther partial element (R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS): " + str(totalNumberOfOtherPartialElements) +
        "\n\n# Composite elements" +
        "\nTotal nested elements (partial or complete): " + str(totalNumberOfNestedElementsInSegments) +
        "\nHost element: " + str(totalNumberOfHostElements) +
        "\nGuest element: " + str(totalNumberOfGuestElements) +
        "\nElement that are both host and guest: " + str(totalNumberOfHostAndGuestElements) +
        "\n\n################### Signature proteins detection statistics ####################" +
        "\n# Total SP detected" +
        "\nTotal SP: " + str(totalNumberSP) +
        "\nTotal Integrase: " + str(totalNumberIntegrase) +
        "\nTotal Relaxase: " + str(totalNumberRelaxase) +
        "\nTotal Coupling protein: " + str(totalNumberCoupling) +
        "\nTotal VirB4: " + str(totalNumberVirb4) +
        "\n\n# Unassigned SP" +
        "\nUnassigned SP: " + str(totalNumberUnassignedIntegrase + totalNumberUnassignedRelaxase + totalNumberUnassignedCoupling + totalNumberUnassignedVirb4) +
        "\nUnassigned Integrase: " + str(totalNumberUnassignedIntegrase) +
        "\nUnassigned Relaxase: " + str(totalNumberUnassignedRelaxase) +
        "\nUnassigned Coupling protein: " + str(totalNumberUnassignedCoupling) +
        "\nUnassigned VirB4: " + str(totalNumberUnassignedVirb4) +
        "\n\n############################# Segments statistics ##############################" +
        "\nNumber of segments: " + str(totalNumberOfSegments) +
        "\nSegment with one element: " + str(totalNumberOfSegmentsWithOneElement) +
        "\nSegment with several elements: " + str(totalNumberOfSegmentsWithSeveralElements) +
        "\nSegment with nested elements: " + str(totalNumberOfSegmentsWithNestedElements) +
        "\nSegment with no element: " + str(totalNumberOfSegmentsWithNoElement) +
        "\n\n############################### Parameters used ################################" +
        "\nMaximum CDS numbers in initial segments: " + str(maxNumberCDSForSplitSPsByColocalizion) +
        "\nMaximum CDS numbers in an IME element: " + str(maxNumberCDSForFilterIMESize)
    )

    print(strTotalToPrint, file=summaryFile)
