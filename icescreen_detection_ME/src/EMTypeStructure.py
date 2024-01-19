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
@author: thomas Lacroix
"""
# import specific class OO for this script
import EMStructure
import hit
import icescreen_OO
import commonMethods

#import datetime

# Once an ICE or IME structure have been established, assign its type
def assignPutativeTypeStructure(
        ICEsIMEsStructureSent
        ):

    if ICEsIMEsStructureSent.delMerging_idxListUpstreamStructure >= 0:
        strTypeToSet = "deleted after merging event"
        ICEsIMEsStructureSent.catStructConjModule = strTypeToSet
    else:

        greenLightDistanceIME = ICEsIMEsStructureSent.isFilterIMESizeOk()

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

        ICEsIMEsStructureContainsAllTypesOfSPConjModule = True
        for typeSPConjModuleIT in icescreen_OO.listTypeSPConjModule:
            if len(ICEsIMEsStructureSent.TypeSPConjModule2listSP[typeSPConjModuleIT]) == 0:
                ICEsIMEsStructureContainsAllTypesOfSPConjModule = False
                break
        if ICEsIMEsStructureContainsAllTypesOfSPConjModule:
            strTypeToSet = "ICE"
        elif len(ICEsIMEsStructureSent.TypeSPConjModule2listSP["VirB4"]) >= 1 or putativePartialICEVirB4:
            strTypeToSet = "partial ICE VirB4"
        elif len(ICEsIMEsStructureSent.TypeSPConjModule2listSP["VirB4"]) == 0 and len(ICEsIMEsStructureSent.listOrderedSPs) >= 2 and greenLightDistanceIME:
            strTypeToSet = "IME"
        elif len(ICEsIMEsStructureSent.TypeSPConjModule2listSP["Relaxase"]) >= 1 or len(ICEsIMEsStructureSent.TypeSPConjModule2listSP["Coupling protein"]) >= 1:
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
                strTypeToSet = "one integrase Tyr"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseUpstream)
            elif ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType == "Serine integrase":
                strTypeToSet = "one integrase Ser"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseUpstream)
            elif ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType == "DDE transposase":
                strTypeToSet = "one integrase DDE"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseUpstream)
            else:
                raise RuntimeError("Error in assignPutativeTypeStructure: unrecognized listIntegraseUpstream[0].SPType {} for locus tag = {}".format(
                        ICEsIMEsStructureSent.listIntegraseUpstream[0].SPType,
                        ICEsIMEsStructureSent.listIntegraseUpstream[0].locusTag))
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 1:
            if ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "Tyrosine integrase":
                strTypeToSet = "one integrase Tyr"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseDownstream)
            elif ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "Serine integrase":
                strTypeToSet = "one integrase Ser"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseDownstream)
            elif ICEsIMEsStructureSent.listIntegraseDownstream[0].SPType == "DDE transposase":
                strTypeToSet = "one integrase DDE"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseDownstream)
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
                strTypeToSet = "multiple adjacent integrases " + setIntegraseType.pop()+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseUpstream)
            elif len(setIntegraseType) == 0:
                raise RuntimeError("Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseUpstream) > 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0: len (setIntegraseType) == 0 for EM structure = {}".format(
                        str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "multiple adjacent integrases various types"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseUpstream)
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
                strTypeToSet = "multiple adjacent integrases " + setIntegraseType.pop()+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseDownstream)
            elif len(setIntegraseType) == 0:
                raise RuntimeError("Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseDownstream) > 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0: len (setIntegraseType) == 0 for EM structure = {}".format(
                        str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "multiple adjacent integrases various types"+hit.ListSPs.stringListSPHasPseudo(ICEsIMEsStructureSent.listIntegraseDownstream)
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
            
            listAllUpAndDownIntegrasesToSend = []
            listAllUpAndDownIntegrasesToSend.extend(ICEsIMEsStructureSent.listIntegraseUpstream)
            listAllUpAndDownIntegrasesToSend.extend(ICEsIMEsStructureSent.listIntegraseDownstream)
            if len(setIntegraseType) == 1:
                strTypeToSet = "both upstream or downstream integrases " + setIntegraseType.pop()+hit.ListSPs.stringListSPHasPseudo(listAllUpAndDownIntegrasesToSend)
            elif len(setIntegraseType) == 0:
                raise RuntimeError(
                        "Error in assignPutativeTypeStructure len(ICEsIMEsStructureSent.listIntegraseUpstream) >= 1 and len(ICEsIMEsStructureSent.listIntegraseDownstream) >= 1: len (setIntegraseType) == 0 for EM structure = {}".format(
                                str(ICEsIMEsStructureSent.internalIdentifier)))
            else:
                strTypeToSet = "both upstream or downstream integrases various types"+hit.ListSPs.stringListSPHasPseudo(listAllUpAndDownIntegrasesToSend)
        elif len(ICEsIMEsStructureSent.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureSent.listIntegraseDownstream) == 0:
            strTypeToSet = "no integrase assigned"
            if len(ICEsIMEsStructureSent.setIntegraseToManuallyCheck) > 0:
                strTypeToSet += " for sure, "+str(len(ICEsIMEsStructureSent.setIntegraseToManuallyCheck))+" to check manually"
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
         is01_Unsure) = getcatStructWholeElem(
             ICEsIMEsStructureSent
             )
        ICEsIMEsStructureSent.catStructWholeElem = strTypeToSet


def getcatStructWholeElem(
        ICEIMEStructureSent
        ):
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
        if ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Tyr"):
            is01_ICETyr = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Ser"):
            is01_ICESer = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase DDE"):
            is01_ICEDDE = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Tyr"):
            is01_ICETyr = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Ser"):
            is01_ICESer = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases DDE"):
            is01_ICEDDE = 1
            strTypeToSet = "Complete ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #"Complete ICE (4 types of SP)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Tyr"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Ser"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases DDE"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases various types"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile ICE: both upstream or downstream integrases various types found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))

        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("no integrase assigned"):
            is01_ConjugationModulesNoIntegrase = 1
            strTypeToSet = "Conjugation module"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #R+C+V)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))
    elif ICEIMEStructureSent.catStructConjModule == "partial ICE VirB4":
        if ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Tyr"):
            is01_PartialICETyr = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Ser"):
            is01_PartialICESer = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase DDE"):
            is01_PartialICEDDE = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Tyr"):
            is01_PartialICETyr = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Ser"):
            is01_PartialICESer = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases DDE"):
            is01_PartialICEDDE = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Tyr"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Ser"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases DDE"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases various types"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial ICE VirB4: both upstream or downstream integrases various types found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("no integrase assigned"):
            is01_PartialICENoIntegrase = 1
            strTypeToSet = "Partial ICE"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, "")+")" #at least V)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))

    elif ICEIMEStructureSent.catStructConjModule == "IME":
        if ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Tyr"):
            is01_IMETyr = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Ser"):
            is01_IMESer = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase DDE"):
            is01_IMEDDE = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Tyr"):
            is01_IMETyr = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Ser"):
            is01_IMESer = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases DDE"):
            is01_IMEDDE = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Complete IME"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Tyr"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Ser"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases DDE"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases various types"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile IME: both upstream or downstream integrases various types found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))

        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("no integrase assigned"):
            is01_MobilizableElementsNoIntegrase = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Mobilizable element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #R+I, R+C+I)"(maxNumberCDSForFilterIMESize) + " CDS)"
        else:
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile: unrecognized ICEIMEStructureSent.categoryRegardingIntegrase = {} for ICE IME {}".format(
                    ICEIMEStructureSent.categoryRegardingIntegrase,
                    ICEIMEStructureSent.internalIdentifier))

    elif ICEIMEStructureSent.catStructConjModule == "partial conj module":
        if ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Tyr"):
            is01_PartialConjModuleTyr = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase Ser"):
            is01_PartialConjModuleSer = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("one integrase DDE"):
            is01_PartialConjModuleDDE = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Tyr"):
            is01_PartialConjModuleTyr = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases Ser"):
            is01_PartialConjModuleSer = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("multiple adjacent integrases DDE"):
            is01_PartialConjModuleDDE = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Tyr"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases Tyr found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases Ser"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases Ser found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases DDE"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases DDE found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("both upstream or downstream integrases various types"):
            raise RuntimeError("Error in printCountCategoriesICEsIMEsStructuresToOutputFile partial conj module: both upstream or downstream integrases various types found for ICE IME {}".format(
                    ICEIMEStructureSent.internalIdentifier))
        elif ICEIMEStructureSent.categoryRegardingIntegrase.startswith("no integrase assigned"):
            is01_PartialConjModuleNoIntegrase = 1
            strToAppendIfMoreThanOneSP = " with # CDS between SP > "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)
            strTypeToSet = "Other partial element"+hit.ListSPs.stringListSPHasPseudo(ICEIMEStructureSent.listOrderedSPs)+" ("+hit.ListSPs.PrintICElineFormat(ICEIMEStructureSent.listOrderedSPs, True, strToAppendIfMoreThanOneSP)+")" #(R+C, R+V, V+C with size > " + str(maxNumberCDSForFilterIMESize) + " CDS)"
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
        dictGenomeAccnum2totalNumberSP
        , dictGenomeAccnum2totalNumberIntegrase
        , dictGenomeAccnum2totalNumberUnassignedIntegrase
        , dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP
        , dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP_unaffected
        , dictGenomeAccnum2totalsetFragmentedSP
        , listOfListAllICEsIMEsStructure
        , summaryFile
        , listGenomeAccessionFromGenbankFile
        ):
    
    dictGenomeAccnum2totalNumberOfSegmentsWithOneElement = {}
    dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements = {}
    dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements = {}
    dictGenomeAccnum2totalNumberOfSegmentsWithNoElement = {}
    dictGenomeAccnum2countICETyrTotal = {}
    dictGenomeAccnum2countICESerTotal = {}
    dictGenomeAccnum2countICEDDETotal = {}
    dictGenomeAccnum2countConjugationModulesNoIntegraseTotal = {}
    dictGenomeAccnum2countPartialICETyrTotal = {}
    dictGenomeAccnum2countPartialICESerTotal = {}
    dictGenomeAccnum2countPartialICEDDETotal = {}
    dictGenomeAccnum2countPartialICENoIntegraseTotal = {}
    dictGenomeAccnum2countIMETyrTotal = {}
    dictGenomeAccnum2countIMESerTotal = {}
    dictGenomeAccnum2countIMEDDETotal = {}
    dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal = {}
    dictGenomeAccnum2countPartialConjModuleTyrTotal = {}
    dictGenomeAccnum2countPartialConjModuleSerTotal = {}
    dictGenomeAccnum2countPartialConjModuleDDETotal = {}
    dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal = {}
    dictGenomeAccnum2countUnsureTotal = {}
    dictGenomeAccnum2countIntegrasesNotAttributedTyrTotal = {}
    dictGenomeAccnum2countIntegrasesNotAttributedSerTotal = {}
    dictGenomeAccnum2countIntegrasesNotAttributedDDETotal = {}
    dictGenomeAccnum2countEMHostOnlyTotal = {}
    dictGenomeAccnum2countEMGuestOnlyTotal = {}
    dictGenomeAccnum2countEMHostAndGuestTotal = {}

    dictGenomeAccnum2HasStatToPrint = {}

    for idx, currListAllICEsIMEsStructure in enumerate(listOfListAllICEsIMEsStructure):

        if len(currListAllICEsIMEsStructure) == 0 :
            continue

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

        #Remark: listOrderedSPs can be empty : The SP TEST1_RS03125 seems to be an isolated relaxase and is therefore not reported as an ICE / IME structure by ICEscreen. Please manually check if it is a false positive or for the potential presence of other undetected SP in its genomic neighborhood.
        genomeAccnumOfCurrListAllICEsIMEsStructureIT = ""
        if len(currListAllICEsIMEsStructure[0].listOrderedSPs) > 0:
            genomeAccnumOfCurrListAllICEsIMEsStructureIT = currListAllICEsIMEsStructure[0].listOrderedSPs[0].genomeAccession
        else:
            if len(currListAllICEsIMEsStructure[0].setSPConjModuleToManuallyCheck) > 0:
                genomeAccnumOfCurrListAllICEsIMEsStructureIT = next(iter(currListAllICEsIMEsStructure[0].setSPConjModuleToManuallyCheck)).genomeAccession
        if len(genomeAccnumOfCurrListAllICEsIMEsStructureIT) == 0:
            raise RuntimeError('Error printOverallStatsToSummaryFile: len(genomeAccnumOfCurrListAllICEsIMEsStructureIT) == 0 ; currListAllICEsIMEsStructure = {} ; len = {} ; elet 0 = {} ; listOrderedSPs = {} of len = {} ; comment = {}'.format(str(currListAllICEsIMEsStructure), str(len(currListAllICEsIMEsStructure)), currListAllICEsIMEsStructure[0].internalIdentifier, str(currListAllICEsIMEsStructure[0].listOrderedSPs), str(len(currListAllICEsIMEsStructure[0].listOrderedSPs)), currListAllICEsIMEsStructure[0].comment))


        if genomeAccnumOfCurrListAllICEsIMEsStructureIT not in dictGenomeAccnum2HasStatToPrint:
            dictGenomeAccnum2HasStatToPrint[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = 1

        for currICEIMEStructure in currListAllICEsIMEsStructure:
            # print stat for each genomic segment and total

            # check that all ICEsIMEsStructure in currListAllICEsIMEsStructure have the same .listOrderedSPs .genomeAccession
            for SPsInlistOrderedSPsIT in currICEIMEStructure.listOrderedSPs:
                if SPsInlistOrderedSPsIT.genomeAccession != genomeAccnumOfCurrListAllICEsIMEsStructureIT:
                    raise RuntimeError('printOverallStatsToSummaryFile error: SPsInlistOrderedSPsIT.genomeAccession {} != genomeAccnumOfCurrListAllICEsIMEsStructureIT {}'.format(str(SPsInlistOrderedSPsIT.genomeAccession), str(genomeAccnumOfCurrListAllICEsIMEsStructureIT)))
            # check for SP in setSPConjModuleToManuallyCheck
            for SPsInsetSPConjModuleToManuallyCheck in currICEIMEStructure.setSPConjModuleToManuallyCheck:
                if SPsInsetSPConjModuleToManuallyCheck.genomeAccession != genomeAccnumOfCurrListAllICEsIMEsStructureIT:
                    raise RuntimeError('printOverallStatsToSummaryFile error: SPsInsetSPConjModuleToManuallyCheck.genomeAccession {} != genomeAccnumOfCurrListAllICEsIMEsStructureIT {}'.format(str(SPsInsetSPConjModuleToManuallyCheck.genomeAccession), str(genomeAccnumOfCurrListAllICEsIMEsStructureIT)))

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
             is01_Unsure) = getcatStructWholeElem(
                 currICEIMEStructure
                 )
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
            if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2totalNumberOfSegmentsWithNoElement:
                dictGenomeAccnum2totalNumberOfSegmentsWithNoElement[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += 1
            else:
                dictGenomeAccnum2totalNumberOfSegmentsWithNoElement[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = 1
        elif (totalElementStructsInSegment) == 1:
            if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2totalNumberOfSegmentsWithOneElement:
                dictGenomeAccnum2totalNumberOfSegmentsWithOneElement[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += 1
            else:
                dictGenomeAccnum2totalNumberOfSegmentsWithOneElement[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = 1
        else:
            if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements:
                dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += 1
            else:
                dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = 1

        if countEMHostOnlySegment > 0:
            if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements:
                dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += 1
            else:
                dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = 1

        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countICETyrTotal:
            dictGenomeAccnum2countICETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countICETyrSegment
        else:
            dictGenomeAccnum2countICETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countICETyrSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countICESerTotal:
            dictGenomeAccnum2countICESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countICESerSegment
        else:
            dictGenomeAccnum2countICESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countICESerSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countICEDDETotal:
            dictGenomeAccnum2countICEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countICEDDESegment
        else:
            dictGenomeAccnum2countICEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countICEDDESegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countConjugationModulesNoIntegraseTotal:
            dictGenomeAccnum2countConjugationModulesNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countConjugationModulesNoIntegraseSegment
        else:
            dictGenomeAccnum2countConjugationModulesNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countConjugationModulesNoIntegraseSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialICETyrTotal:
            dictGenomeAccnum2countPartialICETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialICETyrSegment
        else:
            dictGenomeAccnum2countPartialICETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialICETyrSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialICESerTotal:
            dictGenomeAccnum2countPartialICESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialICESerSegment
        else:
            dictGenomeAccnum2countPartialICESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialICESerSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialICEDDETotal:
            dictGenomeAccnum2countPartialICEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialICEDDESegment
        else:
            dictGenomeAccnum2countPartialICEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialICEDDESegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialICENoIntegraseTotal:
            dictGenomeAccnum2countPartialICENoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialICENoIntegraseSegment
        else:
            dictGenomeAccnum2countPartialICENoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialICENoIntegraseSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIMETyrTotal:
            dictGenomeAccnum2countIMETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIMETyrSegment
        else:
            dictGenomeAccnum2countIMETyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIMETyrSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIMESerTotal:
            dictGenomeAccnum2countIMESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIMESerSegment
        else:
            dictGenomeAccnum2countIMESerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIMESerSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIMEDDETotal:
            dictGenomeAccnum2countIMEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIMEDDESegment
        else:
            dictGenomeAccnum2countIMEDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIMEDDESegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal:
            dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countMobilizableElementsNoIntegraseSegment
        else:
            dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countMobilizableElementsNoIntegraseSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialConjModuleTyrTotal:
            dictGenomeAccnum2countPartialConjModuleTyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialConjModuleTyrSegment
        else:
            dictGenomeAccnum2countPartialConjModuleTyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialConjModuleTyrSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialConjModuleSerTotal:
            dictGenomeAccnum2countPartialConjModuleSerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialConjModuleSerSegment
        else:
            dictGenomeAccnum2countPartialConjModuleSerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialConjModuleSerSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialConjModuleDDETotal:
            dictGenomeAccnum2countPartialConjModuleDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialConjModuleDDESegment
        else:
            dictGenomeAccnum2countPartialConjModuleDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialConjModuleDDESegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal:
            dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countPartialConjModuleNoIntegraseSegment
        else:
            dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countPartialConjModuleNoIntegraseSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countUnsureTotal:
            dictGenomeAccnum2countUnsureTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countUnsureSegment
        else:
            dictGenomeAccnum2countUnsureTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countUnsureSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIntegrasesNotAttributedTyrTotal:
            dictGenomeAccnum2countIntegrasesNotAttributedTyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIntegrasesNotAttributedTyrSegment
        else:
            dictGenomeAccnum2countIntegrasesNotAttributedTyrTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIntegrasesNotAttributedTyrSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIntegrasesNotAttributedSerTotal:
            dictGenomeAccnum2countIntegrasesNotAttributedSerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIntegrasesNotAttributedSerSegment
        else:
            dictGenomeAccnum2countIntegrasesNotAttributedSerTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIntegrasesNotAttributedSerSegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countIntegrasesNotAttributedDDETotal:
            dictGenomeAccnum2countIntegrasesNotAttributedDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countIntegrasesNotAttributedDDESegment
        else:
            dictGenomeAccnum2countIntegrasesNotAttributedDDETotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countIntegrasesNotAttributedDDESegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countEMHostOnlyTotal:
            dictGenomeAccnum2countEMHostOnlyTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countEMHostOnlySegment
        else:
            dictGenomeAccnum2countEMHostOnlyTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countEMHostOnlySegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countEMGuestOnlyTotal:
            dictGenomeAccnum2countEMGuestOnlyTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countEMGuestOnlySegment
        else:
            dictGenomeAccnum2countEMGuestOnlyTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countEMGuestOnlySegment
        if genomeAccnumOfCurrListAllICEsIMEsStructureIT in dictGenomeAccnum2countEMHostAndGuestTotal:
            dictGenomeAccnum2countEMHostAndGuestTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] += countEMHostAndGuestSegment
        else:
            dictGenomeAccnum2countEMHostAndGuestTotal[genomeAccnumOfCurrListAllICEsIMEsStructureIT] = countEMHostAndGuestSegment
    
    strHeaderToPrint = (
        "\n##################### Main parameters of the ICEscreen analysis #####################" +
        "\nMaximum distance in CDS between consecutive SPs within a segments:\t" + str(commonMethods.ConfigParams.maxNumberCDSForSplitSPsByColocalizion) +
        "\nMaximum distance in CDS between consecutive SPs within an IME element:\t" + str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize) +
        "\n\n"
        )
    print(strHeaderToPrint, file=summaryFile)

    # for genomeAccnumIT in dictGenomeAccnum2HasStatToPrint:
    for genomeAccnumIT in listGenomeAccessionFromGenbankFile:

        totalNumberOfCompleteICE = 0
        if genomeAccnumIT in dictGenomeAccnum2countICETyrTotal and genomeAccnumIT in dictGenomeAccnum2countICESerTotal and genomeAccnumIT in dictGenomeAccnum2countICEDDETotal:
            totalNumberOfCompleteICE = dictGenomeAccnum2countICETyrTotal[genomeAccnumIT] + dictGenomeAccnum2countICESerTotal[genomeAccnumIT] + dictGenomeAccnum2countICEDDETotal[genomeAccnumIT]
        totalNumberOfConjugationModules = 0
        if genomeAccnumIT in dictGenomeAccnum2countConjugationModulesNoIntegraseTotal:
            totalNumberOfConjugationModules = dictGenomeAccnum2countConjugationModulesNoIntegraseTotal[genomeAccnumIT]
        totalNumberOfPartialICE = 0
        if genomeAccnumIT in dictGenomeAccnum2countPartialICETyrTotal and genomeAccnumIT in dictGenomeAccnum2countPartialICESerTotal and genomeAccnumIT in dictGenomeAccnum2countPartialICEDDETotal and genomeAccnumIT in dictGenomeAccnum2countPartialICENoIntegraseTotal:
            totalNumberOfPartialICE = dictGenomeAccnum2countPartialICETyrTotal[genomeAccnumIT] + dictGenomeAccnum2countPartialICESerTotal[genomeAccnumIT] + dictGenomeAccnum2countPartialICEDDETotal[genomeAccnumIT] + dictGenomeAccnum2countPartialICENoIntegraseTotal[genomeAccnumIT]
        totalNumberOfCompleteIME = 0
        if genomeAccnumIT in dictGenomeAccnum2countIMETyrTotal and genomeAccnumIT in dictGenomeAccnum2countIMESerTotal and genomeAccnumIT in dictGenomeAccnum2countIMEDDETotal:
            totalNumberOfCompleteIME = dictGenomeAccnum2countIMETyrTotal[genomeAccnumIT] + dictGenomeAccnum2countIMESerTotal[genomeAccnumIT] + dictGenomeAccnum2countIMEDDETotal[genomeAccnumIT]
        totalNumberOfMobilizableElements = 0
        if genomeAccnumIT in dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal:
            totalNumberOfMobilizableElements = dictGenomeAccnum2countMobilizableElementsNoIntegraseTotal[genomeAccnumIT]
        totalNumberOfOtherPartialElements = 0
        if genomeAccnumIT in dictGenomeAccnum2countPartialConjModuleTyrTotal and genomeAccnumIT in dictGenomeAccnum2countPartialConjModuleSerTotal and genomeAccnumIT in dictGenomeAccnum2countPartialConjModuleDDETotal and genomeAccnumIT in dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal:
            totalNumberOfOtherPartialElements = dictGenomeAccnum2countPartialConjModuleTyrTotal[genomeAccnumIT] + dictGenomeAccnum2countPartialConjModuleSerTotal[genomeAccnumIT] + dictGenomeAccnum2countPartialConjModuleDDETotal[genomeAccnumIT] + dictGenomeAccnum2countPartialConjModuleNoIntegraseTotal[genomeAccnumIT]
        totalNumberOfNestedElementsInSegments = 0
        if genomeAccnumIT in dictGenomeAccnum2countEMHostOnlyTotal and genomeAccnumIT in dictGenomeAccnum2countEMGuestOnlyTotal and genomeAccnumIT in dictGenomeAccnum2countEMHostAndGuestTotal:
            totalNumberOfNestedElementsInSegments = dictGenomeAccnum2countEMHostOnlyTotal[genomeAccnumIT] + dictGenomeAccnum2countEMGuestOnlyTotal[genomeAccnumIT] + dictGenomeAccnum2countEMHostAndGuestTotal[genomeAccnumIT]
        totalNumberOfHostElements = 0
        if genomeAccnumIT in dictGenomeAccnum2countEMHostOnlyTotal:
            totalNumberOfHostElements = dictGenomeAccnum2countEMHostOnlyTotal[genomeAccnumIT]
        totalNumberOfGuestElements = 0
        if genomeAccnumIT in dictGenomeAccnum2countEMGuestOnlyTotal:
            totalNumberOfGuestElements = dictGenomeAccnum2countEMGuestOnlyTotal[genomeAccnumIT]
        totalNumberOfHostAndGuestElements = 0
        if genomeAccnumIT in dictGenomeAccnum2countEMHostAndGuestTotal:
            totalNumberOfHostAndGuestElements = dictGenomeAccnum2countEMHostAndGuestTotal[genomeAccnumIT]

        countFromdictGenomeAccnum2totalNumberSP = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberSP:
            countFromdictGenomeAccnum2totalNumberSP = dictGenomeAccnum2totalNumberSP[genomeAccnumIT]
        countFromdictGenomeAccnum2totalNumberIntegrase = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberIntegrase:
            countFromdictGenomeAccnum2totalNumberIntegrase = dictGenomeAccnum2totalNumberIntegrase[genomeAccnumIT]
        countFromdictGenomeAccnum2totalNumberUnassignedIntegrase = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberUnassignedIntegrase:
            countFromdictGenomeAccnum2totalNumberUnassignedIntegrase = dictGenomeAccnum2totalNumberUnassignedIntegrase[genomeAccnumIT]

        countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber = {}
        countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned = {}
        for TypeSPConjModuleIT in icescreen_OO.listTypeSPConjModule:
            if TypeSPConjModuleIT in countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber:
                raise RuntimeError('error : TypeSPConjModuleIT {} in countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber {} '.format(TypeSPConjModuleIT))
            if genomeAccnumIT in dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP:
                TypeSPConjModule_2totalNumberSP = dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP[genomeAccnumIT]
                if TypeSPConjModuleIT in TypeSPConjModule_2totalNumberSP:
                    countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber[TypeSPConjModuleIT] = TypeSPConjModule_2totalNumberSP[TypeSPConjModuleIT]
                else:
                    countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber[TypeSPConjModuleIT] = 0
            else:
                countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber[TypeSPConjModuleIT] = 0
            if TypeSPConjModuleIT in countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned:
                raise RuntimeError('error : TypeSPConjModuleIT {} in countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned {} '.format(TypeSPConjModuleIT))
            if genomeAccnumIT in dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP_unaffected:
                TypeSPConjModule_2totalNumberSP = dictGenomeAccnum_2TypeSPConjModule_2totalNumberSP_unaffected[genomeAccnumIT]
                if TypeSPConjModuleIT in TypeSPConjModule_2totalNumberSP:
                    countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned[TypeSPConjModuleIT] = TypeSPConjModule_2totalNumberSP[TypeSPConjModuleIT]
                else:
                    countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned[TypeSPConjModuleIT] = 0
            else:
                countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned[TypeSPConjModuleIT] = 0
        countFromdictGenomeAccnum2totalNumberFragmentedSP = 0
        if genomeAccnumIT in dictGenomeAccnum2totalsetFragmentedSP:
            countFromdictGenomeAccnum2totalNumberFragmentedSP = len(dictGenomeAccnum2totalsetFragmentedSP[genomeAccnumIT])
        countFromdictGenomeAccnum2totalNumberOfSegmentsWithOneElement = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberOfSegmentsWithOneElement:
            countFromdictGenomeAccnum2totalNumberOfSegmentsWithOneElement = dictGenomeAccnum2totalNumberOfSegmentsWithOneElement[genomeAccnumIT]
        countFromdictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements:
            countFromdictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements = dictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements[genomeAccnumIT]
        countFromdictGenomeAccnum2totalNumberOfSegmentsWithNestedElements = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements:
            countFromdictGenomeAccnum2totalNumberOfSegmentsWithNestedElements = dictGenomeAccnum2totalNumberOfSegmentsWithNestedElements[genomeAccnumIT]
        countFromdictGenomeAccnum2totalNumberOfSegmentsWithNoElement = 0
        if genomeAccnumIT in dictGenomeAccnum2totalNumberOfSegmentsWithNoElement:
            countFromdictGenomeAccnum2totalNumberOfSegmentsWithNoElement = dictGenomeAccnum2totalNumberOfSegmentsWithNoElement[genomeAccnumIT]


        strTotalToPrint = (
            "##################### ICEscreen statistics for genome accession "+genomeAccnumIT+" #####################" +
            "\n\n##### Mobile elements #####" +
            "\n## Number of complete elements" +
            "\nComplete ICE (4 types of SP):\t" + str(totalNumberOfCompleteICE) +
            "\nComplete IME (R+I, R+C+I with distance between consecutive SPs <= "+str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize)+" CDSs):\t" + str(totalNumberOfCompleteIME) +
            "\n## Number of complete modules" +
            "\nConjugation module (R+C+V):\t" + str(totalNumberOfConjugationModules) +
            "\nMobilizable element (R+C with distance between consecutive SPs <= " + str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize) + " CDS):\t" + str(totalNumberOfMobilizableElements) +
            "\n## Number of partial elements" +
            "\nPartial ICE (at least V):\t" + str(totalNumberOfPartialICE) +
            "\nOther partial element (R+C, R+V, V+C with distance between consecutive SPs > " + str(commonMethods.ConfigParams.maxNumberCDSForFilterIMESize) + " CDS):\t" + str(totalNumberOfOtherPartialElements) +
            "\n## Composite elements" +
            "\nTotal nested elements (partial or complete):\t" + str(totalNumberOfNestedElementsInSegments) +
            "\nHost element:\t" + str(totalNumberOfHostElements) +
            "\nGuest element:\t" + str(totalNumberOfGuestElements) +
            "\nElement that are both host and guest:\t" + str(totalNumberOfHostAndGuestElements)
        )
        strTotalToPrint += (
            "\n\n##### Signature proteins (SPs) #####" +
            "\n## Total SPs detected" +
            "\nTotal SPs:\t" + str(countFromdictGenomeAccnum2totalNumberSP) +
            "\nTotal Integrase:\t" + str(countFromdictGenomeAccnum2totalNumberIntegrase)
        )
        for TypeSPConjModuleIT in icescreen_OO.listTypeSPConjModule:
            strTotalToPrint += (
                "\nTotal {}:\t".format(TypeSPConjModuleIT) + str(countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumber[TypeSPConjModuleIT])
            )
        strTotalToPrint += (
            "\nFragmented SPs:\t" + str(countFromdictGenomeAccnum2totalNumberFragmentedSP) +
            "\n## Unassigned SPs"
        )
        countFromdictGenomeAccnum2totalNumberUnassignedSPModuleConj = 0
        for TypeSPConjModuleIT in icescreen_OO.listTypeSPConjModule:
            countFromdictGenomeAccnum2totalNumberUnassignedSPModuleConj += countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned[TypeSPConjModuleIT]
        strTotalToPrint += (
            "\nUnassigned SPs:\t" + str(countFromdictGenomeAccnum2totalNumberUnassignedIntegrase + countFromdictGenomeAccnum2totalNumberUnassignedSPModuleConj)
        )
        strTotalToPrint += (
            "\nUnassigned Integrase:\t" + str(countFromdictGenomeAccnum2totalNumberUnassignedIntegrase)
        )
        for TypeSPConjModuleIT in icescreen_OO.listTypeSPConjModule:
            strTotalToPrint += (
                "\nUnassigned {}:\t".format(TypeSPConjModuleIT) + str(countFromdictGenomeAccnum_2TypeSPConjModule_2totalNumberUnassigned[TypeSPConjModuleIT])
            )
        strTotalToPrint += (
            "\n\n##### Segments #####" +
            "\nNumber of segments with one element:\t" + str(countFromdictGenomeAccnum2totalNumberOfSegmentsWithOneElement) +
            "\nNumber of segments with several elements:\t" + str(countFromdictGenomeAccnum2totalNumberOfSegmentsWithSeveralElements) +
            "\nNumber of segments with nested elements:\t" + str(countFromdictGenomeAccnum2totalNumberOfSegmentsWithNestedElements) +
            "\nNumber of segments with exclusively isolated SPs to manually verify:\t" + str(countFromdictGenomeAccnum2totalNumberOfSegmentsWithNoElement) +
            "\n"
        )
        print(strTotalToPrint, file=summaryFile)
