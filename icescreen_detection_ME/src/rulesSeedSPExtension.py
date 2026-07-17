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
from icescreen_detection_ME.src import icescreen_OO
from icescreen_detection_ME.src import commonMethods
from icescreen_detection_ME.src import hit

# return true if the EM structure already contains sibling of this fragmented SP
def isSPtoTestForAddSiblingFragmentOfASPAlreadyInEMStructureSent(EMStructureSent, SPtoTestForFragmentSibling):
    for SPAlreadyInEMStructureSent in EMStructureSent.listOrderedSPs :
        if SPtoTestForFragmentSibling in SPAlreadyInEMStructureSent.listSiblingFragmentedSP :
            return True
    return False


# return: boolean whether this method give the green light to add a SP of the conjugaison module to this BasicEMStructure
def tryAddingSPToConjugaisonModuleEMStructure(
        EMStructureSent,
        SPtoTestForAdd,
        setAllowCheckingForMultipleDistantSameSPType
        ):
    
    debugtryAddingSPToConjugaisonModuleEMStructure = False

    # consider testing SPtoTestForAdd no matter if it is in conflict or not
    # WARNING there can be some SP that are in conflict in EMStructureSent.listXXX (could be attributed to 2 ICE IME struct), treat them as not there

    if isSPtoTestForAddSiblingFragmentOfASPAlreadyInEMStructureSent(EMStructureSent, SPtoTestForAdd):
        if debugtryAddingSPToConjugaisonModuleEMStructure:
            print("isSPtoTestForAddSiblingFragmentOfASPAlreadyInEMStructureSent {} ({}) for SP {}".format(EMStructureSent.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(EMStructureSent.listOrderedSPs), SPtoTestForAdd.locusTag))
        return True
    

    # exit false if either:
    # (1) the SPtoTestForAdd is an integrase
    if (SPtoTestForAdd.SPType in icescreen_OO.setIntegraseNames):
        if debugtryAddingSPToConjugaisonModuleEMStructure:
            print("tryAddingSPToConjugaisonModuleEMStructure for {}: SPtoTestForAdd.SPType is integrase".format(SPtoTestForAdd.locusTag))
        return False
    # (2) the SPtoTestForAdd is conjugaison module (Relaxase, VirB4, Couplage) && this BasicEMStructure has already at least one of this type of SP registred && the SPtoTestForAdd is not adjacent in the genome to the last SP registred


    elif (SPtoTestForAdd.SPType in icescreen_OO.listTypeSPConjModule): # "Relaxase", "Coupling protein", "VirB4"
        if len(EMStructureSent.TypeSPConjModule2listSP[SPtoTestForAdd.SPType])>0:  # list is not empty
            # get the last SP added to this list that is not in conflict
            for lastSPAdded in reversed(EMStructureSent.TypeSPConjModule2listSP[SPtoTestForAdd.SPType]):
                if len(lastSPAdded.setICEsIMEsStructureInConflict) != 0:  # in conflict, continue loop
                    if debugtryAddingSPToConjugaisonModuleEMStructure:
                        print("in conflict, continue loop")
                    continue
                else:
                    if (abs(lastSPAdded.CDSPositionInGenome - SPtoTestForAdd.CDSPositionInGenome) > icescreen_OO.distanceCDSPositionInGenomeForSPToBeConsideredNotAdjacentInTheGenome[SPtoTestForAdd.SPType]):  # SPs are not adjacent in the genome
                        if debugtryAddingSPToConjugaisonModuleEMStructure:
                            print("tryAddingSPToConjugaisonModuleEMStructure for {}: SPtoTestForAdd.SPType is {} but not adjacent in the genome to already added {}(s). setAllowCheckingForMultipleDistantSameSPType is {}".format(SPtoTestForAdd.locusTag, SPtoTestForAdd.SPType, SPtoTestForAdd.SPType, str(setAllowCheckingForMultipleDistantSameSPType)))
                        if SPtoTestForAdd.SPType not in setAllowCheckingForMultipleDistantSameSPType:
                            return False
                    break
    else:
        raise RuntimeError(
                "Error in tryAddingSPToConjugaisonModuleEMStructure: unrecognized SPtoTestForAdd.SPType = {}".format(SPtoTestForAdd.SPType))

    # (3) test family compatibility depending on the option groupListSPintoICEsIMEsUsingFamilyInfo
    greenLightFamilyInfoTest = testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
            SPtoTestForAdd
            , EMStructureSent
            )
    if not greenLightFamilyInfoTest:
        return False

    return True


def testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
        SPToTest
        , ICEsIMEsStructureToTest
        ):

    debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure = False

    if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
        print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for SPToTest {} and ICEsIMEsStructureToTest {}".format(SPToTest.locusTag, ICEsIMEsStructureToTest.internalIdentifier))

    if len(ICEsIMEsStructureToTest.listOrderedSPs) == 0:  # list is empty
        if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
            print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family, ICEsIMEsStructureToTest.listOrderedSPs is empty".format(SPToTest.locusTag))
        return True
    elif (commonMethods.ConfigParams.groupListSPintoICEsIMEsUsingFamilyInfo == "STRICT"):
        # "STRICT": in order for a SP of the conjugation module to extend a given ICE structure, its ICE super-family (information from Blast matches) needs to match the one of the SP in the given ICE structure. Empty family (no data on family) is considered a super-family of its own and can not be associated with non-empty super-families. Multiple families can only be associated with identical multiple families.
        if SPToTest.SPDetectedByBlast == 1:
            if len(SPToTest.setSPICESuperFamilyFromBlast) == 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) != 0:  # set SP is empty, set ICEsIMEsStructureToTest is not empty
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, set SP is empty, set ICEsIMEsStructureToTest is not empty".format(SPToTest.locusTag))
                return False
            elif len(SPToTest.setSPICESuperFamilyFromBlast) != 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) != 0:  # set SP is not empty, set ICEsIMEsStructureToTest is not empty
                if SPToTest.setSPICESuperFamilyFromBlast == ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule:
                    if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                        print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, set SP is not empty, set ICEsIMEsStructureToTest is not empty, and they overlap".format(SPToTest.locusTag))
                    return True
                else:
                    if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                        print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, set SP is not empty, set ICEsIMEsStructureToTest is not empty, but they are different".format(SPToTest.locusTag))
                    return False
            elif len(SPToTest.setSPICESuperFamilyFromBlast) == 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) == 0:  # set SP is empty, set ICEsIMEsStructureToTest is empty
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, set SP is empty, set ICEsIMEsStructureToTest is empty".format(SPToTest.locusTag))
                return True  # do nothing
            else:  # set SP is not empty, set ICEsIMEsStructureToTest is empty
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, set SP is not empty, set ICEsIMEsStructureToTest is empty".format(SPToTest.locusTag))
                return False
        elif SPToTest.SPDetectedByHMM == 1:
            if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family STRICT, SPDetectedByHMM".format(SPToTest.locusTag))
            return True  # can not use "Best_hmmprofile_JL" information to group protein of a same family
        else:
            raise RuntimeError(
                    "Error in testSPFamilyInfoCompatibilityWithICEsIMEsStructure: SPToTest.SPDetectedByBlast and SPToTest.SPDetectedByHMM are both 0 for SPToTest.locusTag = {}".format(SPToTest.locusTag))
    elif (commonMethods.ConfigParams.groupListSPintoICEsIMEsUsingFamilyInfo == "LOOSE"):
        # "LOOSE": in order for a SP of the conjugation module to extend a given ICE structure, its ICE super-family (information from Blast matches) needs to match the one of the SP in the given ICE structure. Empty ICE super-family (no data) is considered compatible to any ICE super-family. If a SP has multiple ICE super-family, a match of at least one ICE super-family is required, then the multiple ICE super-family of the SP will be registered for the ICE structure.
        if SPToTest.SPDetectedByBlast == 1:
            if len(SPToTest.setSPICESuperFamilyFromBlast) == 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) != 0:  # set SP is empty, set ICEsIMEsStructureToTest is not empty
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family LOOSE, set SP is empty, set ICEsIMEsStructureToTest is not empty".format(SPToTest.locusTag))
                return True  # do nothing
            elif len(SPToTest.setSPICESuperFamilyFromBlast) != 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) != 0:  # set SP is not empty, set ICEsIMEsStructureToTest is not empty
                setInteresctSPFamilies = SPToTest.setSPICESuperFamilyFromBlast.intersection(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule)
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family LOOSE, set SP is not empty, set ICEsIMEsStructureToTest is not empty. Intersect lenght {}: {}".format(
                            SPToTest.locusTag, str(len(setInteresctSPFamilies)), setInteresctSPFamilies))
                if len(setInteresctSPFamilies) != 0:
                    return True
                else:
                    return False
            elif len(SPToTest.setSPICESuperFamilyFromBlast) == 0 and len(ICEsIMEsStructureToTest.setICESuperFamilyFromBlastOfSPConjModule) == 0:  # set SP is empty, set ICEsIMEsStructureToTest is empty
                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family LOOSE, set SP is empty, set ICEsIMEsStructureToTest is empty".format(SPToTest.locusTag))
                return True  # do nothing
            else:  # set SP is not empty, set ICEsIMEsStructureToTest is empty

                if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                    print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family LOOSE, set SP is not empty, set ICEsIMEsStructureToTest is empty".format(SPToTest.locusTag))
                return True  # do nothing
        elif SPToTest.SPDetectedByHMM == 1:
            if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
                print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: test family LOOSE, SPDetectedByHMM".format(SPToTest.locusTag))
            return True  # can not use "Best_hmmprofile_JL" information to group protein of a same family
        else:
            raise RuntimeError(
                    "Error in testSPFamilyInfoCompatibilityWithICEsIMEsStructure: SPToTest.SPDetectedByBlast and SPToTest.SPDetectedByHMM are both 0 for SPToTest.locusTag = {}".format(SPToTest.locusTag))
    elif (commonMethods.ConfigParams.groupListSPintoICEsIMEsUsingFamilyInfo == "NO"):  # "NO" means that the SP of the conjugaison module can be be affilitated with any family

        if debugtestSPFamilyInfoCompatibilityWithICEsIMEsStructure:
            print("testSPFamilyInfoCompatibilityWithICEsIMEsStructure for {}: the SP of the conjugaison module can be be affilitated with any family")
        return True  # do nothing
    else:
        raise RuntimeError(
                "Error in testSPFamilyInfoCompatibilityWithICEsIMEsStructure: unrecognized icescreen_OO.groupListSPintoICEsIMEsUsingFamilyInfo = {}".format(
                        str(commonMethods.ConfigParams.groupListSPintoICEsIMEsUsingFamilyInfo)))
