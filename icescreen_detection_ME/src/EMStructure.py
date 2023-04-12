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
import hit
import icescreen_OO
from EMTypeStructure import assignPutativeTypeStructure
import rulesSeedSPExtension
import re
import rulesAddIntegrases

def returnDashIfEmptyStringOrNone(stringTocheck):
    if len(stringTocheck) == 0 or stringTocheck is None :
        return "-"
    else :
        return stringTocheck

# BasicEMStructure contains attributes that are core to the ICE and
# IME structures and other mobile genetic elements
class BasicEMStructure():
    countInternalIdentifier = 1
    threshold_blast_ali_identity_perc_transfert_ICEFamily_to_structure_module_conj = -1

    def __init__(self, doIncrementCountInternalIdentifier):
        if doIncrementCountInternalIdentifier is True:
            BasicEMStructure.countInternalIdentifier += 1
        # integer, for referencing the structures in the output files
        # and comments. The numbering is not contiguous.
        #self.internalIdentifier = BasicEMStructure.countInternalIdentifier
        self.internalIdentifier = "__IMEICEID_BasicEMStructure__" + str(BasicEMStructure.countInternalIdentifier) + "__"
        # [SP] can be a list if multiple SP of the same type are adjacent
        # in genome
        self.listIntegraseUpstream = []
        # [SP] can be a list if multiple SP of the same type are adjacent
        # in genome
        self.listIntegraseDownstream = []
        # [SP] leave the possibility of a list if changes allow multiple SP
        # of the same type are adjacent in genome
        self.listRelaxase = []
        # [SP] leave the possibility of a list if changes allow multiple SP
        # of the same type are adjacent in genome
        self.listCouplingProtein = []
        # [SP] leave the possibility of a list if changes allow multiple SP
        # of the same type are adjacent in genome
        self.listVirB4 = []
        # list only the different setSPICESuperFamilyFromBlast
        self.setICESuperFamilyFromBlastOfSPConjModule = set()
        # list only the different setSPICEFamilyFromBlast
        self.setICEFamilyFromBlastOfSPConjModule = set()
        # list only the different setSPIMESuperFamilyFromBlast
        self.setIMESuperFamilyFromBlastOfSPConjModule = set()
        # list only the different Relaxase_family_domain_of_most_similar_ref_SPFromBlast
        self.setRelaxase_family_domain_of_most_similar_ref_SPFromBlastOfSPConjModule = set()
        # list only the different Relaxase_family_MOB_of_most_similar_ref_SPFromBlast
        self.setRelaxase_family_MOB_of_most_similar_ref_SPFromBlastOfSPConjModule = set()
        # list only the different Coupling_type_of_most_similar_ref_SPFromBlast
        self.setCoupling_type_of_most_similar_ref_SPFromBlastOfSPConjModule = set()
        # list only the different familyFromHMM
        self.setFamilyFromHMMOfSPConjModule = set()
        # attrbute for ease of navigation in lists
        # [SP] all the SPs of the structure ordered on the genome
        self.listOrderedSPs = []
        # list only the different SPs that can be part of multiple anchors
        self.setSPInConflict = set()

        # list all SPConjModuleLocusTagsToManuallyCheck
        self.setSPConjModuleToManuallyCheck = set()
        self.setIntegraseToManuallyCheck = set()
        # int ; index of the most upstream SP from the conjugation module
        # in the anchor
        self.idxUpstreamConjugationModuleSPInListSPs = -1
        # int ; index of the most downstream SP from the conjugation module
        # in the anchor
        self.idxDownstreamConjugationModuleSPInListSPs = -1
        # str, for exemple if unresolvedNestedOrOverlapStructures, why ?
        # Are there some SP of the conjugaison module that could be
        # attributed to another ICE / IME structure, ect.
        self.comment = ""

    def __hash__(self):
        return hash((self.internalIdentifier))

    def __eq__(self, other):
        if not isinstance(other, BasicEMStructure):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return (self.internalIdentifier) == (other.internalIdentifier)

    def __ne__(self, other):
        if not isinstance(other, BasicEMStructure):
            # don't attempt to compare against unrelated types
            return NotImplemented
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)


    def structureHasThoseIntegraseRegistred(self, listIntegraseToCheckForPresence):
        foundInlistIntegraseUpstream = any(check in listIntegraseToCheckForPresence for check in self.listIntegraseUpstream)
        foundInlistIntegraseDownstream = any(check in listIntegraseToCheckForPresence for check in self.listIntegraseDownstream)
        if foundInlistIntegraseUpstream or foundInlistIntegraseDownstream :
            return True
        else:
            return False



    # findMostUpstreamSpThatIsNotAnIntegrase
    def findMostUpstrSpNotIntegr(self):
        for currSP in self.listOrderedSPs:
            if (currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4"):
                return currSP
        raise RuntimeError("Error in findMostUpstrSpNotIntegr:" + " no non-inegrase SP detected in list" + " {} for structure {}".format(
                                   hit.ListSPs.GetListProtIdsFromListSP(
                                           self.listOrderedSPs),
                                   self.internalIdentifier))

    # findMostDownstreamSpThatIsNotAnIntegrase
    def findMostDownstrSpNotIntegr(self):
        for currSP in reversed(self.listOrderedSPs):
            if (currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4"):
                return currSP
        raise RuntimeError("Error in findMostUpstrSpNotIntegr: no non-inegrase SP detected in list {} for structure {}".format(
                hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs), self.internalIdentifier))

    def clearAllIntegraseUpstream(self):
        for currIntegrase in self.listIntegraseUpstream:
            for currSP in reversed(self.listOrderedSPs):
                if currIntegrase == currSP:
                    self.listOrderedSPs.remove(currSP)
                    break
        self.listIntegraseUpstream.clear()
        self.refreshListIdxOrderedSPs()

    def clearAllIntegraseDownstream(self):
        for currIntegrase in self.listIntegraseDownstream:
            for currSP in reversed(self.listOrderedSPs):
                if currIntegrase == currSP:
                    self.listOrderedSPs.remove(currSP)
                    break
        self.listIntegraseDownstream.clear()
        self.refreshListIdxOrderedSPs()

    def addListIntegraseUpstream(self, listIntegraseUpstreamToAdd):
        for integraseUpstreamToAddIT in listIntegraseUpstreamToAdd :
            if integraseUpstreamToAddIT not in self.listIntegraseUpstream:
                self.listIntegraseUpstream.append(integraseUpstreamToAddIT)
            if integraseUpstreamToAddIT not in self.listOrderedSPs :
                self.listOrderedSPs.append(integraseUpstreamToAddIT)
        self.refreshListIdxOrderedSPs()

    def addListIntegraseDownstream(self, listIntegraseDownstreamToAdd):
        for integraseDownstreamToAddIT in listIntegraseDownstreamToAdd :
            if integraseDownstreamToAddIT not in self.listIntegraseDownstream :
                self.listIntegraseDownstream.append(integraseDownstreamToAddIT)
            if integraseDownstreamToAddIT not in self.listOrderedSPs :
                self.listOrderedSPs.append(integraseDownstreamToAddIT)
        self.refreshListIdxOrderedSPs()


    def transferAllIntegrasesToManuallyCheck(self):
        # for integraseIT in self.listIntegraseUpstream:
        #     if integraseIT.locusTag not in self.setIntegraseToManuallyCheck :
        #         self.setIntegraseToManuallyCheck.add(integraseIT.locusTag)
        # for integraseIT in self.listIntegraseDownstream:
        #     if integraseIT.locusTag not in self.setIntegraseToManuallyCheck :
        #         self.setIntegraseToManuallyCheck.add(integraseIT.locusTag)
        for integraseIT in self.listIntegraseUpstream:
            if integraseIT not in self.setIntegraseToManuallyCheck :
                self.setIntegraseToManuallyCheck.add(integraseIT)
        for integraseIT in self.listIntegraseDownstream:
            if integraseIT not in self.setIntegraseToManuallyCheck :
                self.setIntegraseToManuallyCheck.add(integraseIT)
        self.clearAllIntegraseUpstream()
        self.clearAllIntegraseDownstream()


    def inferICEFamilyFromBlastOfSPConjModule(self):
        self.setICEFamilyFromBlastOfSPConjModule.clear()
        ICEFamilyFromBlastOfSPConjModuleIT = ""
        coherentSPToinferICEFamilyFromBlastOfSPConjModule = True
        for currSP in self.listOrderedSPs:
            #print("** inferICEFamilyFromBlastOfSPConjModule pour structure {} : CDS {}".format(self.internalIdentifier, currSP.locusTag))
            if (currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4"):
                #print("YEP SPType")
                if currSP.SPDetectedByBlast == 1:
                    #print("YEP SPDetectedByBlast : currSP.Blast_ali_identity_perc = {}".format(str(currSP.Blast_ali_identity_perc)))
                    # if len(currSP.setSPICEFamilyFromBlast) != 0:  # add only if not empty or null # OLD METHOD
                    #     self.setICEFamilyFromBlastOfSPConjModule.update(currSP.setSPICEFamilyFromBlast)  # OLD METHOD

                    # threshold_blast_ali_identity_perc_transfert_ICEFamily_to_structure_module_conj 60 as an argument rather than in hard coded
                    # if currSP.Blast_ali_identity_perc >= 60: # original
                    if currSP.Blast_ali_identity_perc >= BasicEMStructure.threshold_blast_ali_identity_perc_transfert_ICEFamily_to_structure_module_conj :
                        if len(currSP.setSPICEFamilyFromBlast) != 0:  # add only if not empty or null
                            #print("inferICEFamilyFromBlastOfSPConjModule pour structure {} : CDS {} with percent id {} has setSPICEFamilyFromBlast = {}".format(self.internalIdentifier, currSP.locusTag, currSP.Blast_ali_identity_perc, currSP.setSPICEFamilyFromBlast))

                            if len(ICEFamilyFromBlastOfSPConjModuleIT) == 0 :
                                if coherentSPToinferICEFamilyFromBlastOfSPConjModule is True :
                                    if len(currSP.setSPICEFamilyFromBlast) == 1:
                                        for currSPSPICEFamilyFromBlastIT in currSP.setSPICEFamilyFromBlast :
                                            ICEFamilyFromBlastOfSPConjModuleIT = currSPSPICEFamilyFromBlastIT
                                    else :
                                        coherentSPToinferICEFamilyFromBlastOfSPConjModule = False   
                            else :
                                if coherentSPToinferICEFamilyFromBlastOfSPConjModule is False :
                                    raise RuntimeError("Error in inferICEFamilyFromBlastOfSPConjModule: len(ICEFamilyFromBlastOfSPConjModuleIT) >= 0 {} and coherentSPToinferICEFamilyFromBlastOfSPConjModule is False".format(ICEFamilyFromBlastOfSPConjModuleIT))

                                if len(currSP.setSPICEFamilyFromBlast) == 1 :
                                    for currSPSPICEFamilyFromBlastIT in currSP.setSPICEFamilyFromBlast :
                                        if currSPSPICEFamilyFromBlastIT == ICEFamilyFromBlastOfSPConjModuleIT :
                                            pass
                                        else :
                                            ICEFamilyFromBlastOfSPConjModuleIT = ""
                                            coherentSPToinferICEFamilyFromBlastOfSPConjModule = False
                                else :
                                    ICEFamilyFromBlastOfSPConjModuleIT = ""
                                    coherentSPToinferICEFamilyFromBlastOfSPConjModule = False

        if len(ICEFamilyFromBlastOfSPConjModuleIT) != 0:  # add only if not empty or null
            self.setICEFamilyFromBlastOfSPConjModule.add(ICEFamilyFromBlastOfSPConjModuleIT)



    # update families in the anchor after modification of member SPs
    def resetSPFamiliesAfterDeletion(self):
        # self.setFamilyFromBlastOfSPConjModule.clear()
        self.setICESuperFamilyFromBlastOfSPConjModule.clear()
        #self.setICEFamilyFromBlastOfSPConjModule.clear()
        self.setIMESuperFamilyFromBlastOfSPConjModule.clear()
        self.setFamilyFromHMMOfSPConjModule.clear()

        for currSP in self.listOrderedSPs:
            if (currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4"):
                if currSP.SPDetectedByBlast == 1:
                    # if len(currSP.setSPFamilyFromBlast) != 0:# add only if not empty or null
                    #    for currFamilyFromBlast in currSP.setSPFamilyFromBlast: #there can be SP with multiple families (i.e. PF01719-PF00910)
                    #        self.setFamilyFromBlastOfSPConjModule.add(currFamilyFromBlast)
                    if len(currSP.setSPICESuperFamilyFromBlast) != 0:  # add only if not empty or null
                        self.setICESuperFamilyFromBlastOfSPConjModule.update(currSP.setSPICESuperFamilyFromBlast)
                    #if len(currSP.setSPICEFamilyFromBlast) != 0:  # add only if not empty or null
                    #    self.setICEFamilyFromBlastOfSPConjModule.update(currSP.setSPICEFamilyFromBlast) # see inferICEFamilyFromBlastOfSPConjModule
                    if len(currSP.setSPIMESuperFamilyFromBlast) != 0:  # add only if not empty or null
                        self.setIMESuperFamilyFromBlastOfSPConjModule.update(currSP.setSPIMESuperFamilyFromBlast)

                    if len(currSP.Relaxase_family_domain_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                        self.setRelaxase_family_domain_of_most_similar_ref_SPFromBlastOfSPConjModule.add(currSP.Relaxase_family_domain_of_most_similar_ref_SPFromBlast)

                    if len(currSP.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                        self.setRelaxase_family_MOB_of_most_similar_ref_SPFromBlastOfSPConjModule.add(currSP.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast)

                    if len(currSP.Coupling_type_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                        self.setCoupling_type_of_most_similar_ref_SPFromBlastOfSPConjModule.add(currSP.Coupling_type_of_most_similar_ref_SPFromBlast)
                        
                # elif currSP.SPDetectedByHMM == 1:
                if currSP.SPDetectedByHMM == 1:

                    if len(currSP.setSPFamilyFromHMM) != 0:  # add only if not empty or null
                        for currFamilyFromHMM in currSP.setSPFamilyFromHMM:
                            # self.setFamilyFromHMMOfSPConjModule.add(currFamilyFromHMM)
                            if currFamilyFromHMM:  # add only if not empty or null
                                # self.setFamilyFromHMMOfSPConjModule.add(currFamilyFromHMM)
                                if currSP.SPType == "Relaxase":
                                    self.setFamilyFromHMMOfSPConjModule.add(f'R:{currFamilyFromHMM}') #OLD
                                if currSP.SPType == "Coupling protein":
                                    self.setFamilyFromHMMOfSPConjModule.add(f'C:{currFamilyFromHMM}')
                                if currSP.SPType == "VirB4":
                                    self.setFamilyFromHMMOfSPConjModule.add(f'V:{currFamilyFromHMM}')

                if currSP.SPDetectedByBlast == 0 and currSP.SPDetectedByHMM == 0:
                    raise RuntimeError("Error in resetSPFamiliesAfterDeletion: currSP.SPDetectedByBlast and currSP.SPDetectedByHMM are both 0 for currSP.locusTag = {}".format(currSP.locusTag))

        self.inferICEFamilyFromBlastOfSPConjModule()

    def addSPToConjugaisonModule(self, SPtoAdd):
        # update the fields:
        # self.listRelaxase
        # self.listCouplingProtein
        # self.listVirB4
        # self.setFamilyFromBlastOfSPConjModule
        # self.setFamilyFromHMMOfSPConjModule
        # self.listOrderedSPs = [] #[SP]

        if (SPtoAdd.SPType == "Relaxase"):
            if SPtoAdd not in self.listRelaxase :
                self.listRelaxase.append(SPtoAdd)
        elif (SPtoAdd.SPType == "Coupling protein"):
            if SPtoAdd not in self.listCouplingProtein :
                self.listCouplingProtein.append(SPtoAdd)
        elif (SPtoAdd.SPType == "VirB4"):
            if SPtoAdd not in self.listVirB4 :
                self.listVirB4.append(SPtoAdd)
        else:
            raise RuntimeError("Error in addSPToConjugaisonModule: unrecognized SPtoAdd.SPType = {} for SPtoAdd = {}".format(SPtoAdd.SPType, SPtoAdd.locusTag))

        # print("adding {}".format(SPtoAdd.locusTag))
        if SPtoAdd not in self.listOrderedSPs :
            self.listOrderedSPs.append(SPtoAdd)
        self.refreshListIdxOrderedSPs()

        if SPtoAdd.SPDetectedByBlast == 1 and not SPtoAdd.pseudo and len(SPtoAdd.listSiblingFragmentedSP) == 0 :
            # if len(SPtoAdd.setSPFamilyFromBlast) != 0:# add only if not empty or null
            #    for currFamilyFromBlast in SPtoAdd.setSPFamilyFromBlast: #there can be SP with multiple families (i.e. PF01719-PF00910)
            #        if currFamilyFromBlast: # add only if not empty or null
            #            self.setFamilyFromBlastOfSPConjModule.add(currFamilyFromBlast)
            if len(SPtoAdd.setSPICESuperFamilyFromBlast) != 0:  # add only if not empty or null
                self.setICESuperFamilyFromBlastOfSPConjModule.update(SPtoAdd.setSPICESuperFamilyFromBlast)
            #if len(SPtoAdd.setSPICEFamilyFromBlast) != 0:  # add only if not empty or null
            #    self.setICEFamilyFromBlastOfSPConjModule.update(SPtoAdd.setSPICEFamilyFromBlast) # see inferICEFamilyFromBlastOfSPConjModule
            if len(SPtoAdd.setSPIMESuperFamilyFromBlast) != 0:  # add only if not empty or null
                self.setIMESuperFamilyFromBlastOfSPConjModule.update(SPtoAdd.setSPIMESuperFamilyFromBlast)
            #print("{} should add Relaxase_family_domain_of_most_similar_ref_SPFromBlast \"{}\"".format(SPtoAdd.locusTag, SPtoAdd.Relaxase_family_domain_of_most_similar_ref_SPFromBlast))
            if len(SPtoAdd.Relaxase_family_domain_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                #print("ADDED !!!")
                self.setRelaxase_family_domain_of_most_similar_ref_SPFromBlastOfSPConjModule.add(SPtoAdd.Relaxase_family_domain_of_most_similar_ref_SPFromBlast)

            if len(SPtoAdd.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                self.setRelaxase_family_MOB_of_most_similar_ref_SPFromBlastOfSPConjModule.add(SPtoAdd.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast)

            if len(SPtoAdd.Coupling_type_of_most_similar_ref_SPFromBlast) != 0:  # add only if not empty or null
                self.setCoupling_type_of_most_similar_ref_SPFromBlastOfSPConjModule.add(SPtoAdd.Coupling_type_of_most_similar_ref_SPFromBlast)

        self.inferICEFamilyFromBlastOfSPConjModule()

        # elif SPtoAdd.SPDetectedByHMM == 1:
        if SPtoAdd.SPDetectedByHMM == 1 and not SPtoAdd.pseudo and len(SPtoAdd.listSiblingFragmentedSP) == 0 :
            if len(SPtoAdd.setSPFamilyFromHMM) != 0:  # add only if not empty or null
                for currFamilyFromHMM in SPtoAdd.setSPFamilyFromHMM:
                    if currFamilyFromHMM:  # add only if not empty or null
                        # self.setFamilyFromHMMOfSPConjModule.add(currFamilyFromHMM)
                        if SPtoAdd.SPType == "Relaxase":
                            self.setFamilyFromHMMOfSPConjModule.add(f'R:{currFamilyFromHMM}')
                        if SPtoAdd.SPType == "Coupling protein":
                            self.setFamilyFromHMMOfSPConjModule.add(f'C:{currFamilyFromHMM}')
                        if SPtoAdd.SPType == "VirB4":
                            self.setFamilyFromHMMOfSPConjModule.add(f'V:{currFamilyFromHMM}')

        if SPtoAdd.SPDetectedByBlast == 0 and SPtoAdd.SPDetectedByHMM == 0:
            raise RuntimeError("Error in addSPToConjugaisonModule: SPtoAdd.SPDetectedByBlast and SPtoAdd.SPDetectedByHMM are both 0 for SPtoAdd.locusTag = {}".format(SPtoAdd.locusTag))


    def refreshListIdxOrderedSPs(self):
        # update the fields:
        # self.idxUpstreamConjugationModuleSPInListSPs = "" # int
        # self.idxDownstreamConjugationModuleSPInListSPs = "" # int
        #print ("refreshListIdxOrderedSPs {}".format(hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)))
        #self.listOrderedSPs.sort(key=lambda x: x.start, reverse=False)
        self.listOrderedSPs.sort(key=lambda x: (x.genomeAccessionRank, x.start), reverse=False)
        

        # for testing
        # strlistOrderedSPsLocusTag = ""
        # for currSP in self.listOrderedSPs:
        #    if strlistOrderedSPsLocusTag:
        #        strlistOrderedSPsLocusTag += ", "
        #    strlistOrderedSPsLocusTag += currSP.locusTag
        # print("refreshListIdxOrderedSPs = {}".format(strlistOrderedSPsLocusTag))

        self.idxUpstreamConjugationModuleSPInListSPs = -1
        self.idxDownstreamConjugationModuleSPInListSPs = -1
        for idx, currSP in enumerate(self.listOrderedSPs):
            if currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4":
                self.idxUpstreamConjugationModuleSPInListSPs = idx
                break
        for idx in range(len(self.listOrderedSPs) - 1, -1, -1):
            currSP = self.listOrderedSPs[idx]
            if currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4":
                self.idxDownstreamConjugationModuleSPInListSPs = idx
                break
        if self.idxUpstreamConjugationModuleSPInListSPs == -1:
            raise RuntimeError("Error in refreshListIdxOrderedSPs: self.idxUpstreamConjugationModuleSPInListSPs == -1 for ICE/IME internalIdentifier = {}".format(str(self.internalIdentifier)))
        if self.idxDownstreamConjugationModuleSPInListSPs == -1:
            raise RuntimeError("Error in refreshListIdxOrderedSPs: self.idxDownstreamConjugationModuleSPInListSPs == -1 for ICE/IME internalIdentifier = {}".format(str(self.internalIdentifier)))

    # check if we need to add upstream subsequent SP of the conjugation module to an anchor. This may result in SP attributed to multiple anchors.
    def addPotentialUpstreamSPInConflict(
            self
            , listSPs
            # , idxCurrentSP
            , listICEsIMEsStructures
            , groupListSPintoICEsIMEsUsingFamilyInfo
            , useFamilyInfoToTryToResolveSPModuleConjConflict
            , useDistanceCDSInfoToTryToResolveSPModuleConjConflict
            , SPsInSameFamilyMergeStructures2SameFamilyMergeStructure
            # , listSameFamilyMergeStructures
            ):
        
        # idxSPsUpstreamToCheck = idxCurrentSP - len(self.listOrderedSPs) # OLD
        idxSPsUpstreamToCheck = self.findMostUpstrSpNotIntegr().idxInListSP - 1

        if idxSPsUpstreamToCheck >= 0:
            listSPsUpstream = listSPs[:(idxSPsUpstreamToCheck + 1)]
            for currentSP in reversed(listSPsUpstream):
                # print("HERE up {}".format(currentSP.locusTag))

                greenLightAddSPConjugaisonModule = False
                # if currentSP is part of SPsInSameFamilyMergeStructures2SameFamilyMergeStructure
                if currentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP]
                    greenLightAddSPConjugaisonModule = self.listSPsIsContainedWithinOtherStructure(currentSP, sameFamilyMergeStructureToCheck)
                    if greenLightAddSPConjugaisonModule:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(self, currentSP, groupListSPintoICEsIMEsUsingFamilyInfo, False, False, False)
                # if self is in SameFamilyMergeStructures but not currentSP
                elif len(self.listOrderedSPs) > 0 and self.listOrderedSPs[0] in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[self.listOrderedSPs[0]]
                    greenLightAddSPConjugaisonModule = self.listSPsIsContainedWithinOtherStructure(currentSP, sameFamilyMergeStructureToCheck)
                    if greenLightAddSPConjugaisonModule:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                            self,
                            currentSP,
                            groupListSPintoICEsIMEsUsingFamilyInfo,
                            False,
                            False,
                            False)
                else:
                    greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                            self,
                            currentSP,
                            groupListSPintoICEsIMEsUsingFamilyInfo,
                            False,
                            False,
                            False)
                if greenLightAddSPConjugaisonModule:

                    # listICEsIMEsStructures can now be the full list
                    # upstreamICEsIMEsStructure = listICEsIMEsStructures[-1] # OLD WAY, listICEsIMEsStructures.append(currentICEsIMEsStructure) CHECK was after addPotentialUpstreamSPInConflict initially
                    upstreamICEsIMEsStructure = ICEsIMEsStructure.getICEsIMEsStructureUpstreamOfICEsIMEsStructure(self, listICEsIMEsStructures, True)

                    if listICEsIMEsStructures:
                        pass
                    else:
                        raise RuntimeError("Error in addPotentialUpstreamSPInConflict: listICEsIMEsStructures is empty for internalIdentifier = {}".format(str(self.internalIdentifier)))

                    if upstreamICEsIMEsStructure is None:
                        break
                    else:
                        resolveStatus = self.tryResolveSPConjugaisonModuleConflictWithUpstreamICEsIMEsStruct(
                                currentSP,
                                upstreamICEsIMEsStructure,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                useFamilyInfoToTryToResolveSPModuleConjConflict,
                                useDistanceCDSInfoToTryToResolveSPModuleConjConflict)
                        if resolveStatus == 0:
                            # the conflict could not be resolved, both ICEsIMEsStructure keep the currentSP
                            # print("HERE 0 the conflict could not be resolved, both ICEsIMEsStructure keep the currentSP = {} ; {} ({}) ; {} ({})" \
                            #      .format(currentSP.locusTag,\
                            #              self.internalIdentifier, \
                            #              hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs), \
                            #              upstreamICEsIMEsStructure.internalIdentifier, \
                            #              hit.ListSPs.GetListProtIdsFromListSP(upstreamICEsIMEsStructure.listOrderedSPs)))
                            self.addSPToConjugaisonModule(currentSP)
                            currentSP.setICEsIMEsStructureInConflict.add(self)
                            currentSP.setICEsIMEsStructureInConflict.add(upstreamICEsIMEsStructure)
                            self.setSPInConflict.add(currentSP)
                            upstreamICEsIMEsStructure.setSPInConflict.add(currentSP)
                        elif resolveStatus == 1:
                            # the currentSP was attributed to this (self) ICEsIMEsStructure
                            self.addSPToConjugaisonModule(currentSP)
                            upstreamICEsIMEsStructure.removeSPConjugaisonModule(currentSP, self, True, False)
                        elif resolveStatus == 2:
                            # the SPInConflict was attributed to upstreamICEsIMEsStructure
                            break
                else:  # not greenLightAddSPConjugaisonModule
                    break

    # return resolveStatus
    # resolveStatus == 0: the conflict could not be resolved, both ICEsIMEsStructure keep the SPInConflict
    # resolveStatus == 1: the SPInConflict was attributed to this (self) ICEsIMEsStructure
    # resolveStatus == 2: the SPInConflict was attributed to upstreamICEsIMEsStructure
    # in case of SP that could be attributed to 2 different anchors, check if the SP could be attributed to one anchor preferably.
    def tryResolveSPConjugaisonModuleConflictWithUpstreamICEsIMEsStruct(
            self,
            SPInConflict,
            upstreamICEsIMEsStructure,
            groupListSPintoICEsIMEsUsingFamilyInfo,
            useFamilyInfoToTryToResolveSPModuleConjConflict,
            useDistanceCDSInfoToTryToResolveSPModuleConjConflict):

        resolveStatusFamilyInfo = 0
        if useFamilyInfoToTryToResolveSPModuleConjConflict == "YES":
            greenLightFamilyInfoTest_self = rulesSeedSPExtension.testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
                    SPInConflict,
                    self,
                    groupListSPintoICEsIMEsUsingFamilyInfo)
            greenLightFamilyInfoTest_upstreamICEsIMEsStructure = rulesSeedSPExtension.testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
                    SPInConflict,
                    upstreamICEsIMEsStructure,
                    groupListSPintoICEsIMEsUsingFamilyInfo)
            if greenLightFamilyInfoTest_self and greenLightFamilyInfoTest_upstreamICEsIMEsStructure:
                resolveStatusFamilyInfo = 0
            elif greenLightFamilyInfoTest_self:
                resolveStatusFamilyInfo = 1
            elif greenLightFamilyInfoTest_upstreamICEsIMEsStructure:
                resolveStatusFamilyInfo = 2
            else:
                resolveStatusFamilyInfo = 0
        elif useFamilyInfoToTryToResolveSPModuleConjConflict == "NO":
            pass
        else:
            raise RuntimeError(
                    "Error in tryResolveSPConjugaisonModuleConflictWithUpstreamICEsIMEsStruct: useFamilyInfoToTryToResolveSPModuleConjConflict is neither YES or NO = {}".format(
                            useFamilyInfoToTryToResolveSPModuleConjConflict))

        resolveStatusDistanceInfo = 0
        if useDistanceCDSInfoToTryToResolveSPModuleConjConflict == "YES":
            distanceCDSsTest_self = abs(SPInConflict.CDSPositionInGenome - self.listOrderedSPs[0].CDSPositionInGenome)
            SPToMesure_upstreamICEsIMEsStructure = upstreamICEsIMEsStructure.listOrderedSPs[upstreamICEsIMEsStructure.idxDownstreamConjugationModuleSPInListSPs]
            distanceCDSsTest_upstreamICEsIMEsStructure = abs(SPInConflict.CDSPositionInGenome - SPToMesure_upstreamICEsIMEsStructure.CDSPositionInGenome)
            if distanceCDSsTest_self < distanceCDSsTest_upstreamICEsIMEsStructure:
                resolveStatusDistanceInfo = 1
            elif distanceCDSsTest_upstreamICEsIMEsStructure < distanceCDSsTest_self:
                resolveStatusDistanceInfo = 2
            else:
                resolveStatusDistanceInfo = 0
        elif useDistanceCDSInfoToTryToResolveSPModuleConjConflict == "NO":
            pass
        else:
            raise RuntimeError(
                    "Error in tryResolveSPConjugaisonModuleConflictWithUpstreamICEsIMEsStruct: useDistanceCDSInfoToTryToResolveSPModuleConjConflict is neither YES or NO = {}".format(
                            useDistanceCDSInfoToTryToResolveSPModuleConjConflict))

        if resolveStatusFamilyInfo == 0 and resolveStatusDistanceInfo == 0:
            return 0
        elif resolveStatusFamilyInfo == 0 and resolveStatusDistanceInfo == 1:
            return 1
        elif resolveStatusFamilyInfo == 0 and resolveStatusDistanceInfo == 2:
            return 2
        elif resolveStatusFamilyInfo == 1 and resolveStatusDistanceInfo == 0:
            return 1
        elif resolveStatusFamilyInfo == 1 and resolveStatusDistanceInfo == 1:
            return 1
        elif resolveStatusFamilyInfo == 1 and resolveStatusDistanceInfo == 2:
            return 0
        elif resolveStatusFamilyInfo == 2 and resolveStatusDistanceInfo == 0:
            return 2
        elif resolveStatusFamilyInfo == 2 and resolveStatusDistanceInfo == 1:
            return 0
        elif resolveStatusFamilyInfo == 2 and resolveStatusDistanceInfo == 2:
            return 2
        else:
            raise RuntimeError(
                    "Error in tryResolveSPConjugaisonModuleConflictWithUpstreamICEsIMEsStruct: unrecognized case for resolveStatusFamilyInfo = {} and resolveStatusDistanceInfo = {}".format(
                            resolveStatusFamilyInfo, resolveStatusDistanceInfo))

    # return boolean True/False
    # check if the SP registred in this EMStructure are contained in another EMStructure, for example FamilyMergeStructure that was built previously
    def listSPsIsContainedWithinOtherStructure(self, currentSP, EMStructureToCompareSent):
        # print("listSPsIsContainedWithinOtherStructure self {}: {}".format(self.internalIdentifier, self.GetObjectAsJson(True, "")))
        # print("listSPsIsContainedWithinOtherStructure EMStructureToCompareSent {}: {}".format(EMStructureToCompareSent.internalIdentifier, EMStructureToCompareSent.GetObjectAsJson(True, "")))

        if currentSP not in EMStructureToCompareSent.listOrderedSPs:
            # The SP we want to add is not in the modele structure
            return False

        if len(self.listOrderedSPs) == 0:
            return True
        else:
            boolToReturn = set(self.listOrderedSPs).issubset(set(EMStructureToCompareSent.listOrderedSPs))
            # print("boolToReturn: {}".format(boolToReturn))
            return boolToReturn



    @staticmethod
    def buildFamily2SetICEsIMEsStructures(listICEsIMEsStructures):
        family2SetICEsIMEsStructures = {}
        for currICEsIMEsStructure in listICEsIMEsStructures:
            for currFamilyFromBlast in currICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule:
                if currFamilyFromBlast in family2SetICEsIMEsStructures:  # key already there
                    currSetICEsIMEsStructures = family2SetICEsIMEsStructures[currFamilyFromBlast]
                    currSetICEsIMEsStructures.add(currICEsIMEsStructure)
                else:  # key not there
                    currSetICEsIMEsStructures = set()
                    currSetICEsIMEsStructures.add(currICEsIMEsStructure)
                    family2SetICEsIMEsStructures[currFamilyFromBlast] = currSetICEsIMEsStructures
            for currFamilyFromBlast in currICEsIMEsStructure.setIMESuperFamilyFromBlastOfSPConjModule:
                if currFamilyFromBlast in family2SetICEsIMEsStructures:  # key already there
                    currSetICEsIMEsStructures = family2SetICEsIMEsStructures[currFamilyFromBlast]
                    currSetICEsIMEsStructures.add(currICEsIMEsStructure)
                else:  # key not there
                    currSetICEsIMEsStructures = set()
                    currSetICEsIMEsStructures.add(currICEsIMEsStructure)
                    family2SetICEsIMEsStructures[currFamilyFromBlast] = currSetICEsIMEsStructures

        return family2SetICEsIMEsStructures


    @staticmethod
    def listStucturesHasAtLeastOneStructureWithoutRegistredIntergase(listEMStructures):
        for EMStructureIT in listEMStructures :
            if len(EMStructureIT.listIntegraseUpstream) == 0 and len(EMStructureIT.listIntegraseDownstream) == 0 :
                return True
        return False



    @staticmethod
    def GetListInternIdFromListEMStructure(listEMStructures):
        strToReturn = ""
        for currentICEsIMEsStructure in listEMStructures:
            if strToReturn:  # not empty
                strToReturn += ", "
            strToReturn += str(currentICEsIMEsStructure.internalIdentifier)
        return strToReturn

    @staticmethod
    def GetListInternIdFromSetEMStructure(setEMStructures):
        strToReturn = ""
        for currentICEsIMEsStructure in sorted(setEMStructures, key=lambda x: x.internalIdentifier, reverse=False):
            if strToReturn:  # not empty
                strToReturn += ", "
            strToReturn += str(currentICEsIMEsStructure.internalIdentifier)
        return strToReturn


    @staticmethod
    def getDigitInEMstructureInternalIdentifier(internalIdentifierSent):
        #RQ : self.internalIdentifier = "__IMEICEID_BasicEMStructure__" + str(BasicEMStructure.countInternalIdentifier) + "__"
        matchInternalIdentifierSent = re.match(
                r'^__IMEICEID_BasicEMStructure__(\d+)__$', internalIdentifierSent)
        if matchInternalIdentifierSent:
            digitInInternalIdentifierSent = int(matchInternalIdentifierSent.group(1))
            return digitInInternalIdentifierSent
        else:
            raise RuntimeError(
                    "Error in getDigitInEMstructureInternalIdentifier: internalIdentifierSent = {} does not match the regex".format(
                                internalIdentifierSent))


    # return the name of the attributes in this object
    @staticmethod
    def GetSummaryObjectHeaderAsTsv():

        printSegmentNumber = True
        segmentNumberHeaderStr = ""
        if printSegmentNumber:
            segmentNumberHeaderStr = "\tSegment_number"

        listStrToReturn = []

        stToReturnPrefix = ""
        stToReturnPrefix += "ICE_IME_id"
        stToReturnPrefix += segmentNumberHeaderStr
        stToReturnPrefix += "\tGenome_accession"

        listStrToReturn.append(stToReturnPrefix)

        stToReturnPostfix = ""
        stToReturnPostfix += "\tICEline_format"
        # stToReturnPostfix += "\tfamily SP conj module Blast"
        stToReturnPostfix += "\tICE_consensus_superfamily_SP_conj_module"
        stToReturnPostfix += "\tICE_consensus_family_SP_conj_module"
        #stToReturnPostfix += "\tIME SuperFamily From Blast Of SP Conj Module"
        stToReturnPostfix += "\tIME_relaxase_family_domains_blast"
        stToReturnPostfix += "\tHMM_family_SP_conj_module"

        stToReturnPostfix += "\tIntegrase_upstream"
        stToReturnPostfix += "\tIntegrase_downstream"
        stToReturnPostfix += "\tRelaxase"
        stToReturnPostfix += "\tCoupling_protein"
        stToReturnPostfix += "\tVirB4"
        stToReturnPostfix += "\tList_SP_ordered_genomic_position"

        stToReturnPostfix += "\tOther_potential_SP_conj_module_need_manual_curation_and_review"
        stToReturnPostfix += "\tOther_potential_integrase_need_manual_curation_and_review"
        # stToReturnPostfix += "\tSPs may belong to another ICE / IME"
        stToReturnPostfix += "\tComments_regarding_structure"
        listStrToReturn.append(stToReturnPostfix)

        return listStrToReturn

    # return the value of the attributes in this object
    def GetSummaryObjectAsTsv(self, 
        # outputPrintCDSNumberInsteadOfProteinIdAndStart,
        segmentNumberSent,
        catStructConjModuleSendIT,
        dictIMEICEID2humanReadableIMEICEIIdentifier):

        printSegmentNumber = True
        segmentNumberStr = ""
        if printSegmentNumber:
            segmentNumberStr = "\t" + str(segmentNumberSent)

        listStrToReturn = []

        stToReturnPrefix = ""
        #stToReturnPrefix += str(self.internalIdentifier)
        stToReturnPrefix += returnDashIfEmptyStringOrNone(dictIMEICEID2humanReadableIMEICEIIdentifier[self.internalIdentifier])
        stToReturnPrefix += returnDashIfEmptyStringOrNone(segmentNumberStr)
        genome_accessionStr = ""
        for spToExtractGenomeAccessionFrom in self.listOrderedSPs:
            if len(genome_accessionStr) == 0:
                genome_accessionStr = spToExtractGenomeAccessionFrom.genomeAccession
            else :
                if genome_accessionStr != spToExtractGenomeAccessionFrom.genomeAccession:
                    raise RuntimeError('GetSummaryObjectAsTsv error: genome_accessionStr {} != spToExtractGenomeAccessionFrom.genomeAccession {}'.format(genome_accessionStr, spToExtractGenomeAccessionFrom.genomeAccession))
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(genome_accessionStr)
        listStrToReturn.append(stToReturnPrefix)

        stToReturnPostfix = ""
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.PrintICElineFormat(self.listOrderedSPs, False, ""))

        if catStructConjModuleSendIT == "ICE" or catStructConjModuleSendIT == "partial ICE VirB4" or catStructConjModuleSendIT == "partial conj module" or catStructConjModuleSendIT == "unsure" :
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setICESuperFamilyFromBlastOfSPConjModule)))
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setICEFamilyFromBlastOfSPConjModule)))
        else :
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone("")
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone("")

        if catStructConjModuleSendIT == "IME" or catStructConjModuleSendIT == "partial conj module"  :
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setRelaxase_family_domain_of_most_similar_ref_SPFromBlastOfSPConjModule)))
        else :
            stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone("")

        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setFamilyFromHMMOfSPConjModule)))
        # if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(sorted(self.listIntegraseUpstream, key=lambda x: (x.genomeAccessionRank, x.start), reverse=False))) #key=lambda x: x.start
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(sorted(self.listIntegraseDownstream, key=lambda x: (x.genomeAccessionRank, x.start), reverse=False))) #key=lambda x: x.start
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(sorted(self.listRelaxase, key=lambda x: (x.genomeAccessionRank, x.start), reverse=False))) #key=lambda x: x.start
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(sorted(self.listCouplingProtein, key=lambda x: (x.genomeAccessionRank, x.start), reverse=False))) #key=lambda x: x.start
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(sorted(self.listVirB4, key=lambda x: (x.genomeAccessionRank, x.start), reverse=False))) #key=lambda x: x.start
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs))
        # elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(sorted(self.listIntegraseUpstream, key=lambda x: x.start, reverse=False)))
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(sorted(self.listIntegraseDownstream, key=lambda x: x.start, reverse=False)))
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(sorted(self.listRelaxase, key=lambda x: x.start, reverse=False)))
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(sorted(self.listCouplingProtein, key=lambda x: x.start, reverse=False)))
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(sorted(self.listVirB4, key=lambda x: x.start, reverse=False)))
        #     stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListCDSNumberFromListSP(self.listOrderedSPs))
        # else:
        #     raise RuntimeError("Error in GetSummaryObjectAsTsv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(outputPrintCDSNumberInsteadOfProteinIdAndStart))

        # stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setSPConjModuleToManuallyCheck)))
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromSetSP(self.setSPConjModuleToManuallyCheck))
        # stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(", ".join(str(i) for i in sorted(self.setIntegraseToManuallyCheck)))
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(hit.ListSPs.GetListProtIdsFromSetSP(self.setIntegraseToManuallyCheck))

        strCommentToPrint = self.comment
        for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
            strCommentToPrint = strCommentToPrint.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)
        stToReturnPostfix += "\t" + returnDashIfEmptyStringOrNone(strCommentToPrint)

        listStrToReturn.append(stToReturnPostfix)

        return listStrToReturn

    # return the value of the attributes in this object in Jason format
    def GetObjectAsJson(self, printCurlyBraces, prefix):
        stToReturn = ""
        if printCurlyBraces:
            stToReturn = "{\n"
        stToReturn += "\t" + prefix + "\"internalIdentifier\": \"" + str(self.internalIdentifier) + "\"\n"
        stToReturn += "\t" + prefix + "\"listIntegraseUpstream\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listIntegraseUpstream) + "\"\n"
        stToReturn += "\t" + prefix + "\"listIntegraseDownstream\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listIntegraseDownstream) + "\"\n"
        stToReturn += "\t" + prefix + "\"listRelaxase\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listRelaxase) + "\"\n"
        stToReturn += "\t" + prefix + "\"listCouplingProtein\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listCouplingProtein) + "\"\n"
        stToReturn += "\t" + prefix + "\"listVirB4\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listVirB4) + "\"\n"
        # stToReturn += "\t"+prefix+"\"setFamilyFromBlastOfSPConjModule\": \"" + ", ".join(str(i) for i in sorted(self.setFamilyFromBlastOfSPConjModule)) + "\"\n"
        stToReturn += "\t" + prefix + "\"setICESuperFamilyFromBlastOfSPConjModule\": \"" + ", ".join(str(i) for i in sorted(self.setICESuperFamilyFromBlastOfSPConjModule)) + "\"\n"
        stToReturn += "\t" + prefix + "\"setICEFamilyFromBlastOfSPConjModule\": \"" + ", ".join(str(i) for i in sorted(self.setICEFamilyFromBlastOfSPConjModule)) + "\"\n"
        stToReturn += "\t" + prefix + "\"setIMESuperFamilyFromBlastOfSPConjModule\": \"" + ", ".join(str(i) for i in sorted(self.setIMESuperFamilyFromBlastOfSPConjModule)) + "\"\n"

        stToReturn += "\t" + prefix + "\"setFamilyFromHMMOfSPConjModule\": \"" + ", ".join(str(i) for i in sorted(self.setFamilyFromHMMOfSPConjModule)) + "\"\n"
        stToReturn += "\t" + prefix + "\"listOrderedSPs\": \"" + hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs) + "\"\n"
        stToReturn += "\t" + prefix + "\"listOrderedSPs as ICElineFormat\": \"" + hit.ListSPs.PrintICElineFormat(self.listOrderedSPs, False, "") + "\"\n"
        # stToReturn += "\t" + prefix + "\"setSPConjModuleToManuallyCheck\": \"" + ", ".join(str(i) for i in sorted(self.setSPConjModuleToManuallyCheck)) + "\"\n"
        stToReturn += "\t" + prefix + "\"setSPConjModuleToManuallyCheck\": \"" + hit.ListSPs.GetListProtIdsFromSetSP(self.setSPConjModuleToManuallyCheck) + "\"\n"
        # stToReturn += "\t" + prefix + "\"setIntegraseToManuallyCheck\": \"" + ", ".join(str(i) for i in sorted(self.setIntegraseToManuallyCheck)) + "\"\n"
        stToReturn += "\t" + prefix + "\"setIntegraseToManuallyCheck\": \"" + hit.ListSPs.GetListProtIdsFromSetSP(self.setIntegraseToManuallyCheck) + "\"\n"
        stToReturn += "\t" + prefix + "\"setSPInConflict\": \"" + hit.ListSPs.GetListProtIdsFromSetSP(self.setSPInConflict) + "\"\n"
        stToReturn += "\t" + prefix + "\"idxUpstreamConjugationModuleSPInListSPs\": \"" + str(self.idxUpstreamConjugationModuleSPInListSPs) + "\"\n"
        stToReturn += "\t" + prefix + "\"idxDownstreamConjugationModuleSPInListSPs\": \"" + str(self.idxDownstreamConjugationModuleSPInListSPs) + "\"\n"
        stToReturn += "\t" + prefix + "\"comment\": \"" + self.comment + "\"\n"
        if printCurlyBraces:
            stToReturn += "}\n"
        return stToReturn


# ICEsIMEsStructure contains attributes that are specific to
# the ICE and IME structures and not mobile genetic elements
class ICEsIMEsStructure(BasicEMStructure):
    # [SP] list of (lowerBoundMergedInternalIdentifierIT,
    # higherBoundMergedInternalIdentifierIT) of structures
    # that have already been merged.
    listLowerBoundHigherBoundStructureMergedInternalIdentifier = []

    def __init__(self, doIncrementCountInternalIdentifier):
        super(ICEsIMEsStructure, self).__init__(
                doIncrementCountInternalIdentifier)
        # see EMTypeStructure.py
        # categoryOfICEsIMEsStructureRegardingWholeElement
        self.catStructWholeElem = ""
        # categoryOfICEsIMEsStructureRegardingConjModule
        self.catStructConjModule = ""
        self.categoryRegardingIntegrase = ""
        self.idxInSeedList = -1  # used for reference for merge events
        # listIdxInSeedListOfDownstreamICEsIMEsStructureMerged
        # empty for not merged
        self.idxListDownstrStructMerged = []
        # deletedAfterMerging_idxInSeedListOfUpstreamICEsIMEsStructureMerged
        # -1 for not merged
        self.delMerging_idxListUpstreamStructure = -1
        self.setHostNestedICEsIMEsStructure = set()  # [ICEsIMEsStructure]
        self.setGuestsNestedICEsIMEsStructure = set()  # [ICEsIMEsStructure]
        self.setOtherICEsIMEsStructureColocalized = set()

    def __hash__(self):
        return super(ICEsIMEsStructure, self).__hash__()

    def __eq__(self, other):
        if not isinstance(other, ICEsIMEsStructure):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return super(ICEsIMEsStructure, self).__eq__(other)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        if not isinstance(other, ICEsIMEsStructure):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return super(ICEsIMEsStructure, self).__ne__(other)

    def isFilterIMESizeOk(self, maxNumberCDSForFilterIMESize):
        if len(self.listOrderedSPs) < 2:
            return False
        distance = abs(self.listOrderedSPs[0].CDSPositionInGenome - self.listOrderedSPs[-1].CDSPositionInGenome)
        if distance <= maxNumberCDSForFilterIMESize:
            return True
        else:
            return False

    def getIterativeMasterMergeStructure(self, listICEsIMEsStructures, levelIteration):
        if levelIteration > 15 :
            raise RuntimeError("Error in getIterativeMasterMergeStructure: levelIteration > 15 for structure {}".format(
                    self.internalIdentifier))
        if self.delMerging_idxListUpstreamStructure == -1 :
            return self
        else :
            masterMergeStructureIT = listICEsIMEsStructures[self.delMerging_idxListUpstreamStructure]
            if self.idxInSeedList in masterMergeStructureIT.idxListDownstrStructMerged :
                if masterMergeStructureIT.delMerging_idxListUpstreamStructure == -1 :
                    return masterMergeStructureIT
                else :
                    levelIteration += 1
                    return masterMergeStructureIT.getIterativeMasterMergeStructure(listICEsIMEsStructures, levelIteration)
            else :
                raise RuntimeError("Error in getIterativeMasterMergeStructure: self.idxInSeedList {} (structure {}) NOT in masterMergeStructureIT.idxListDownstrStructMerged (structure {})".format(
                    self.idxInSeedList,
                    self.internalIdentifier,
                    masterMergeStructureIT.internalIdentifier
                    ))

    # This method merge otherICEsIMEsStructureToMerge in this ICEsIMEsStructure
    def mergeWith(self,
                otherICEsIMEsStructureToMerge,
                groupListSPintoICEsIMEsUsingFamilyInfo,
                locusTagMerge2Comment,
                locusTagIntegrase2Comment,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                allowCheckingForMultipleDistantRelaxase,
                allowCheckingForMultipleDistantCoupling,
                allowCheckingForMultipleDistantVirB4):

        digitInSelfInternalIdentifier = BasicEMStructure.getDigitInEMstructureInternalIdentifier(self.internalIdentifier)
        digitInotherICEsIMEsStructureToMergeInternalIdentifier = BasicEMStructure.getDigitInEMstructureInternalIdentifier(otherICEsIMEsStructureToMerge.internalIdentifier)

        #if self.internalIdentifier >= otherICEsIMEsStructureToMerge.internalIdentifier:
        if digitInSelfInternalIdentifier >= digitInotherICEsIMEsStructureToMergeInternalIdentifier:
            raise RuntimeError("Error in mergeWith: self.internalIdentifier {} >= otherICEsIMEsStructureToMerge.internalIdentifier {}".format(
                    self.internalIdentifier,
                    otherICEsIMEsStructureToMerge.internalIdentifier))

        # print("HERE mergeWith for {} and {}".format(str(self.internalIdentifier), str(otherICEsIMEsStructureToMerge.internalIdentifier)))
        # print("idxInSeedList {} and {}".format(self.idxInSeedList, otherICEsIMEsStructureToMerge.idxInSeedList))
        # print("{}".format(self.GetObjectAsJson(True, "")))
        # print("{}".format(otherICEsIMEsStructureToMerge.GetObjectAsJson(True, "")))

        for (lowerBoundMergedInternalIdentifierIT, higherBoundMergedInternalIdentifierIT) in ICEsIMEsStructure.listLowerBoundHigherBoundStructureMergedInternalIdentifier:
            if self.internalIdentifier == lowerBoundMergedInternalIdentifierIT and otherICEsIMEsStructureToMerge.internalIdentifier == higherBoundMergedInternalIdentifierIT:
                raise RuntimeError("Error in mergeWith: self.internalIdentifier == lowerBoundMergedInternalIdentifierIT {} and otherICEsIMEsStructureToMerge.internalIdentifier == higherBoundMergedInternalIdentifierIT {}".format(
                        self.internalIdentifier,
                        otherICEsIMEsStructureToMerge.internalIdentifier))
            if BasicEMStructure.getDigitInEMstructureInternalIdentifier(lowerBoundMergedInternalIdentifierIT) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(self.internalIdentifier) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(higherBoundMergedInternalIdentifierIT) and BasicEMStructure.getDigitInEMstructureInternalIdentifier(otherICEsIMEsStructureToMerge.internalIdentifier) > BasicEMStructure.getDigitInEMstructureInternalIdentifier(higherBoundMergedInternalIdentifierIT):
                commitITToAdd = "Cannot merge the ICE / IME structure {} ({}) into the ICE / IME structure {} ({}) because this merge will partly overlap with the already merged structures {} and {}. ".format(
                        str(otherICEsIMEsStructureToMerge.internalIdentifier),
                        hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                        str(self.internalIdentifier),
                        hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
                        lowerBoundMergedInternalIdentifierIT,
                        higherBoundMergedInternalIdentifierIT)
                if commitITToAdd not in self.comment:
                    self.comment += commitITToAdd
                if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                    otherICEsIMEsStructureToMerge.comment += commitITToAdd
                for currSP in self.listOrderedSPs:
                    icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                for currSP in otherICEsIMEsStructureToMerge.listOrderedSPs:
                    icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                return False
            if BasicEMStructure.getDigitInEMstructureInternalIdentifier(lowerBoundMergedInternalIdentifierIT) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(otherICEsIMEsStructureToMerge.internalIdentifier) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(higherBoundMergedInternalIdentifierIT) and BasicEMStructure.getDigitInEMstructureInternalIdentifier(self.internalIdentifier) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(lowerBoundMergedInternalIdentifierIT):
                commitITToAdd = "Cannot merge the ICE / IME structure {} ({}) into the ICE / IME structure {} ({}) because this merge will partly overlap with the already merged structures {} and {}. ".format(
                        str(otherICEsIMEsStructureToMerge.internalIdentifier),
                        hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                        str(self.internalIdentifier),
                        hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
                        lowerBoundMergedInternalIdentifierIT,
                        higherBoundMergedInternalIdentifierIT)
                if commitITToAdd not in self.comment:
                    self.comment += commitITToAdd
                if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                    otherICEsIMEsStructureToMerge.comment += commitITToAdd
                for currSP in self.listOrderedSPs:
                    icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                for currSP in otherICEsIMEsStructureToMerge.listOrderedSPs:
                    icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                return False

        if self.delMerging_idxListUpstreamStructure >= 0:  # self has already been merged, do not merge it again
            commitITToAdd = "The ICE / IME structure {} has already been merged with another structure but in theory it could have been also merged with the ICE / IME structure {}. ".format(
                    str(self.internalIdentifier),
                    str(otherICEsIMEsStructureToMerge.internalIdentifier))
            if commitITToAdd not in self.comment:
                self.comment += commitITToAdd
            if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                otherICEsIMEsStructureToMerge.comment += commitITToAdd
            for currSP in self.listOrderedSPs:
                icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
            for currSP in otherICEsIMEsStructureToMerge.listOrderedSPs:
                icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
            return False
        if otherICEsIMEsStructureToMerge.delMerging_idxListUpstreamStructure >= 0:  # otherICEsIMEsStructureToMerge has already been merged, do not merge it again
            commitITToAdd = "The ICE / IME structure {} has already been merged with another structure but in theory it could have been also merged with the ICE / IME structure {}. ".format(
                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                    str(self.internalIdentifier))
            if commitITToAdd not in self.comment:
                self.comment += commitITToAdd
            if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                otherICEsIMEsStructureToMerge.comment += commitITToAdd
            for currSP in self.listOrderedSPs:
                icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
            for currSP in otherICEsIMEsStructureToMerge.listOrderedSPs:
                icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
            return False
            # otherICEsIMEsStructureToMerge.comment += "This ICE / IME structure has been merged already but in theory it could have been also merged with ICE / IME structure " + str(self.internalIdentifier) + ". "
            # self.comment += "In theory this ICE / IME structure could have been merged with ICE / IME structure " + str(otherICEsIMEsStructureToMerge.internalIdentifier) + ", but the latter has been merged to another ICE / IME structure already. "

        commitITToAdd = "The now former ICE / IME structure {} ({}) is being merged into the ICE / IME structure {} ({}). ".format(
                str(otherICEsIMEsStructureToMerge.internalIdentifier),
                hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                str(self.internalIdentifier),
                hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs))
        if self.idxInSeedList >= otherICEsIMEsStructureToMerge.idxInSeedList:
            raise RuntimeError("Error in mergeWith: self {} idxInSeedList {} **NOT** >= otherICEsIMEsStructureToMerge {} idxInSeedList {}".format(
                    str(self.internalIdentifier),
                    str(self.idxInSeedList), str(otherICEsIMEsStructureToMerge.internalIdentifier),
                    str(otherICEsIMEsStructureToMerge.idxInSeedList)))
        if otherICEsIMEsStructureToMerge.idxInSeedList - self.idxInSeedList == 1:
            commitITToAdd += "Please note however that there is no other ICE / IME structure found between the 2 structures that are being merged. It can be due to an integrase from another non ICE / IME mobile element, or a degenerate guest ICE / IME element, or a false positive conjugaison module SP for example. "

        if commitITToAdd not in self.comment:
            self.comment += commitITToAdd
        if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
            otherICEsIMEsStructureToMerge.comment += commitITToAdd
        for currSP in self.listOrderedSPs:
            icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
        for currSP in otherICEsIMEsStructureToMerge.listOrderedSPs:
            icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)

        # deal with listIntegraseUpstream
        if otherICEsIMEsStructureToMerge.listIntegraseUpstream:
            if self.listIntegraseUpstream:
                # check which Integrase is more upstream
                if otherICEsIMEsStructureToMerge.listIntegraseUpstream[0].CDSPositionInGenome < self.listIntegraseUpstream[0].CDSPositionInGenome:  # otherICEsIMEsStructureToMerge is more upstream

                    commitITToAdd = "Both ICE / IME structure merged have registred upstream integrase(s)" \
                        + ", keeping only the most upstream integrase(s) {} and not taking into account {}. ".format(
                                ", ".join(str(i.locusTag) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseUpstream)),
                                ", ".join(str(i.locusTag) for i in sorted(self.listIntegraseUpstream)))
                    if commitITToAdd not in self.comment:
                        self.comment += commitITToAdd
                    if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                        otherICEsIMEsStructureToMerge.comment += commitITToAdd
                    for currSP in self.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    for currSP in otherICEsIMEsStructureToMerge.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    # self.comment += commitITToAdd
                    # self.listIntegraseUpstream = otherICEsIMEsStructureToMerge.listIntegraseUpstream
                    self.clearAllIntegraseUpstream()
                    self.addListIntegraseUpstream(otherICEsIMEsStructureToMerge.listIntegraseUpstream)

                else:  # self is more upstream

                    commitITToAdd = "Both ICE / IME structure merged have registred upstream integrase(s)" \
                        + ", keeping only the most upstream integrase(s) {} and not taking into account {}. ".format(
                                ", ".join(str(i.locusTag) for i in sorted(self.listIntegraseUpstream)),
                                ", ".join(str(i.locusTag) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseUpstream)))
                    if commitITToAdd not in self.comment:
                        self.comment += commitITToAdd
                    if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                        otherICEsIMEsStructureToMerge.comment += commitITToAdd
                    for currSP in self.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    for currSP in otherICEsIMEsStructureToMerge.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    # self.comment += "Both ICE / IME structure merging have registred upstream integrase(s)" \
                    #    + ", keeping only the most upstream integrase(s) for this ICE / IME structure ({}) and discarding the previsouly registred upstream integrase(s) for the ICE / IME structure to merge ({}). " \
                    #    .format( ", ".join(str(i) for i in sorted(self.listIntegraseUpstream)), ", ".join(str(i) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseUpstream)) )
            else:
                # self.comment += "Adding the upstream integrase(s) from the ICE / IME structure to merge ({})" \
                #    .format( ", ".join(str(i) for i in otherICEsIMEsStructureToMerge.listIntegraseUpstream) )
                # self.listIntegraseUpstream = otherICEsIMEsStructureToMerge.listIntegraseUpstream
                self.clearAllIntegraseUpstream()
                self.addListIntegraseUpstream(otherICEsIMEsStructureToMerge.listIntegraseUpstream)

        # deal with listIntegraseDownstream
        if otherICEsIMEsStructureToMerge.listIntegraseDownstream:

            if self.listIntegraseDownstream:
                # check which Integrase is more downstream
                if otherICEsIMEsStructureToMerge.listIntegraseDownstream[0].CDSPositionInGenome > self.listIntegraseDownstream[0].CDSPositionInGenome:  # otherICEsIMEsStructureToMerge is more downstream

                    commitITToAdd = "Both ICE / IME structure merged have registred downstream integrase(s)" \
                        + ", keeping only the most downstream integrase(s) {} and not taking into account {}. ".format(
                                ", ".join(str(i.locusTag) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseDownstream)),
                                ", ".join(str(i.locusTag) for i in sorted(self.listIntegraseDownstream)))
                    if commitITToAdd not in self.comment:
                        self.comment += commitITToAdd
                    if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                        otherICEsIMEsStructureToMerge.comment += commitITToAdd
                    for currSP in self.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    for currSP in otherICEsIMEsStructureToMerge.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)

                    # self.comment += "Both ICE / IME structure merging have registred downstream integrase(s)" \
                    #    + ", keeping only the most downstream integrase(s) of the ICE / IME structure to merge ({}) and discarding the previsouly registred downstream integrase(s) of this ICE / IME structure ({}). " \
                    #    .format( ", ".join(str(i) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseDownstream)), ", ".join(str(i) for i in sorted(self.listIntegraseDownstream)) )
                    # self.listIntegraseDownstream = otherICEsIMEsStructureToMerge.listIntegraseDownstream
                    self.clearAllIntegraseDownstream()
                    self.addListIntegraseDownstream(otherICEsIMEsStructureToMerge.listIntegraseDownstream)

                else:  # self is more downstream

                    commitITToAdd = "Both ICE / IME structure merged have registred downstream integrase(s), keeping only the most downstream integrase(s) {} and not taking into account {}. ".format(
                                ", ".join(str(i.locusTag) for i in sorted(self.listIntegraseDownstream)),
                                ", ".join(str(i.locusTag) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseDownstream)))
                    if commitITToAdd not in self.comment:
                        self.comment += commitITToAdd
                    if commitITToAdd not in otherICEsIMEsStructureToMerge.comment:
                        otherICEsIMEsStructureToMerge.comment += commitITToAdd
                    for currSP in self.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    for currSP in otherICEsIMEsStructureToMerge.listIntegraseUpstream:
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)
                    # self.comment += "Both ICE / IME structure merging have registred downstream integrase(s)" \
                    #    + ", keeping only the most downstream integrase(s) of this ICE / IME structure ({}) and discarding the previsouly registred downstream integrase(s) of the ICE / IME structure to merge ({}). " \
                    #    .format( ", ".join(str(i) for i in sorted(self.listIntegraseDownstream)), ", ".join(str(i) for i in sorted(otherICEsIMEsStructureToMerge.listIntegraseDownstream)) )
            else:
                # self.comment += "Adding the downstream integrase(s) from the ICE / IME structure to merge ({})" \
                #    .format( ", ".join(str(i) for i in otherICEsIMEsStructureToMerge.listIntegraseDownstream) )
                # self.listIntegraseDownstream = otherICEsIMEsStructureToMerge.listIntegraseDownstream
                self.clearAllIntegraseDownstream()
                self.addListIntegraseDownstream(otherICEsIMEsStructureToMerge.listIntegraseDownstream)

        if self.listIntegraseUpstream and self.listIntegraseDownstream:

            (listUpstreamIntChanged, listDownstreamIntChanged) = rulesAddIntegrases.changeObviousIntegraseAttributionToUnsureBecauseOfBothUpstreamAndDownstreamEqualyPossible(
                    self, otherICEsIMEsStructureToMerge, locusTagIntegrase2Comment)

            rulesAddIntegrases.useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase(
                    listUpstreamIntChanged,
                    listDownstreamIntChanged,
                    self,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                    True,
                    locusTagIntegrase2Comment)

            if listUpstreamIntChanged and listDownstreamIntChanged:
                for currSp in listUpstreamIntChanged:
                    # self.setIntegraseToManuallyCheck.add(currSp.locusTag)
                    self.setIntegraseToManuallyCheck.add(currSp)
                for currSp in listDownstreamIntChanged:
                    # self.setIntegraseToManuallyCheck.add(currSp.locusTag)
                    self.setIntegraseToManuallyCheck.add(currSp)
            elif listUpstreamIntChanged:
                for currSp in listUpstreamIntChanged:
                    self.listIntegraseUpstream.append(currSp)
                    self.listOrderedSPs.append(currSp)
                self.refreshListIdxOrderedSPs()
            elif listDownstreamIntChanged:
                for currSp in listDownstreamIntChanged:
                    self.listIntegraseDownstream.append(currSp)
                    self.listOrderedSPs.append(currSp)
                self.refreshListIdxOrderedSPs()

        # deal with SP of conj module, we have to consider that there may be some SP among nested structures to merge that may be in conflict (could be attributed to 2 structures, listICEsIMEsStructureInConflict.isNotEmpty), if so structures with SP in conflict can be merge multiple times and a comment will be made to notify it
        # merge SP of conjugaison module (Relaxase, VirB4, Couplage) from otherICEsIMEsStructureToMerge into this ICEsIMEsStructure, updating listIntegraseDownstream, Relaxase, ect. accordingly
        # deal with listRelaxase
        if otherICEsIMEsStructureToMerge.listRelaxase:
            for currSP in otherICEsIMEsStructureToMerge.listRelaxase:
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(self, currSP, groupListSPintoICEsIMEsUsingFamilyInfo, allowCheckingForMultipleDistantRelaxase, allowCheckingForMultipleDistantCoupling, allowCheckingForMultipleDistantVirB4)
                if greenLightAddSPConjugaisonModule:
                    self.addSPToConjugaisonModule(currSP)
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:  # transfert conflict to merged structure
                        if otherICEsIMEsStructureToMerge in currSP.setICEsIMEsStructureInConflict:
                            currSP.setICEsIMEsStructureInConflict.remove(otherICEsIMEsStructureToMerge)
                            currSP.setICEsIMEsStructureInConflict.add(self)
                            self.setSPInConflict.add(currSP)
                            otherICEsIMEsStructureToMerge.setSPInConflict.remove(currSP)
                        else:
                            raise RuntimeError("Error in mergeWith listRelaxase: currSP.setICEsIMEsStructureInConflict for {} but not otherICEsIMEsStructureToMerge in currSP.setICEsIMEsStructureInConflict ({} ({}) not in [{}]) during merge of structure {} ({}) to structure {} ({})".format(
                                    currSP.locusTag,
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict),
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    str(self.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                    ))
                else:
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:  # SP in conflict were not checked in the method scoreMergeTwoICEsIMEsStructures

                        commitITToAdd = "The relaxase {} could have been attributed to multiple ICE / IME structures ({}) but this \"conflict\" was resolved because of the merge event. ".format(
                                currSP.locusTag,
                                BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict))  # ", ".join(str(i) for i in sorted(currSP.setICEsIMEsStructureInConflict)
                        # self.comment += commitITToAdd
                        if commitITToAdd not in self.comment:
                            self.comment += commitITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)

                    else:
                        raise RuntimeError("Error in mergeWith listRelaxase: not greenLightAddSPConjugaisonModule and not currSP.setICEsIMEsStructureInConflict for merging locus tag {} from structure {} ({}) to structure {} ({})".format(
                                currSP.locusTag,
                                str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                str(self.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                ))
        # deal with listCouplingProtein
        if otherICEsIMEsStructureToMerge.listCouplingProtein:
            for currSP in otherICEsIMEsStructureToMerge.listCouplingProtein:

                # print("{} ({}) and {} : groupListSPintoICEsIMEsUsingFamilyInfo {} ; allowCheckingForMultipleDistantRelaxase {} ; allowCheckingForMultipleDistantCoupling {} ; allowCheckingForMultipleDistantVirB4 {}".format(self.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs), currSP.locusTag, groupListSPintoICEsIMEsUsingFamilyInfo, allowCheckingForMultipleDistantRelaxase, allowCheckingForMultipleDistantCoupling, allowCheckingForMultipleDistantVirB4))

                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(self, currSP, groupListSPintoICEsIMEsUsingFamilyInfo, allowCheckingForMultipleDistantRelaxase, allowCheckingForMultipleDistantCoupling, allowCheckingForMultipleDistantVirB4)
                if greenLightAddSPConjugaisonModule:
                    self.addSPToConjugaisonModule(currSP)
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:  # transfert conflict to merged structure
                        if otherICEsIMEsStructureToMerge in currSP.setICEsIMEsStructureInConflict:
                            currSP.setICEsIMEsStructureInConflict.remove(otherICEsIMEsStructureToMerge)
                            currSP.setICEsIMEsStructureInConflict.add(self)
                            self.setSPInConflict.add(currSP)
                            otherICEsIMEsStructureToMerge.setSPInConflict.remove(currSP)
                        else:
                            raise RuntimeError("Error in mergeWith listCouplingProtein : currSP.setICEsIMEsStructureInConflict for {} but not otherICEsIMEsStructureToMerge.internalIdentifier in currSP.setICEsIMEsStructureInConflict ({} ({}) not in [{}]) during merge of structure {} ({}) to structure {} ({})".format(
                                    currSP.locusTag,
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict),
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    str(self.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                    ))
                else:
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:

                        commitITToAdd = "The coupling protein {} could have been attributed to multiple ICE / IME structures ({}) but this \"conflict\" was resolved because of the merge event. ".format(
                                currSP.locusTag,
                                BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict))  # ", ".join(str(i) for i in sorted(currSP.setICEsIMEsStructureInConflict)
                        # self.comment += commitITToAdd
                        if commitITToAdd not in self.comment:
                            self.comment += commitITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)

                        # self.comment += "The coupling protein {} of the ICE / IME structure to merge could not be added during the merge event to this ICE / IME structure, but this coupling protein could have been attributed to multiple ICE / IME structures anyway: {}. " \
                        #    .format( currSP.locusTag, BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict) )
                    else:
                        # do not work for structure that are made of merged SP (SP that were fusioned because split be an element)
                        raise RuntimeError("Error in mergeWith listCouplingProtein: not greenLightAddSPConjugaisonModule and not currSP.setICEsIMEsStructureInConflict for merging locus tag {} from structure {} ({}) to structure {} ({})".format(
                                currSP.locusTag,
                                str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                str(self.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                ))
        # deal with listVirB4
        if otherICEsIMEsStructureToMerge.listVirB4:
            for currSP in otherICEsIMEsStructureToMerge.listVirB4:
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(self, currSP, groupListSPintoICEsIMEsUsingFamilyInfo, allowCheckingForMultipleDistantRelaxase, allowCheckingForMultipleDistantCoupling, allowCheckingForMultipleDistantVirB4)
                if greenLightAddSPConjugaisonModule:
                    self.addSPToConjugaisonModule(currSP)
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:  # transfert conflict to merged structure
                        if otherICEsIMEsStructureToMerge in currSP.setICEsIMEsStructureInConflict:
                            currSP.setICEsIMEsStructureInConflict.remove(otherICEsIMEsStructureToMerge)
                            currSP.setICEsIMEsStructureInConflict.add(self)
                            self.setSPInConflict.add(currSP)
                            otherICEsIMEsStructureToMerge.setSPInConflict.remove(currSP)
                        else:
                            raise RuntimeError("Error in mergeWith listVirB4: currSP.setICEsIMEsStructureInConflict for {} but not otherICEsIMEsStructureToMerge.internalIdentifier in currSP.setICEsIMEsStructureInConflict ({} ({}) not in [{}]) during merge of structure {} ({}) to structure {} ({})".format(
                                    currSP.locusTag,
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict),
                                    str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                    str(self.internalIdentifier),
                                    hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                    ))

                else:
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:

                        commitITToAdd = "The VirB4 {} could have been attributed to multiple ICE / IME structures ({}) but this \"conflict\" was resolved because of the merge event. ".format(
                                # ", ".join(str(i) for i in sorted(currSP.setICEsIMEsStructureInConflict)
                                currSP.locusTag,
                                BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict))
                        # self.comment += commitITToAdd
                        if commitITToAdd not in self.comment:
                            self.comment += commitITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currSP.locusTag, commitITToAdd, locusTagMerge2Comment)

                        # self.comment += "The VirB4 {} of the ICE / IME structure to merge could not be added during the merge event to this ICE / IME structure, but this VirB4 could have been attributed to multiple ICE / IME structures anyway: {}. " \
                        #    .format( currSP.locusTag, BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict) )
                    else:
                        raise RuntimeError("Error in mergeWith listVirB4: not greenLightAddSPConjugaisonModule and not currSP.setICEsIMEsStructureInConflict for merging locus tag {} from structure {} ({}) to structure {} ({})".format(
                                currSP.locusTag,
                                str(otherICEsIMEsStructureToMerge.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(otherICEsIMEsStructureToMerge.listOrderedSPs),
                                str(self.internalIdentifier),
                                hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)
                                ))

        # integrase will be dealt with in method addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure() launched after merging

        # set flags
        self.idxListDownstrStructMerged.append(otherICEsIMEsStructureToMerge.idxInSeedList)
        otherICEsIMEsStructureToMerge.delMerging_idxListUpstreamStructure = self.idxInSeedList

        # register merge in static variable
        ICEsIMEsStructure.listLowerBoundHigherBoundStructureMergedInternalIdentifier.append(
                (self.internalIdentifier, otherICEsIMEsStructureToMerge.internalIdentifier)
        )

        return True

    def removeSPConjugaisonModule(self,
                                  SPToRemove,
                                  otherICEsIMEsStructureToCheckIfResolvedConflict,
                                  printAComment,
                                  doNotRemoveFromConflictState):
        # update the fields:
        # self.listRelaxase
        # self.listCouplingProtein
        # self.listVirB4
        # self.setFamilyFromBlastOfSPConjModule
        # self.setFamilyFromHMMOfSPConjModule
        # self.listOrderedSPs = [] #[SP]
        if (SPToRemove.SPType == "Relaxase"):
            idxToDel = hit.ListSPs.GetIdxSPInList(SPToRemove, self.listRelaxase)
            if idxToDel == -1:
                raise RuntimeError("Error in removeSPConjugaisonModule Relaxase: idxToDel == -1 for SPToRemove = {}".format(SPToRemove.locusTag))
            del self.listRelaxase[idxToDel]
        elif (SPToRemove.SPType == "Coupling protein"):
            idxToDel = hit.ListSPs.GetIdxSPInList(SPToRemove, self.listCouplingProtein)
            if idxToDel == -1:
                raise RuntimeError("Error in removeSPConjugaisonModule Coupling protein: idxToDel == -1 for SPToRemove = {}".format(SPToRemove.locusTag))
            del self.listCouplingProtein[idxToDel]
        elif (SPToRemove.SPType == "VirB4"):
            idxToDel = hit.ListSPs.GetIdxSPInList(SPToRemove, self.listVirB4)
            if idxToDel == -1:
                raise RuntimeError("Error in removeSPConjugaisonModule VirB4: idxToDel == -1 for SPToRemove = {}".format(SPToRemove.locusTag))
            del self.listVirB4[idxToDel]
        else:
            raise RuntimeError("Error in removeSPConjugaisonModule: unrecognized SPToRemove.SPType = {} for SPToRemove = {}".format(
                    SPToRemove.SPType, SPToRemove.locusTag))

        self.resetSPFamiliesAfterDeletion()

        idxToDel = hit.ListSPs.GetIdxSPInList(SPToRemove, self.listOrderedSPs)
        # print("HERE 1: {}".format(str(idxToDel)))
        if idxToDel == -1:
            raise RuntimeError("Error in removeSPConjugaisonModule listOrderedSPs: idxToDel == -1 for SPToRemove = {}".format(SPToRemove.locusTag))
        del self.listOrderedSPs[idxToDel]
        self.refreshListIdxOrderedSPs()

        # deal with prot in conflict that may be resolved that way
        # remove from ICEsIMEsStruct.setSPInConflict
        if not doNotRemoveFromConflictState:
            if SPToRemove in self.setSPInConflict:
                self.setSPInConflict.remove(SPToRemove)
            if otherICEsIMEsStructureToCheckIfResolvedConflict is not None and SPToRemove in otherICEsIMEsStructureToCheckIfResolvedConflict.setSPInConflict:
                otherICEsIMEsStructureToCheckIfResolvedConflict.setSPInConflict.remove(SPToRemove)
            # remove from SPToRemove.setICEsIMEsStructureInConflict
            if self in SPToRemove.setICEsIMEsStructureInConflict:
                SPToRemove.setICEsIMEsStructureInConflict.remove(self)
            if otherICEsIMEsStructureToCheckIfResolvedConflict is not None and otherICEsIMEsStructureToCheckIfResolvedConflict in SPToRemove.setICEsIMEsStructureInConflict:
                SPToRemove.setICEsIMEsStructureInConflict.remove(otherICEsIMEsStructureToCheckIfResolvedConflict)

        if printAComment:
            commentSupp = ""
            if otherICEsIMEsStructureToCheckIfResolvedConflict is not None:
                commentSupp = ", and the other ICE /IME structure {} potentially in conflict regarding this SP was updated accordingly".format(
                        otherICEsIMEsStructureToCheckIfResolvedConflict.internalIdentifier)
            self.comment += "The SP from the conjugaison module {} was removed from this ICE / IME structure{}. ".format(
                    SPToRemove.locusTag, commentSupp)
            if otherICEsIMEsStructureToCheckIfResolvedConflict is not None:
                otherICEsIMEsStructureToCheckIfResolvedConflict.comment += "The SP from the conjugaison module {} was removed from the other ICE / IME structure {} (potentially in conflict regarding this SP), and this ICE /IME structure was updated accordingly. ".format(
                        self.internalIdentifier, SPToRemove.locusTag)

    # in case a merge events added some SP from conj module and there is another SP from conj module in conflict
    def tryResolveSPsInConflictAfterMergeEvents(self):

        if self.delMerging_idxListUpstreamStructure >= 0:
            return

        # print ("tryResolveSPsInConflictAfterMergeEvents {} ; setSPInConflict = {} "\
        #       .format(self.internalIdentifier, hit.ListSPs.GetListProtIdsFromSetSP(self.setSPInConflict)))
        # ", ".join(str(i) for i in sorted(self.setSPInConflict))

        # check for incoherent data
        for currSPInConflictBis in self.setSPInConflict:
            # currLocusTag
            currLocusTagInConflictHasCoherentData = False
            for currSP in self.listOrderedSPs:
                if currSP == currSPInConflictBis:

                    # print ("if currSP == currSPInConflictBis {} ; GetListInternIdFromSetEMStructure = {} "\
                    #       .format(currSP.locusTag, BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict)))

                    if self in currSP.setICEsIMEsStructureInConflict:
                        if len(currSP.setICEsIMEsStructureInConflict) == 2:
                            currLocusTagInConflictHasCoherentData = True
                            break
                        else:
                            raise RuntimeError(
                                    "Error in tryResolveSPsInConflictAfterMergeEvents: incoherent conflict 0 data for {}: len(currSP.setICEsIMEsStructureInConflict) = {} != 2 ; currSP.setICEsIMEsStructureInConflict = {}".format(
                                            currSPInConflictBis.locusTag,
                                            len(currSP.setICEsIMEsStructureInConflict),
                                            BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict))
                            )
                    else:
                        raise RuntimeError("Error in tryResolveSPsInConflictAfterMergeEvents: incoherent conflict 1 data for {}: {} not in [{}] ".format(
                                currSPInConflictBis.locusTag,
                                self.internalIdentifier,
                                BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict)))
            if not currLocusTagInConflictHasCoherentData:
                raise RuntimeError("Error in tryResolveSPsInConflictAfterMergeEvents: incoherent conflict 2 data for {}: not in [{}] ".format(
                        currSPInConflictBis.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs)))
        for currSP in self.listOrderedSPs:
            if len(currSP.setICEsIMEsStructureInConflict) != 0:
                if currSP in self.setSPInConflict:
                    pass  # ok
                else:
                    raise RuntimeError("Error in tryResolveSPsInConflictAfterMergeEvents: SP {} is in conflict ({}) but is not found in setSPInConflict: [{}]".format(
                            currSP.locusTag,
                            BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict),
                            hit.ListSPs.GetListProtIdsFromSetSP(self.setSPInConflict)))

        # try resolve conflict for SP of conj module if any
        if self.listRelaxase:
            SPToRemoveFromConjModule2ICEsIMEsStructure = {}
            if len(self.listRelaxase) > 2:
                for currSP in self.listRelaxase:
                    if len(currSP.setICEsIMEsStructureInConflict) != 0:
                        # found currSP as SP in conflict
                        for currICEIMEInConflict in currSP.setICEsIMEsStructureInConflict:
                            if currICEIMEInConflict == self:
                                pass
                            else:
                                # found currICEIMEInConflict as otherICEsIMEsStructureBeingInConflict
                                SPToRemoveFromConjModule2ICEsIMEsStructure[currSP] = currICEIMEInConflict
                                break

            if len(SPToRemoveFromConjModule2ICEsIMEsStructure) != 0:
                # only remove if we don't remove all SP from this structure, othewise keep as in conflict
                if len(self.listRelaxase) - len(SPToRemoveFromConjModule2ICEsIMEsStructure) > 0:
                    for SPToRemove, otherICEsIMEsStructure in SPToRemoveFromConjModule2ICEsIMEsStructure.items():
                        self.removeSPConjugaisonModule(SPToRemove, otherICEsIMEsStructure, True, False)
                else:
                    self.comment += "Some Relaxase have been found in conflict but were not removed from this ICE / IME structure following a merge event because it would have remove all the relaxase from this ICE / IME structure. "
        if self.listCouplingProtein:
            SPToRemoveFromConjModule2ICEsIMEsStructure = {}
            if len(self.listCouplingProtein) > 2:
                for currSP in self.listCouplingProtein:
                    if currSP.setICEsIMEsStructureInConflict:
                        # found currSP as SP in conflict
                        for currICEIMEInConflict in currSP.setICEsIMEsStructureInConflict:
                            if currICEIMEInConflict == self:
                                pass
                            else:
                                # found currICEIMEInConflict as otherICEsIMEsStructureBeingInConflict
                                SPToRemoveFromConjModule2ICEsIMEsStructure[currSP] = currICEIMEInConflict
                                break

            if len(SPToRemoveFromConjModule2ICEsIMEsStructure) != 0:
                # only remove if we don't remove all SP from this structure, othewise keep as in conflict
                if len(self.listCouplingProtein) - len(SPToRemoveFromConjModule2ICEsIMEsStructure) > 0:
                    for SPToRemove, otherICEsIMEsStructure in SPToRemoveFromConjModule2ICEsIMEsStructure.items():
                        self.removeSPConjugaisonModule(SPToRemove, otherICEsIMEsStructure, True, False)
                else:
                    self.comment += "Some Coupling protein have been found in conflict but were not removed from this ICE / IME structure following a merge event because it would have remove all the relaxase from this ICE / IME structure. "

        if self.listVirB4:
            SPToRemoveFromConjModule2ICEsIMEsStructure = {}
            if len(self.listVirB4) > 2:
                for currSP in self.listVirB4:
                    if currSP.setICEsIMEsStructureInConflict:
                        # found currSP as SP in conflict
                        for currICEIMEInConflict in currSP.setICEsIMEsStructureInConflict:
                            if currICEIMEInConflict == self:
                                pass
                            else:
                                # found currICEIMEInConflict as otherICEsIMEsStructureBeingInConflict
                                SPToRemoveFromConjModule2ICEsIMEsStructure[currSP] = currICEIMEInConflict
                                break

            if len(SPToRemoveFromConjModule2ICEsIMEsStructure) != 0:
                # only remove if we don't remove all SP from this structure, othewise keep as in conflict
                if len(self.listVirB4) - len(SPToRemoveFromConjModule2ICEsIMEsStructure) > 0:
                    for SPToRemove, otherICEsIMEsStructure in SPToRemoveFromConjModule2ICEsIMEsStructure.items():
                        self.removeSPConjugaisonModule(SPToRemove, otherICEsIMEsStructure, True, False)
                else:
                    self.comment += "Some VirB4 have been found in conflict but were not removed from this ICE / IME structure following a merge event because it would have remove all the relaxase from this ICE / IME structure. "

    # if SP in conflict, int that have a different family than SP proteins
    def finalizeICEIMEStruct(self, listICEsIMEsStructures, maxNumberCDSForFilterIMESize, moveSingleSPToCheck, locusTag2Comment):

        self.fillUpColocalizedOtherICEsIMEsStructures(listICEsIMEsStructures, moveSingleSPToCheck)

        for currSP in self.listOrderedSPs:
            if len(currSP.setICEsIMEsStructureInConflict) != 0:
                # found currSP as SP in conflict

                if currSP not in self.setSPInConflict:
                    raise RuntimeError(
                            "Error in finalizeICEIMEStruct: SP {} is in conflict ({}) but is not found in setSPInConflict: [{}]".format(
                                    currSP.locusTag,
                                    BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict),
                                    hit.ListSPs.GetListProtIdsFromSetSP(self.setSPInConflict)))

                if (currSP.SPType in icescreen_OO.setIntegraseNames):
                    raise RuntimeError("Error in finalizeICEIMEStruct: SP {} is in conflict ({}) but is an integrase: {}".format(
                            currSP.locusTag, BasicEMStructure.GetListInternIdFromSetEMStructure(currSP.setICEsIMEsStructureInConflict), currSP.SPType))

                for currICEIMEInConflict in currSP.setICEsIMEsStructureInConflict:
                    if currICEIMEInConflict == self:
                        pass
                    else:

                        # found currICEIMEInConflict as
                        # otherICEsIMEsStructureBeingInConflict
                        strBothInternalIdentifier = ""
                        if BasicEMStructure.getDigitInEMstructureInternalIdentifier(self.internalIdentifier) < BasicEMStructure.getDigitInEMstructureInternalIdentifier(currICEIMEInConflict.internalIdentifier):
                            strBothInternalIdentifier = str(self.internalIdentifier) + " or " + str(currICEIMEInConflict.internalIdentifier)
                        else:
                            strBothInternalIdentifier = str(currICEIMEInConflict.internalIdentifier) + " or " + str(self.internalIdentifier)
                        commentITToAdd = "The SP {}".format(currSP.locusTag)\
                            + " could be attributed to the ICE / IME"\
                            + " structure {}, please manually check. ".format(
                                    strBothInternalIdentifier)
                        if commentITToAdd not in self.comment:
                            self.comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(
                                currSP.locusTag,
                                commentITToAdd,
                                locusTag2Comment)
                        self.removeSPConjugaisonModule(
                                currSP, None, False, True)
                        # self.setSPConjModuleToManuallyCheck.add(currSP.locusTag)
                        self.setSPConjModuleToManuallyCheck.add(currSP)
                        break

        # add comment in structure if SP has listSiblingFragmentedSP
        for currSP in self.listOrderedSPs:
            if len(currSP.listSiblingFragmentedSP) > 0:
                commentITToAdd = "{} complement each other as fragments to form a single {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(currSP.listSiblingFragmentedSP), currSP.SPType)
                if commentITToAdd not in self.comment:
                    self.comment += commentITToAdd
                # icescreen_OO.addCommentToLocusTag2Comment(
                #         currSP.locusTag,
                #         commentITToAdd,
                #         locusTag2Comment)
        for currSP in self.setSPConjModuleToManuallyCheck:
            if len(currSP.listSiblingFragmentedSP) > 0:
                commentITToAdd = "{} complement each other as fragments to form a single {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(currSP.listSiblingFragmentedSP), currSP.SPType)
                if commentITToAdd not in self.comment:
                    self.comment += commentITToAdd
                # icescreen_OO.addCommentToLocusTag2Comment(
                #         currSP.locusTag,
                #         commentITToAdd,
                #         locusTag2Comment)
        for currSP in self.setIntegraseToManuallyCheck:
            if len(currSP.listSiblingFragmentedSP) > 0:
                commentITToAdd = "{} complement each other as fragments to form a single {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(currSP.listSiblingFragmentedSP), currSP.SPType)
                if commentITToAdd not in self.comment:
                    self.comment += commentITToAdd
                # icescreen_OO.addCommentToLocusTag2Comment(
                #         currSP.locusTag,
                #         commentITToAdd,
                #         locusTag2Comment)

        if moveSingleSPToCheck == "YES":
            totalNumberSPSureAndNotSure = len(self.listOrderedSPs) + len(self.setSPConjModuleToManuallyCheck) + len(self.setIntegraseToManuallyCheck)
            # single SP, move to unsure
            if totalNumberSPSureAndNotSure <= 1 and len(self.listVirB4) == 0:
                currSPLocusTag = ""
                typeCurrSPLocusTagIT = ""
                for currSP in self.listRelaxase:
                    currSPLocusTag = currSP.locusTag
                    typeCurrSPLocusTagIT = "relaxase"
                    # self.setSPConjModuleToManuallyCheck.add(currSP.locusTag)
                    self.setSPConjModuleToManuallyCheck.add(currSP)
                    self.listOrderedSPs.remove(currSP)
                self.listRelaxase.clear()
                for currSP in self.listCouplingProtein:
                    currSPLocusTag = currSP.locusTag
                    typeCurrSPLocusTagIT = "coupling protein"
                    # self.setSPConjModuleToManuallyCheck.add(currSP.locusTag)
                    self.setSPConjModuleToManuallyCheck.add(currSP)
                    self.listOrderedSPs.remove(currSP)
                self.listCouplingProtein.clear()
                for currSP in self.listVirB4:
                    currSPLocusTag = currSP.locusTag
                    typeCurrSPLocusTagIT = "VirB4"
                    # self.setSPConjModuleToManuallyCheck.add(currSP.locusTag)
                    self.setSPConjModuleToManuallyCheck.add(currSP)
                    self.listOrderedSPs.remove(currSP)
                self.listVirB4.clear()
                self.resetSPFamiliesAfterDeletion()

                commentITToAdd = "The SP {}".format(currSPLocusTag)\
                    + " seems to be an isolated "+typeCurrSPLocusTagIT\
                    + " and is therefore not reported as an ICE / IME structure by ICEscreen."\
                    + " Please manually check if it is a false positive or for the potential"\
                    + " presence of other undetected SP in its genomic neighborhood."
                # commentITToAdd = "The SP {} ".format(currSPLocusTag)\
                #     + "is not a VirB4, constitutes a conjugaison module by"\
                #     + " itself, and no integrase has been attributed to the "\
                #     + "element {}".format(self.internalIdentifier)\
                #     + ", please manually check. "
                if commentITToAdd not in self.comment:
                    self.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(currSPLocusTag,
                                                          commentITToAdd,
                                                          locusTag2Comment)

        # last assignPutativeTypeStructure
        assignPutativeTypeStructure(self, maxNumberCDSForFilterIMESize)



    def findIdxInSeedListMostUpstreamICEsIMEsStructureMerged(self):
        idxInSeedListMostUpstreamICEsIMEsStructureMerged = -1
        for idxICEsIMEsStructureOrIntegraseMerged in self.idxListDownstrStructMerged:
            matchObjICEsIMEsStructMerged = re.match(
                    r'^(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            matchObjUpstreamIntegraseMerged = re.match(
                    r'^I<(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            matchObjDownstreamIntegraseMerged = re.match(
                    r'^I>(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            idxICEsIMEsStructureOrIntegraseMerged = -1
            if matchObjICEsIMEsStructMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjICEsIMEsStructMerged.group(1))
            elif matchObjUpstreamIntegraseMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjUpstreamIntegraseMerged.group(1))
            elif matchObjDownstreamIntegraseMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjDownstreamIntegraseMerged.group(1))
            else:
                raise RuntimeError(
                        "Error in findIdxInSeedListMostUpstreamICEsIMEsStructureMerged: unrecognized idxICEsIMEsStructureOrIntegraseMerged = {} for ICEsIMEsStructures {}".format(
                                    str(idxICEsIMEsStructureOrIntegraseMerged),
                                    self.internalIdentifier))

            if idxInSeedListMostUpstreamICEsIMEsStructureMerged == -1:
                idxInSeedListMostUpstreamICEsIMEsStructureMerged = idxICEsIMEsStructureOrIntegraseMerged
            elif idxICEsIMEsStructureOrIntegraseMerged < idxInSeedListMostUpstreamICEsIMEsStructureMerged:
                idxInSeedListMostUpstreamICEsIMEsStructureMerged = idxICEsIMEsStructureOrIntegraseMerged
        if idxInSeedListMostUpstreamICEsIMEsStructureMerged == -1:
            idxInSeedListMostUpstreamICEsIMEsStructureMerged = self.idxInSeedList
        elif idxInSeedListMostUpstreamICEsIMEsStructureMerged > self.idxInSeedList:
            idxInSeedListMostUpstreamICEsIMEsStructureMerged = self.idxInSeedList
            # raise RuntimeError("Error in fillUpColocalizedOtherICEsIMEs
            # Structures: idxInSeedListMostUpstreamICEsIMEsStructureMerged
            # == -1 for ICEsIMEsStructures {}"\
            #                   .format(self.internalIdentifier))
        return idxInSeedListMostUpstreamICEsIMEsStructureMerged

    def findIdxInSeedListMostDownstreamICEsIMEsStructureMerged(self):
        idxInSeedListMostDownstreamICEsIMEsStructureMerged = -1
        for idxICEsIMEsStructureOrIntegraseMerged in self.idxListDownstrStructMerged:
            matchObjICEsIMEsStructMerged = re.match(
                    r'^(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            matchObjUpstreamIntegraseMerged = re.match(
                    r'^I<(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            matchObjDownstreamIntegraseMerged = re.match(
                    r'^I>(\d+)$', str(idxICEsIMEsStructureOrIntegraseMerged))
            idxICEsIMEsStructureOrIntegraseMerged = -1
            if matchObjICEsIMEsStructMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjICEsIMEsStructMerged.group(1))
            elif matchObjUpstreamIntegraseMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjUpstreamIntegraseMerged.group(1))
            elif matchObjDownstreamIntegraseMerged:
                idxICEsIMEsStructureOrIntegraseMerged = int(matchObjDownstreamIntegraseMerged.group(1))
            else:
                raise RuntimeError(
                        "Error in findIdxInSeedListMostDownstreamICEsIMEsStructureMerged: unrecognized idxICEsIMEsStructureOrIntegraseMerged = {} for ICEsIMEsStructures {}".format(
                            str(idxICEsIMEsStructureOrIntegraseMerged),
                            self.internalIdentifier))

            if idxInSeedListMostDownstreamICEsIMEsStructureMerged == -1:
                idxInSeedListMostDownstreamICEsIMEsStructureMerged = idxICEsIMEsStructureOrIntegraseMerged
            elif idxICEsIMEsStructureOrIntegraseMerged > idxInSeedListMostDownstreamICEsIMEsStructureMerged:
                idxInSeedListMostDownstreamICEsIMEsStructureMerged = idxICEsIMEsStructureOrIntegraseMerged

        if idxInSeedListMostDownstreamICEsIMEsStructureMerged < self.idxInSeedList:
            idxInSeedListMostDownstreamICEsIMEsStructureMerged = self.idxInSeedList

        return idxInSeedListMostDownstreamICEsIMEsStructureMerged



    def fillUpColocalizedOtherICEsIMEsStructures(self, listICEsIMEsStructure, moveSingleSPToCheck):

        # try to fill up:
        # ICEsIMEsStructures.setHostNestedICEsIMEsStructure
        # ICEsIMEsStructures.setGuestsNestedICEsIMEsStructure
        # ICEsIMEsStructures.setOtherICEsIMEsStructureColocalized
        if self.delMerging_idxListUpstreamStructure >= 0:
            return  # skip If delMerging_idxListUpstreamStructure

        # moveSingleSPToCheck, do not report
        if moveSingleSPToCheck == "YES":
            totalNumberSPSureAndNotSure = len(self.listOrderedSPs) + len(self.setSPConjModuleToManuallyCheck) + len(self.setIntegraseToManuallyCheck)
            if totalNumberSPSureAndNotSure <= 1 and len(self.listVirB4) == 0:
            # single SP
                return

        self_mostUpstreamStart = -1
        for SPsIT in self.listOrderedSPs:
            if len(SPsIT.setICEsIMEsStructureInConflict) > 0:
                continue
            self_mostUpstreamStart = SPsIT.CDSPositionInGenome
            break
        if self_mostUpstreamStart == -1:
            raise RuntimeError(
                "Error in fillUpColocalizedOtherICEsIMEsStructures: self_mostUpstreamStart == -1 for ICEsIMEsStructure {} ({}) (start {})".format(
                    self.internalIdentifier,
                    hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
                    str(self_mostUpstreamStart)
                    ))
        #self_mostUpstreamStart = self.listOrderedSPs[0].CDSPositionInGenome
        self_mostUpstreamStop = -1
        for SPsIT in reversed(self.listOrderedSPs):
            if len(SPsIT.setICEsIMEsStructureInConflict) > 0:
                continue
            self_mostUpstreamStop = SPsIT.CDSPositionInGenome
            break
        if self_mostUpstreamStop == -1:
            raise RuntimeError(
                "Error in fillUpColocalizedOtherICEsIMEsStructures: self_mostUpstreamStop == -1 for ICEsIMEsStructure {} ({}) (start {} - stop {})".format(
                    self.internalIdentifier,
                    hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
                    str(self_mostUpstreamStart),
                    str(self_mostUpstreamStop)
                    ))
        #self_mostUpstreamStop = self.listOrderedSPs[-1].CDSPositionInGenome

        for ICEsIMEsStructureIT in listICEsIMEsStructure:

            if moveSingleSPToCheck == "YES":
                totalNumberSPSureAndNotSure = len(ICEsIMEsStructureIT.listOrderedSPs) + len(ICEsIMEsStructureIT.setSPConjModuleToManuallyCheck) + len(ICEsIMEsStructureIT.setIntegraseToManuallyCheck)
                if totalNumberSPSureAndNotSure <= 1 and len(ICEsIMEsStructureIT.listVirB4) == 0:
                # single SP
                    continue

            #ICEsIMEsStructureIT_mostUpstreamStart = ICEsIMEsStructureIT.listOrderedSPs[0].CDSPositionInGenome
            ICEsIMEsStructureIT_mostUpstreamStart = -1
            for SPsIT in ICEsIMEsStructureIT.listOrderedSPs:
                if len(SPsIT.setICEsIMEsStructureInConflict) > 0:
                    continue
                ICEsIMEsStructureIT_mostUpstreamStart = SPsIT.CDSPositionInGenome
                break
            if ICEsIMEsStructureIT_mostUpstreamStart == -1:
                raise RuntimeError(
                    "Error in fillUpColocalizedOtherICEsIMEsStructures: ICEsIMEsStructureIT_mostUpstreamStart == -1 for ICEsIMEsStructure {} ({}) (start {})".format(
                        ICEsIMEsStructureIT.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureIT.listOrderedSPs),
                        str(ICEsIMEsStructureIT_mostUpstreamStart)
                        ))
            #ICEsIMEsStructureIT_mostUpstreamStop = ICEsIMEsStructureIT.listOrderedSPs[-1].CDSPositionInGenome
            ICEsIMEsStructureIT_mostUpstreamStop = -1
            for SPsIT in reversed(ICEsIMEsStructureIT.listOrderedSPs):
                if len(SPsIT.setICEsIMEsStructureInConflict) > 0:
                    continue
                ICEsIMEsStructureIT_mostUpstreamStop = SPsIT.CDSPositionInGenome
                break
            if ICEsIMEsStructureIT_mostUpstreamStop == -1:
                raise RuntimeError(
                    "Error in fillUpColocalizedOtherICEsIMEsStructures: ICEsIMEsStructureIT_mostUpstreamStop == -1 for ICEsIMEsStructure {} ({}) (start {} - stop {})".format(
                        ICEsIMEsStructureIT.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureIT.listOrderedSPs),
                        str(ICEsIMEsStructureIT_mostUpstreamStart),
                        str(ICEsIMEsStructureIT_mostUpstreamStop)
                        ))

            # print("ICEsIMEsStructure {} ({}) (start {} - stop {}) and ICEsIMEsStructure {} ({}) (start {} - stop {}).".format(
            #                 self.internalIdentifier,
            #                 hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
            #                 str(self_mostUpstreamStart),
            #                 str(self_mostUpstreamStop),
            #                 ICEsIMEsStructureIT.internalIdentifier,
            #                 hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureIT.listOrderedSPs),
            #                 str(ICEsIMEsStructureIT_mostUpstreamStart),
            #                 str(ICEsIMEsStructureIT_mostUpstreamStop)
            #                 ))

            if ICEsIMEsStructureIT == self:
                # print(" -> ICEsIMEsStructureIT == self")
                continue
            elif ICEsIMEsStructureIT.delMerging_idxListUpstreamStructure >= 0:
                # print(" -> ICEsIMEsStructureIT.delMerging_idxListUpstreamStructure >= 0")
                continue
            elif self_mostUpstreamStop <= ICEsIMEsStructureIT_mostUpstreamStart or ICEsIMEsStructureIT_mostUpstreamStop <= self_mostUpstreamStart :
                # print(" -> setOtherICEsIMEsStructureColocalized")
                self.setOtherICEsIMEsStructureColocalized.add(ICEsIMEsStructureIT)
            elif self_mostUpstreamStart <= ICEsIMEsStructureIT_mostUpstreamStart and self_mostUpstreamStop >= ICEsIMEsStructureIT_mostUpstreamStop :
                # print(" -> setGuestsNestedICEsIMEsStructure")
                self.setGuestsNestedICEsIMEsStructure.add(ICEsIMEsStructureIT)
            elif ICEsIMEsStructureIT_mostUpstreamStart <= self_mostUpstreamStart and ICEsIMEsStructureIT_mostUpstreamStop >= self_mostUpstreamStop:
                # print(" -> setHostNestedICEsIMEsStructure")
                self.setHostNestedICEsIMEsStructure.add(ICEsIMEsStructureIT)
            else :
                raise RuntimeError(
                        "Error in fillUpColocalizedOtherICEsIMEsStructures: unrecognized case for ICEsIMEsStructure {} ({}) (start {} - stop {}) and ICEsIMEsStructure {} ({}) (start {} - stop {}).".format(
                            self.internalIdentifier,
                            hit.ListSPs.GetListProtIdsFromListSP(self.listOrderedSPs),
                            str(self_mostUpstreamStart),
                            str(self_mostUpstreamStop),
                            ICEsIMEsStructureIT.internalIdentifier,
                            hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureIT.listOrderedSPs),
                            str(ICEsIMEsStructureIT_mostUpstreamStart),
                            str(ICEsIMEsStructureIT_mostUpstreamStop)
                            ))




    @staticmethod
    def getICEsIMEsStructureUpstreamOfICEsIMEsStructure(ICEsIMEsStructureSent, listICEsIMEsStructures, returnNoneIfDelMerging):
        # for currICEsIMEsStructureToTest in listICEsIMEsStructures:
        #     if currICEsIMEsStructureToTest == ICEsIMEsStructureSent:
        #         #found our ICEsIMEsStructureSent

        if ICEsIMEsStructureSent not in listICEsIMEsStructures :
            raise RuntimeError(
                        "Error in getICEsIMEsStructureUpstreamOfICEsIMEsStructure: ICEsIMEsStructureSent = {} ({}) not found in list listICEsIMEsStructures {}".format(
                            ICEsIMEsStructureSent.internalIdentifier,
                            hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureSent.listOrderedSPs),
                            BasicEMStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures)))

        idxCurrICEsIMEsStructure = ICEsIMEsStructureSent.findIdxInSeedListMostUpstreamICEsIMEsStructureMerged()

        if idxCurrICEsIMEsStructure == 0:
            return None
        else:
            upstreamOfICEsIMEsStructure = listICEsIMEsStructures[idxCurrICEsIMEsStructure-1]
            if upstreamOfICEsIMEsStructure.delMerging_idxListUpstreamStructure == -1 :
                return upstreamOfICEsIMEsStructure
            else :
                if returnNoneIfDelMerging :
                    return None
                else :
                    return upstreamOfICEsIMEsStructure.getIterativeMasterMergeStructure(listICEsIMEsStructures, 0)
        # raise RuntimeError(
        #                 "Error in getICEsIMEsStructureUpstreamOfICEsIMEsStructure: ICEsIMEsStructureSent = {} not found in list listICEsIMEsStructures {}".format(
        #                     ICEsIMEsStructureSent.internalIdentifier,
        #                     BasicEMStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures)))


    @staticmethod
    def IsThereAnIntegraseBetweenThoseTwoConjModule(ICEsIMEsStructureOne, ICEsIMEsStructureTwo, listSPs, setIntegraseTypeToCheck) :
        idx_StructureOne_MostUpstrSp = ICEsIMEsStructureOne.findMostUpstrSpNotIntegr().idxInListSP
        idx_StructureOne_MostDownstrSp = ICEsIMEsStructureOne.findMostDownstrSpNotIntegr().idxInListSP
        idx_StructureTwo_MostUpstrSp = ICEsIMEsStructureTwo.findMostUpstrSpNotIntegr().idxInListSP
        idx_StructureTwo_MostDownstrSp = ICEsIMEsStructureTwo.findMostDownstrSpNotIntegr().idxInListSP
        if idx_StructureOne_MostUpstrSp > idx_StructureTwo_MostUpstrSp and idx_StructureOne_MostDownstrSp < idx_StructureTwo_MostDownstrSp :
            # ICEsIMEsStructureOne is guest of ICEsIMEsStructureTwo
            return False
        elif idx_StructureOne_MostUpstrSp < idx_StructureTwo_MostUpstrSp and idx_StructureOne_MostDownstrSp > idx_StructureTwo_MostDownstrSp :
            # ICEsIMEsStructureTwo is guest of ICEsIMEsStructureOne
            return False
        elif idx_StructureOne_MostUpstrSp < idx_StructureTwo_MostUpstrSp and idx_StructureOne_MostDownstrSp >= idx_StructureTwo_MostUpstrSp :
            # ICEsIMEsStructureOne is upstream of ICEsIMEsStructureTwo but they share one or more SP (SP in conflict)
            if len(ICEsIMEsStructureOne.findMostDownstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0:
                raise RuntimeError(
                    "Error in IsThereAnIntegraseBetweenThoseTwoConjModule: ICEsIMEsStructureOne is upstream of ICEsIMEsStructureTwo but they share one or more SP (SP in conflict) and len(ICEsIMEsStructureOne.findMostDownstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0 ofICEsIMEsStructureOne = {} (list genes = {}) (idx_StructureOne_MostUpstrSp {} - idx_StructureOne_MostDownstrSp {}) and ICEsIMEsStructureTwo = {} (list genes = {}) (idx_StructureTwo_MostUpstrSp {} - idx_StructureTwo_MostDownstrSp {})".format(
                        ICEsIMEsStructureOne.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOne.listOrderedSPs),
                        str(idx_StructureOne_MostUpstrSp),
                        str(idx_StructureOne_MostDownstrSp),
                        ICEsIMEsStructureTwo.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureTwo.listOrderedSPs),
                        str(idx_StructureTwo_MostUpstrSp),
                        str(idx_StructureTwo_MostDownstrSp)
                        ))
            if len(ICEsIMEsStructureTwo.findMostUpstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0:
                raise RuntimeError(
                    "Error in IsThereAnIntegraseBetweenThoseTwoConjModule: ICEsIMEsStructureOne is upstream of ICEsIMEsStructureTwo but they share one or more SP (SP in conflict) and len(ICEsIMEsStructureTwo.findMostUpstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0 ofICEsIMEsStructureOne = {} (list genes = {}) (idx_StructureOne_MostUpstrSp {} - idx_StructureOne_MostDownstrSp {}) and ICEsIMEsStructureTwo = {} (list genes = {}) (idx_StructureTwo_MostUpstrSp {} - idx_StructureTwo_MostDownstrSp {})".format(
                        ICEsIMEsStructureOne.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOne.listOrderedSPs),
                        str(idx_StructureOne_MostUpstrSp),
                        str(idx_StructureOne_MostDownstrSp),
                        ICEsIMEsStructureTwo.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureTwo.listOrderedSPs),
                        str(idx_StructureTwo_MostUpstrSp),
                        str(idx_StructureTwo_MostDownstrSp)
                        ))
            return False
        elif idx_StructureTwo_MostUpstrSp < idx_StructureOne_MostUpstrSp and idx_StructureTwo_MostDownstrSp >= idx_StructureOne_MostUpstrSp :
            # ICEsIMEsStructureTwo is upstream of ICEsIMEsStructureOne but they share one or more SP (SP in conflict)
            if len(ICEsIMEsStructureTwo.findMostDownstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0:
                raise RuntimeError(
                    "Error in IsThereAnIntegraseBetweenThoseTwoConjModule: ICEsIMEsStructureTwo is upstream of ICEsIMEsStructureOne but they share one or more SP (SP in conflict) and len(ICEsIMEsStructureTwo.findMostDownstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0 ofICEsIMEsStructureOne = {} (list genes = {}) (idx_StructureOne_MostUpstrSp {} - idx_StructureOne_MostDownstrSp {}) and ICEsIMEsStructureTwo = {} (list genes = {}) (idx_StructureTwo_MostUpstrSp {} - idx_StructureTwo_MostDownstrSp {})".format(
                        ICEsIMEsStructureOne.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOne.listOrderedSPs),
                        str(idx_StructureOne_MostUpstrSp),
                        str(idx_StructureOne_MostDownstrSp),
                        ICEsIMEsStructureTwo.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureTwo.listOrderedSPs),
                        str(idx_StructureTwo_MostUpstrSp),
                        str(idx_StructureTwo_MostDownstrSp)
                        ))
            if len(ICEsIMEsStructureOne.findMostUpstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0:
                raise RuntimeError(
                    "Error in IsThereAnIntegraseBetweenThoseTwoConjModule: ICEsIMEsStructureTwo is upstream of ICEsIMEsStructureOne but they share one or more SP (SP in conflict) and len(ICEsIMEsStructureOne.findMostUpstrSpNotIntegr().setICEsIMEsStructureInConflict) == 0 ofICEsIMEsStructureOne = {} (list genes = {}) (idx_StructureOne_MostUpstrSp {} - idx_StructureOne_MostDownstrSp {}) and ICEsIMEsStructureTwo = {} (list genes = {}) (idx_StructureTwo_MostUpstrSp {} - idx_StructureTwo_MostDownstrSp {})".format(
                        ICEsIMEsStructureOne.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOne.listOrderedSPs),
                        str(idx_StructureOne_MostUpstrSp),
                        str(idx_StructureOne_MostDownstrSp),
                        ICEsIMEsStructureTwo.internalIdentifier,
                        hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureTwo.listOrderedSPs),
                        str(idx_StructureTwo_MostUpstrSp),
                        str(idx_StructureTwo_MostDownstrSp)
                        ))
            return False
        elif idx_StructureOne_MostDownstrSp < idx_StructureTwo_MostUpstrSp :
            # ICEsIMEsStructureOne is upstream of ICEsIMEsStructureTwo
            sliceListSPsToConsider = listSPs[idx_StructureOne_MostDownstrSp+1:idx_StructureTwo_MostUpstrSp]
            for SPToCheckIfIntegrase in sliceListSPsToConsider :
                if SPToCheckIfIntegrase.SPType in setIntegraseTypeToCheck:
                    return True
            return False
        elif idx_StructureTwo_MostDownstrSp < idx_StructureOne_MostUpstrSp :
            # ICEsIMEsStructureTwo is upstream of ICEsIMEsStructureOne
            sliceListSPsToConsider = listSPs[idx_StructureTwo_MostDownstrSp+1:idx_StructureOne_MostUpstrSp]
            for SPToCheckIfIntegrase in sliceListSPsToConsider :
                if SPToCheckIfIntegrase.SPType in setIntegraseTypeToCheck:
                    return True
            return False
        else :
            raise RuntimeError(
                "Error in IsThereAnIntegraseBetweenThoseTwoConjModule: unrecognized positioning ofICEsIMEsStructureOne = {} (list genes = {}) (idx_StructureOne_MostUpstrSp {} - idx_StructureOne_MostDownstrSp {}) and ICEsIMEsStructureTwo = {} (list genes = {}) (idx_StructureTwo_MostUpstrSp {} - idx_StructureTwo_MostDownstrSp {})".format(
                    ICEsIMEsStructureOne.internalIdentifier,
                    hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOne.listOrderedSPs),
                    str(idx_StructureOne_MostUpstrSp),
                    str(idx_StructureOne_MostDownstrSp),
                    ICEsIMEsStructureTwo.internalIdentifier,
                    hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureTwo.listOrderedSPs),
                    str(idx_StructureTwo_MostUpstrSp),
                    str(idx_StructureTwo_MostDownstrSp)
                    ))




    @staticmethod
    def getICEsIMEsStructureDownstreamOfICEsIMEsStructure(ICEsIMEsStructureSent, listICEsIMEsStructures, returnNoneIfDelMerging):
        # for idxCurrICEsIMEsStructure, currICEsIMEsStructureToTest in enumerate(listICEsIMEsStructures):
        #     if currICEsIMEsStructureToTest == ICEsIMEsStructureSent:
        #         #found our ICEsIMEsStructureSent
        if ICEsIMEsStructureSent not in listICEsIMEsStructures :
            raise RuntimeError(
                "Error in getICEsIMEsStructureDownstreamOfICEsIMEsStructure: ICEsIMEsStructureSent = {} not found in list listICEsIMEsStructures {}".format(
                    ICEsIMEsStructureSent.internalIdentifier,
                    BasicEMStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures)))

        idxCurrICEsIMEsStructure = ICEsIMEsStructureSent.findIdxInSeedListMostDownstreamICEsIMEsStructureMerged()

        if idxCurrICEsIMEsStructure == len(listICEsIMEsStructures)-1:
            return None
        else:

            downstreamOfICEsIMEsStructureIT = listICEsIMEsStructures[idxCurrICEsIMEsStructure+1]
            if downstreamOfICEsIMEsStructureIT.delMerging_idxListUpstreamStructure == -1 :
                return downstreamOfICEsIMEsStructureIT
            else :
                if returnNoneIfDelMerging :
                    return None
                else :
                    return downstreamOfICEsIMEsStructureIT.getIterativeMasterMergeStructure(listICEsIMEsStructures, 0)
        # raise RuntimeError(
        #                 "Error in getICEsIMEsStructureDownstreamOfICEsIMEsStructure: ICEsIMEsStructureSent = {} not found in list listICEsIMEsStructures {}".format(
        #                     ICEsIMEsStructureSent.internalIdentifier,
        #                     BasicEMStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures)))

    @staticmethod
    def GetSummaryObjectHeaderAsTsv():
        listStrFromSuper = BasicEMStructure.GetSummaryObjectHeaderAsTsv()
        strFromSuperPrefix = listStrFromSuper[0]
        strFromSuperPostfix = listStrFromSuper[1]

        listStrToReturn = []
        stToReturnPrefix = ""
        stToReturnPrefix += strFromSuperPrefix
        stToReturnPrefix += "\tCategory_of_element"
        # stToReturnPrefix += "\tcategory of conj module"
        stToReturnPrefix += "\tCategory_of_integrase"
        stToReturnPrefix += "\tHost_ICE_IME_ids"
        stToReturnPrefix += "\tGuest_ICE_IME_ids"
        stToReturnPrefix += "\tColocalized_ICE_IME_ids"
        stToReturnPrefix += strFromSuperPostfix
        listStrToReturn.append(stToReturnPrefix)

        stToReturnPostfix = ""
        listStrToReturn.append(stToReturnPostfix)

        return listStrToReturn

    def GetSummaryObjectAsTsv(self,
                            #   outputPrintCDSNumberInsteadOfProteinIdAndStart,
                              segmentNumberSent,
                              dictIMEICEID2humanReadableIMEICEIIdentifier):

        listStrFromSuper = super(ICEsIMEsStructure, self).GetSummaryObjectAsTsv(
                # outputPrintCDSNumberInsteadOfProteinIdAndStart,
                segmentNumberSent,
                self.catStructConjModule,
                dictIMEICEID2humanReadableIMEICEIIdentifier)
        strFromSuperPrefix = listStrFromSuper[0]
        strFromSuperPostfix = listStrFromSuper[1]

        listStrToReturn = []
        stToReturnPrefix = ""
        stToReturnPrefix += strFromSuperPrefix
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(self.catStructWholeElem)
        # stToReturnPrefix += "\t"+self.catStructConjModule
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(self.categoryRegardingIntegrase)

        strHostNestedICEsIMEsStructureToPrint = BasicEMStructure.GetListInternIdFromSetEMStructure(self.setHostNestedICEsIMEsStructure)
        for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
            strHostNestedICEsIMEsStructureToPrint = strHostNestedICEsIMEsStructureToPrint.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(strHostNestedICEsIMEsStructureToPrint)

        strGuestsNestedICEsIMEsStructureToPrint = BasicEMStructure.GetListInternIdFromSetEMStructure(
                    self.setGuestsNestedICEsIMEsStructure)
        for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
            strGuestsNestedICEsIMEsStructureToPrint = strGuestsNestedICEsIMEsStructureToPrint.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(strGuestsNestedICEsIMEsStructureToPrint)

        strOtherICEsIMEsStructureColocalizedToPrint = BasicEMStructure.GetListInternIdFromSetEMStructure(
                    self.setOtherICEsIMEsStructureColocalized)
        for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
            strOtherICEsIMEsStructureColocalizedToPrint = strOtherICEsIMEsStructureColocalizedToPrint.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)
        stToReturnPrefix += "\t" + returnDashIfEmptyStringOrNone(strOtherICEsIMEsStructureColocalizedToPrint)

        stToReturnPrefix += strFromSuperPostfix
        listStrToReturn.append(stToReturnPrefix)

        stToReturnPostfix = ""
        listStrToReturn.append(stToReturnPostfix)

        return listStrToReturn

    def GetObjectAsJson(self, printCurlyBraces, prefix):
        strFromSuper = super(ICEsIMEsStructure, self).GetObjectAsJson(
                False, prefix)
        stToReturn = ""
        if printCurlyBraces:
            stToReturn = "{\n"
        stToReturn += strFromSuper
        # stToReturn += "\t"+prefix+"\"type\": \"" + self.type + "\"\n"
        stToReturn += "\t" + prefix + "\"catStructWholeElem\": \"" + self.catStructWholeElem + "\"\n"
        stToReturn += "\t" + prefix + "\"catStructConjModule\": \"" + self.catStructConjModule + "\"\n"
        stToReturn += "\t" + prefix + "\"categoryRegardingIntegrase\": \"" + self.categoryRegardingIntegrase + "\"\n"
        stToReturn += "\t" + prefix + "\"idxInSeedList\": \"" + str(self.idxInSeedList) + "\"\n"
        stToReturn += "\t" + prefix + "\"idxListDownstrStructMerged\": \"" + repr(self.idxListDownstrStructMerged) + "\"\n"
        stToReturn += "\t" + prefix + "\"delMerging_idxListUpstreamStructure\": \"" + str(self.delMerging_idxListUpstreamStructure) + "\"\n"
        stToReturn += "\t" + prefix + "\"setHostNestedICEsIMEsStructure\": [" + BasicEMStructure.GetListInternIdFromSetEMStructure(
                    self.setHostNestedICEsIMEsStructure) + "]\n"
        stToReturn += "\t" + prefix + "\"setGuestsNestedICEsIMEsStructure\": [" + BasicEMStructure.GetListInternIdFromSetEMStructure(
                    self.setGuestsNestedICEsIMEsStructure) + "]\n"
        stToReturn += "\t" + prefix + "\"setOtherICEsIMEsStructureColocalized\": [" + BasicEMStructure.GetListInternIdFromSetEMStructure(
                    self.setOtherICEsIMEsStructureColocalized) + "]\n"
        if printCurlyBraces:
            stToReturn += "}\n"
        return stToReturn
