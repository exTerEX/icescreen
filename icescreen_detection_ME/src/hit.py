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
import icescreen_OO
import rulesSeedSPExtension
import rulesAddIntegrases


# The CDS object class mimics the attributes of the biological object CDS
class CDS():
    def __init__(self):
        self.locusTag = ""  # col "CDS_locus_tag" or if absent a combination of "Genome_accession"-"CDS_protein_id"-"CDS_start" ("Genome_accession" is mandatory only if multigenbank, otherwise ignored)
        self.proteinId = ""  # col CDS
        self.genomeAccession = ""  # col "Genome_accession"
        self.genomeAccessionRank = ""  # col "Genome_accession"
        self.start = ""  # col CDS_start
        self.stop = ""  # col CDS_end
        self.strand = ""  # col CDS_strand ; Values: +, -
        self.CDSPositionInGenome = ""  # col CDS_num
        self.idxInListSP = -1
        self.pseudo = False

    def __hash__(self):
        return hash((self.locusTag))

    def __eq__(self, other):
        if not isinstance(other, CDS):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return (self.locusTag) == (other.locusTag)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        if not isinstance(other, CDS):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return not(self == other)

    

# The SP object class extends CDS and add extra attributes to mimic the biological object signature protein
class SP(CDS):
    def __init__(self):
        super(SP, self).__init__()
        self.SPType = ""  # col CDS_Protein_Type ; Values: Coupling, Relaxase, Virb4, icescreen_OO.setIntegrase ; modify in tryAddingSPToConjugaisonModuleEMStructure() if make changes
        self.SPDetectedByBlast = ""  # col Is_hit_blast ; Values: 1, 0
        self.SPDetectedByHMM = ""  # col Is_hit_HMM_JL or Is_hit_HMM_CC depending ; Values: 1, 0
        # self.SPFamilyFromBlast = ""
        # self.setSPFamilyFromBlast = set()
        self.setSPICESuperFamilyFromBlast = set()  # col ICE_superfamily_of_most_similar_ref_SP
        self.setSPICEFamilyFromBlast = set()  # col ICE_family_of_most_similar_ref_SP
        self.setSPIMESuperFamilyFromBlast = set()  # col IME_superfamily_of_most_similar_ref_SP
        self.Relaxase_family_domain_of_most_similar_ref_SPFromBlast = ""
        self.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast = ""
        self.Coupling_type_of_most_similar_ref_SPFromBlast = ""
        self.False_positivesFromBlast = ""
        self.blast_validation = ""
        self.Use_annotationFromBlast = ""
        self.setSPFamilyFromHMM = set()  # col Best_hmmprofile_JL or Best_hmmprofile_CC depending
        self.setICEsIMEsStructureInConflict = set()  # list internalIdentifier, max 2 ICEsIMEsStructure In Conflict
        self.Length_of_blast_most_similar_ref_SP = -1 # -1 if NA
        self.Blast_ali_length = -1 # -1 if NA
        self.Blast_ali_start_CDS = -1 # -1 if NA
        self.Blast_ali_end_CDS = -1 # -1 if NA
        self.Blast_ali_start_Query_blast = -1 # -1 if NA
        self.Blast_ali_end_Query_blast = -1 # -1 if NA
        self.Blast_ali_identity_perc = -1 # -1 if NA
        self.E_value_blast = -1 # -1 if NA
        self.Blast_ali_bitscore = -1 # -1 if NA
        self.CDS_coverage_blast = -1 # -1 if NA
        self.Blast_ali_coverage_most_similar_ref_SP = -1 # -1 if NA


    def __hash__(self):
        return super(SP, self).__hash__()

    def __eq__(self, other):
        if not isinstance(other, SP):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return super(SP, self).__eq__(other)

    def __ne__(self, other):
        if not isinstance(other, SP):
            # don't attempt to compare against unrelated types
            return NotImplemented
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return super(SP, self).__ne__(other)


    def GetObjectAsJson(self):
        stToReturn = "{\n"
        stToReturn += "\t\"locusTag\": \"" + self.locusTag + "\"\n"
        stToReturn += "\t\"proteinId\": \"" + self.proteinId + "\"\n"
        stToReturn += "\t\"genomeAccession\": \"" + self.genomeAccession + "\"\n"
        stToReturn += "\t\"genomeAccessionRank\": \"" + self.genomeAccessionRank + "\"\n"
        stToReturn += "\t\"start\": \"" + str(self.start) + "\"\n"
        stToReturn += "\t\"stop\": \"" + str(self.stop) + "\"\n"
        stToReturn += "\t\"strand\": \"" + self.strand + "\"\n"
        stToReturn += "\t\"CDSPositionInGenome\": \"" + str(self.CDSPositionInGenome) + "\"\n"
        stToReturn += "\t\"pseudo\": \"" + str(self.pseudo) + "\"\n"
        stToReturn += "\t\"SPType\": \"" + self.SPType + "\"\n"
        stToReturn += "\t\"SPDetectedByBlast\": \"" + str(self.SPDetectedByBlast) + "\"\n"
        stToReturn += "\t\"SPDetectedByHMM\": \"" + str(self.SPDetectedByHMM) + "\"\n"
        # stToReturn += "\t\"setSPFamilyFromBlast\": \"" + ", ".join(str(i) for i in sorted(self.setSPFamilyFromBlast)) + "\"\n"
        stToReturn += "\t\"setSPICESuperFamilyFromBlast\": \"" + ", ".join(str(i) for i in sorted(self.setSPICESuperFamilyFromBlast)) + "\"\n"
        stToReturn += "\t\"setSPICEFamilyFromBlast\": \"" + ", ".join(str(i) for i in sorted(self.setSPICEFamilyFromBlast)) + "\"\n"
        stToReturn += "\t\"setSPIMESuperFamilyFromBlast\": \"" + ", ".join(str(i) for i in sorted(self.setSPIMESuperFamilyFromBlast)) + "\"\n"
        stToReturn += "\t\"Relaxase_family_domain_of_most_similar_ref_SPFromBlast\": \"" + str(self.Relaxase_family_domain_of_most_similar_ref_SPFromBlast) + "\"\n"
        stToReturn += "\t\"Relaxase_family_MOB_of_most_similar_ref_SPFromBlast\": \"" + str(self.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast) + "\"\n"
        stToReturn += "\t\"Coupling_type_of_most_similar_ref_SPFromBlast\": \"" + str(self.Coupling_type_of_most_similar_ref_SPFromBlast) + "\"\n"
        stToReturn += "\t\"False_positivesFromBlast\": \"" + str(self.False_positivesFromBlast) + "\"\n"
        stToReturn += "\t\"blast_validation\": \"" + str(self.blast_validation) + "\"\n"
        stToReturn += "\t\"Use_annotationFromBlast\": \"" + str(self.Use_annotationFromBlast) + "\"\n"
        stToReturn += "\t\"setSPFamilyFromHMM\": \"" + ", ".join(str(i) for i in sorted(self.setSPFamilyFromHMM)) + "\"\n"
        stToReturn += "\t\"setICEsIMEsStructureInConflict\": \"" + EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(self.setICEsIMEsStructureInConflict) + "\"\n"  # ", ".join(str(i) for i in sorted(self.setICEsIMEsStructureInConflict))
        stToReturn += "\t\"Length_of_blast_most_similar_ref_SP\": \"" + str(self.Length_of_blast_most_similar_ref_SP) + "\"\n"
        stToReturn += "\t\"Blast_ali_length\": \"" + str(self.Blast_ali_length) + "\"\n"
        stToReturn += "\t\"Blast_ali_start_CDS\": \"" + str(self.Blast_ali_start_CDS) + "\"\n"
        stToReturn += "\t\"Blast_ali_end_CDS\": \"" + str(self.Blast_ali_end_CDS) + "\"\n"
        stToReturn += "\t\"Blast_ali_start_Query_blast\": \"" + str(self.Blast_ali_start_Query_blast) + "\"\n"
        stToReturn += "\t\"Blast_ali_end_Query_blast\": \"" + str(self.Blast_ali_end_Query_blast) + "\"\n"
        stToReturn += "\t\"Blast_ali_identity_perc\": \"" + str(self.Blast_ali_identity_perc) + "\"\n"
        stToReturn += "\t\"E_value_blast\": \"" + str(self.E_value_blast) + "\"\n"
        stToReturn += "\t\"Blast_ali_bitscore\": \"" + str(self.Blast_ali_bitscore) + "\"\n"
        stToReturn += "\t\"CDS_coverage_blast\": \"" + str(self.CDS_coverage_blast) + "\"\n"
        stToReturn += "\t\"Blast_ali_coverage_most_similar_ref_SP\": \"" + str(self.Blast_ali_coverage_most_similar_ref_SP) + "\"\n"
        stToReturn += "}\n"
        return stToReturn


# The ListSPs object class is convenient to gather methods dedicated on manipulating list of SP object
class ListSPs():
    def __init__(self):
        self.list = []  # [SP]

    def sortListSPsByProximityToMEStructure(self, MEStructureSent):
        MostUpstreamSpFromMEStructureSentThatIsNotAnIntegrase = MEStructureSent.findMostUpstrSpNotIntegr()
        MostDownstreamSpFromMEStructureSentThatIsNotAnIntegrase = MEStructureSent.findMostDownstrSpNotIntegr()
        listOfTuplesIT = []
        for currIndex, currSP in enumerate(self.list, start=0):
            distanceCDSToMostUpstreamSpFromMEStructureSent = abs(MostUpstreamSpFromMEStructureSentThatIsNotAnIntegrase.CDSPositionInGenome - currSP.CDSPositionInGenome)
            distanceCDSToMostDownstreamSpFromMEStructureSent = abs(MostDownstreamSpFromMEStructureSentThatIsNotAnIntegrase.CDSPositionInGenome - currSP.CDSPositionInGenome)
            distanceCDSToMEStructureIT = distanceCDSToMostUpstreamSpFromMEStructureSent + distanceCDSToMostDownstreamSpFromMEStructureSent
            listOfTuplesIT.append((distanceCDSToMEStructureIT, currSP))
        # not useful when distanceCDSToMEStructureIT is equal, like for test#4 CIRMBP-1274_01370 and CIRMBP-1274_01372 ; it would be better to sort by another criterion such as evalue of the alignment ? For now the second attribute sorted is the CDSPositionInGenome to account for consenstency, but this is not ideal
        listOfTuplesIT.sort(key=lambda a: (a[0], a[1].CDSPositionInGenome))
        newListSPsIT = []  # [SP]
        for currIndex, (distanceCDSToMEStructureIT, currSP) in enumerate(listOfTuplesIT, start=0):
            newListSPsIT.append(currSP)
        self.list = newListSPsIT

    def sortListSPsByStart(self):
        #self.list.sort(key=lambda x: x.start, reverse=False)
        self.list.sort(key=lambda x: (x.genomeAccessionRank, x.start), reverse=False)

    
    def splitListOrderedSPsByColocalizion(self, maxNumberCDSForSplitSPsByColocalizion):
        listOfListOfColocalizedSPsToReturn = []
        currListSPs = ListSPs()
        for currIndex, currSP in enumerate(self.list, start=0):
            currListSPs.list.append(currSP)
            if (currIndex != (len(self.list) - 1)):  # if not last object in list
                otherSP = self.list[currIndex + 1]
                if otherSP.genomeAccessionRank != currSP.genomeAccessionRank:
                    # start a new currListSPs
                    listOfListOfColocalizedSPsToReturn.append(currListSPs)
                    currListSPs = ListSPs()
                else :
                    # distanceWithNextSp = currSP.GetDistanceCDSsToOtherSP(otherSP)
                    distanceWithNextSp = otherSP.CDSPositionInGenome - currSP.CDSPositionInGenome
                    if distanceWithNextSp <= 0:
                        raise RuntimeError(
                                "Error in splitListOrderedSPsByColocalizion: distanceWithNextSp <= 0: {} for SP {} and {}".format(
                                        str(distanceWithNextSp),
                                        currSP.locusTag,
                                        otherSP.locusTag))
                    if (distanceWithNextSp > maxNumberCDSForSplitSPsByColocalizion):
                        # start a new currListSPs
                        listOfListOfColocalizedSPsToReturn.append(currListSPs)
                        currListSPs = ListSPs()
            else:  # last object in list
                listOfListOfColocalizedSPsToReturn.append(currListSPs)
        return listOfListOfColocalizedSPsToReturn

    def fillUpIdxOfSPsInListSP(self):
        for idxCurrentSP, currentSP in enumerate(self.list):
            currentSP.idxInListSP = idxCurrentSP

    # return: [ICEsIMEsStructure] listICEsIMEsStructures
    # the aim of this method is to try to build ICEs and IMEs not taking into account nested structures. Merge of nested structures will be assessed and carried out at a latter stage in another method. This method first look sequentially for any SP of the conjugaison module (Relaxase, VirB4, Couplage) as a seed to start a new ICEsIMEsStructure. Then it tries to extend this ICEsIMEsStructure by adding other SP of the conjugaison module if relevant. Once the seed can not be extended further, it look for integrase that are appropriate and that could upstream or downstream. If the integrase could be downstream, do a check for potential conflict regarding upstream SPs of the conjugaison module that were attributed to the ICEsIMEsStructure previously built.
    def seedICEsIMEsStructure(
            self,
            groupListSPintoICEsIMEsUsingFamilyInfo,
            useFamilyInfoToTryToResolveSPModuleConjConflict,
            useDistanceCDSInfoToTryToResolveSPModuleConjConflict,
            allowAdjacentIntegraseOnlyForSer,
            locusTagIntegrase2Comment,
            SPsInSameFamilyMergeStructures2SameFamilyMergeStructure,
            listSameFamilyMergeStructures,
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):

        def registerCurrentICEsIMEsStructureIfNeeded(
                currentICEsIMEsStructure,
                listICEsIMEsStructures,
                idxCurrentSP,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):
            if currentICEsIMEsStructure.listOrderedSPs:
                currentICEsIMEsStructure.idxInSeedList = len(listICEsIMEsStructures)
                currentICEsIMEsStructure.addPotentialUpstreamSPInConflict(
                        self.list,
                        idxCurrentSP,
                        listICEsIMEsStructures,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        useFamilyInfoToTryToResolveSPModuleConjConflict,
                        useDistanceCDSInfoToTryToResolveSPModuleConjConflict,
                        SPsInSameFamilyMergeStructures2SameFamilyMergeStructure,
                        listSameFamilyMergeStructures)

                listICEsIMEsStructures.append(currentICEsIMEsStructure)
                # EMStructure.BasicEMStructure.countInternalIdentifier += 1
                currentICEsIMEsStructure = EMStructure.ICEsIMEsStructure(True)
                return currentICEsIMEsStructure
            else:
                return currentICEsIMEsStructure

        listICEsIMEsStructures = []
        currentICEsIMEsStructure = EMStructure.ICEsIMEsStructure(False)
        sameFamilyMergeStructureToCheck = None
        for idxCurrentSP, currentSP in enumerate(self.list):
            # print("dealing with currentSP {}".format(currentSP.locusTag))

            greenLightAddSPConjugaisonModule = False
            transfertCommentFromSameFamilyMergeStructure = False
            if currentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # print("\ncurrentSP: {}".format(currentSP.locusTag))
                sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP]
                greenLightAddSPConjugaisonModule = currentICEsIMEsStructure.listSPsIsContainedWithinOtherStructure(
                    currentSP,
                    sameFamilyMergeStructureToCheck)
                if greenLightAddSPConjugaisonModule:
                    greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currentICEsIMEsStructure,
                        currentSP,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        False,
                        False,
                        False)
                transfertCommentFromSameFamilyMergeStructure = True
            else:
                if sameFamilyMergeStructureToCheck is None:
                    greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currentICEsIMEsStructure,
                        currentSP,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        False,
                        False,
                        False)
                else:
                    # print("\ncurrentSP: {}".format(currentSP.locusTag))
                    greenLightAddSPConjugaisonModule = currentICEsIMEsStructure.listSPsIsContainedWithinOtherStructure(
                            currentSP,
                            sameFamilyMergeStructureToCheck)
                    if greenLightAddSPConjugaisonModule:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                            currentICEsIMEsStructure,
                            currentSP,
                            groupListSPintoICEsIMEsUsingFamilyInfo,
                            False,
                            False,
                            False)

            if greenLightAddSPConjugaisonModule:
                currentICEsIMEsStructure.addSPToConjugaisonModule(currentSP)
                if transfertCommentFromSameFamilyMergeStructure:
                    # transfer comment from merge structure to structure being built if yes
                    if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment not in currentICEsIMEsStructure.comment:
                        currentICEsIMEsStructure.comment += SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment
            else:
                # start new ICEsIMEsStructure and continue main loop
                # print("start new ICEsIMEsStructure and continue main loop")
                currentICEsIMEsStructure = registerCurrentICEsIMEsStructureIfNeeded(
                        currentICEsIMEsStructure,
                        listICEsIMEsStructures,
                        idxCurrentSP - 1,
                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
                # try adding again with new EM struct
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                    currentICEsIMEsStructure,
                    currentSP,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    False,
                    False,
                    False)
                if greenLightAddSPConjugaisonModule:
                    currentICEsIMEsStructure.addSPToConjugaisonModule(currentSP)
                    if currentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                        sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP]
                    else:
                        sameFamilyMergeStructureToCheck = None
                    if transfertCommentFromSameFamilyMergeStructure:
                        # transfer comment from merge structure to structure being built if yes
                        if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment not in currentICEsIMEsStructure.comment:
                            currentICEsIMEsStructure.comment += SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment

        # last iteration
        registerCurrentICEsIMEsStructureIfNeeded(
                currentICEsIMEsStructure,
                listICEsIMEsStructures,
                len(self.list) - 1,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
        return listICEsIMEsStructures

    @staticmethod
    def GetIdxSPInList(SPToFind, listSPs):
        for idx, SPInList in enumerate(listSPs):
            if SPToFind == SPInList: # was if SPToFind.locusTag == SPInList.locusTag:
                return idx
        return -1

    @staticmethod
    def GetListCDSNumberFromListSP(listSPs):
        strToReturn = ""
        for currentSP in listSPs:
            if strToReturn:  # not tempty
                strToReturn += ", "
            strToReturn += str(currentSP.CDSPositionInGenome)
        return strToReturn

    # @staticmethod
    def getListIntegraseGroupJustUpstreamOfThisSP(SPSent, listSPsSent, dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd, setIntegraseTypeToCheck, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment):
        listUpstreamIntegraseToReturn = []
        
        if SPSent not in listSPsSent :
            raise RuntimeError("Error in getListIntegraseGroupJustUpstreamOfThisSP: SPSent {} not in listSPsSent {}".format(SPSent.locusTag, ListSPs.GetListCDSNumberFromListSP(listSPsSent)))

        idxSPsUpstreamToCheck = SPSent.idxInListSP - 1
        if idxSPsUpstreamToCheck >= 0:
            listSPsUpstream = listSPsSent[:(idxSPsUpstreamToCheck + 1)]
            for currSp in reversed(listSPsUpstream):
                if dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd is not None :
                    if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                        continue
                        
                # print("HERE up {}".format(currSp.locusTag))
                if listUpstreamIntegraseToReturn:
                    if (abs(listUpstreamIntegraseToReturn[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                        if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, True)
                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                #setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "Integrase {} can not be associated with the downstream structure {}, it needs to be on the - strand. ".format(currSp.locusTag, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier)
                                if strCommentIT not in downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                    downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                strCommentIT = "If {} is part of the structure {}, the upstream integrase {} would not be associated with this structure as it needs to be on the - strand. ".format(
                                        ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier, currSp.locusTag)
                                if strCommentIT not in downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                    downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                            if currSp.SPType == listUpstreamIntegraseToReturn[-1].SPType and currSp.strand == listUpstreamIntegraseToReturn[-1].strand:

                                listUpstreamIntegraseToReturn.append(currSp)

                                continue
                            else:
                                #setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "The upstream integrases {} and {} are adjacent but not of the same type or strand, {} is not considered as part of a tandem.".format(listUpstreamIntegraseToReturn[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                if strCommentIT not in downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                    downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                        else:
                            break
                    else:  # SP are not adjacent in genome
                        break
                else:  # listUpstreamIntegraseToReturn is empty
                    if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                        valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, True)
                        if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                            #setIntegraseDiscarded.add(currSp.locusTag)
                            strCommentIT = "Integrase {} can not be associated with the downstream structure {}, it needs to be on the - strand. ".format(currSp.locusTag, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier)
                            if strCommentIT not in downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                            icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                            break
                        elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                            strCommentIT = "If {} is part of the structure {}, the upstream integrase {} would not be associated with this structure as it needs to be on the - strand. ".format(
                                    ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier, currSp.locusTag)
                            if strCommentIT not in downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                            for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                        listUpstreamIntegraseToReturn.append(currSp)
                        continue
                    else:
                        break

        return listUpstreamIntegraseToReturn


    @staticmethod
    def getListIntegraseGroupJustDownstreamOfThisSP(SPSent, listSPsSent, dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd, setIntegraseTypeToCheck, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment):
        listDownstreamIntegraseToReturn = []

        if SPSent not in listSPsSent :
            raise RuntimeError("Error in getListIntegraseGroupJustDownstreamOfThisSP: SPSent {} not in listSPsSent {}".format(SPSent.locusTag, ListSPs.GetListCDSNumberFromListSP(listSPsSent)))

        idxSPsDownstreamToCheck = SPSent.idxInListSP + 1
        if idxSPsDownstreamToCheck < len(listSPsSent):
            listSPsDownstream = listSPsSent[idxSPsDownstreamToCheck:]
            for currSp in listSPsDownstream:
                if dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd is not None :
                    if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                        continue

                #  print("HERE down {}".format(currSp.locusTag))
                if listDownstreamIntegraseToReturn:  # listIntegraseDownstream is not empty
                    if (abs(listDownstreamIntegraseToReturn[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                        if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                            if upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented is not None :
                                valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, False)
                                if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                    #setIntegraseDiscarded.add(currSp.locusTag)
                                    strCommentIT = "Integrase {} can not be associated with the upstream structure {}, it needs to be on the + strand. ".format(currSp.locusTag, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier)
                                    if strCommentIT not in upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                        upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    break
                                elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                    strCommentIT = "If {} is part of the structure {}, the downstream integrase {} would not be associated with this structure as it needs to be on the + strand. ".format(
                                            ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier, currSp.locusTag)
                                    if strCommentIT not in upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                        upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                    for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                        icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                            if currSp.SPType == listDownstreamIntegraseToReturn[-1].SPType and currSp.strand == listDownstreamIntegraseToReturn[-1].strand:

                                listDownstreamIntegraseToReturn.append(currSp)

                                continue
                            else:
                                
                                #setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "The downstream integrases {} and {} are adjacent but not of the same type or strand, {} is not considered as part of a tandem. ".format(listDownstreamIntegraseToReturn[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                if upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented is not None :
                                    if strCommentIT not in upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                        upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                        else:
                            break
                    else:  # SP are not adjacent in genome
                        break
                else:  # listDownstreamIntegraseToReturn is empty

                    if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                        if upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented is not None :
                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, False)
                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                #setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "Integrase {} can not be associated with the upstream structure {}, it needs to be on the + strand. ".format(currSp.locusTag, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier)
                                if strCommentIT not in upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                    upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                strCommentIT = "If {} is part of the structure {}, the downstream integrase {} would not be associated with this structure as it needs to be on the + strand. ".format(
                                        ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.internalIdentifier, currSp.locusTag)
                                if strCommentIT not in upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment:
                                    upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented.comment += strCommentIT
                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                        listDownstreamIntegraseToReturn.append(currSp)
                        continue
                    else:
                        break

        return listDownstreamIntegraseToReturn


    @staticmethod
    def GetListProtIdsFromSetSP(setSPs):
        strToReturn = ""
        for currentSP in sorted(setSPs, key=lambda x: x.locusTag, reverse=False):
            if strToReturn:  # not tempty
                strToReturn += ", "
            strToReturn += currentSP.locusTag
            if currentSP.pseudo :
                strToReturn += " [Pseudo]"
        return strToReturn


    @staticmethod
    def GetListProtIdsFromListSP(listSPs):
        strToReturn = ""
        for currentSP in listSPs:
            if strToReturn:  # not tempty
                strToReturn += ", "
            strToReturn += currentSP.locusTag
            if currentSP.pseudo :
                strToReturn += " [Pseudo]"
        return strToReturn


    @staticmethod
    def stringListSPHasPseudo(listSPs):
        stringListSPHasPseudoToReturn = ""
        for currIndex, currSP in enumerate(listSPs, start=0):
            if currSP.pseudo:
                if len(listSPs) > 1 :
                    stringListSPHasPseudoToReturn = " with Pseudo"
                else :
                    stringListSPHasPseudoToReturn = " Pseudo"
        return stringListSPHasPseudoToReturn




    @staticmethod
    def PrintICElineFormat(listSPs, onlyReturnListFirstLetterSPWithoutDistanceBetween, strToAppendIfMoreThanOneSP):
        # example ICElineFormat: C4.V5.C2.R7.V1:I1:R4.I1!i2!i7.R1:v6.c5!R2.I2!R2.I
        # SP type: R = Relaxase, C = Coupling Protein, V = VirB4, I = Integrase, D = DDE
        # SP strand: R = Relaxase on forward strand, r = relaxase on reverse strand
        # CDS distance between each SP:        .=units,:=dozens, !=hundreds
        strToReturn = ""
        for currIndex, currSP in enumerate(listSPs, start=0):
            if (currIndex != (len(listSPs) - 1)):  # if not last object in list
                nextSP = listSPs[currIndex + 1]
                distance = abs(currSP.CDSPositionInGenome - nextSP.CDSPositionInGenome) - 1

                # strCurrSP = currSP.SPType[0]
                if (currSP.SPType == "DDE transposase"):
                    strCurrSP = "D"
                elif (currSP.SPType in icescreen_OO.setIntegraseNames):
                    strCurrSP = "I"
                elif (currSP.SPType == "Relaxase"):
                    strCurrSP = "R"
                elif (currSP.SPType == "Coupling protein"):
                    strCurrSP = "C"
                elif (currSP.SPType == "VirB4"):
                    strCurrSP = "V"
                else:
                    raise RuntimeError("Error in PrintICElineFormat: unrecognized currSP.SPType = {}".format(currSP.SPType))

                strDistanceNumber = ""
                strDistanceLegend = ""
                separatorBetweenSP = ""
                if not onlyReturnListFirstLetterSPWithoutDistanceBetween :

                    if currSP.strand == "+":
                        strCurrSP = strCurrSP.upper()
                    else:
                        strCurrSP = strCurrSP.lower()

                    strDistanceNumber = str(distance)[0]  # first char

                    if distance < 10:
                        strDistanceLegend = "."
                    elif distance < 100:
                        strDistanceLegend = ":"
                    else:
                        strDistanceLegend = "!"
                else :
                    separatorBetweenSP = "+"

                strPseudoOrNot = ""
                if currSP.pseudo:
                    strPseudoOrNot = "[Pseudo]"

                strToReturn += strCurrSP + strPseudoOrNot + strDistanceNumber + strDistanceLegend + separatorBetweenSP

            else:  # if last object in list
                if (currSP.SPType == "DDE transposase"):
                    strCurrSP = "D"
                elif (currSP.SPType in icescreen_OO.setIntegraseNames):
                    strCurrSP = "I"
                elif (currSP.SPType == "Relaxase"):
                    strCurrSP = "R"
                elif (currSP.SPType == "Coupling protein"):
                    strCurrSP = "C"
                elif (currSP.SPType == "VirB4"):
                    strCurrSP = "V"
                else:
                    raise RuntimeError("Error in PrintICElineFormat: unrecognized currSP.SPType = {}".format(currSP.SPType))

                if not onlyReturnListFirstLetterSPWithoutDistanceBetween :
                    if currSP.strand == "+":
                        strCurrSP = strCurrSP.upper()
                    else:
                        strCurrSP = strCurrSP.lower()

                strPseudoOrNot = ""
                if currSP.pseudo:
                    strPseudoOrNot = "[Pseudo]"

                strToReturn += strCurrSP + strPseudoOrNot


        if len(listSPs) > 1 :
            strToReturn += strToAppendIfMoreThanOneSP

        return strToReturn
