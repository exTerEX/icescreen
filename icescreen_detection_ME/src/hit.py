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
import commonMethods
import itertools
import traceback

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
    
    def __lt__(self, other):
        if not isinstance(other, CDS):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.start < other.start
    

# The SP object class extends CDS and add extra attributes to mimic the biological object signature protein
class SP(CDS):
    def __init__(self):
        super(SP, self).__init__()
        self.SPType = ""  # col CDS_Protein_Type ; Values: icescreen_OO.listTypeSPConjModule, icescreen_OO.setIntegrase # WAS Values: Coupling, Relaxase, Virb4
        self.SPDetectedByBlast = ""  # col Is_hit_blast ; Values: 1, 0
        self.SPDetectedByHMM = ""  # col Is_hit_HMM_JL or Is_hit_HMM_CC depending ; Values: 1, 0
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

        self.setICEsIMEsStructureInConflict = set() # list internalIdentifier, max 2 ICEsIMEsStructure In Conflict
        self.listSiblingFragmentedSP = [] # contains all the sibling fragmented SP includng self, sorted by CDSPositionInGenome

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

    def __lt__(self, other):
        if not isinstance(other, SP):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return super(SP, self).__lt__(other)


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

        stToReturn += "\t\"setICEsIMEsStructureInConflict\": \"" + EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(self.setICEsIMEsStructureInConflict) + "\"\n"
        stToReturn += "\t\"listSiblingFragmentedSP\": \"" + ListSPs.GetListProtIdsFromListSP(self.listSiblingFragmentedSP) + "\"\n"

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
        self.list.sort(key=lambda x: (x.genomeAccessionRank, x.start), reverse=False)


    def splitListOrderedSPs_colocalizeByFragmentedSPs(self, locusTagMerge2Comment):
        listOfListOfColocalizedSPsToReturn = []
        tmpListOfListOFSPToReassambleUsingFragmentedSPs = []
        tmpListOfListOFSPToReassambleUsingFragmentedSPs.append(self.list)

        index_tmpListOFSPToReassambleUsingFragmentedSPsIT = 0
        while index_tmpListOFSPToReassambleUsingFragmentedSPsIT < len(tmpListOfListOFSPToReassambleUsingFragmentedSPs):
            tmpListOFSPToReassambleUsingFragmentedSPsIT = tmpListOfListOFSPToReassambleUsingFragmentedSPs[index_tmpListOFSPToReassambleUsingFragmentedSPsIT]
            indexInMainSPList = 0
            builtUpListOfColocalizedSPsToReturn = ListSPs()
            while indexInMainSPList < len(tmpListOFSPToReassambleUsingFragmentedSPsIT):
                SPToReassambleUsingFragmentedSPsIT = tmpListOFSPToReassambleUsingFragmentedSPsIT[indexInMainSPList]
                if len(SPToReassambleUsingFragmentedSPsIT.listSiblingFragmentedSP) > 0:
                    #found a fragmented SP to reassamble
                    for indexFragmentedSPToReassambleIT, fragmentedSPToReassambleIT in enumerate(SPToReassambleUsingFragmentedSPsIT.listSiblingFragmentedSP):
                        builtUpListOfColocalizedSPsToReturn.list.append(fragmentedSPToReassambleIT)
                        if indexFragmentedSPToReassambleIT == 0:
                            if fragmentedSPToReassambleIT != SPToReassambleUsingFragmentedSPsIT:
                                raise RuntimeError("Error in splitListOrderedSPs_colocalizeByFragmentedSPs: fragmentedSPToReassambleIT {} != SPToReassambleUsingFragmentedSPsIT {} ; list SPs = {}; initial self.list = {}".format(
                                    fragmentedSPToReassambleIT.locusTag,
                                    SPToReassambleUsingFragmentedSPsIT.locusTag,
                                    ListSPs.GetListProtIdsFromListSP(tmpListOFSPToReassambleUsingFragmentedSPsIT),
                                    ListSPs.GetListProtIdsFromListSP(self.list)
                                    ))
                        else:
                            try:
                                indexNextFragmentToReassembleInMainSPList = tmpListOFSPToReassambleUsingFragmentedSPsIT.index(fragmentedSPToReassambleIT, indexInMainSPList)
                                sliceOfInnerSP = tmpListOFSPToReassambleUsingFragmentedSPsIT[(indexInMainSPList+1):indexNextFragmentToReassembleInMainSPList]
                                tmpListOfListOFSPToReassambleUsingFragmentedSPs.append(sliceOfInnerSP)
                                indexInMainSPList = indexNextFragmentToReassembleInMainSPList
                            except ValueError:
                                #fragmentedSPToReassambleIT was not found
                                commentITToAdd = "The fragment {} was not found when looking to it assemble it with its sibling fragments ({}) in the subsegment {}. This may lead to inconsistency, merging of fragmented SPs is therefore not considered in this segment. ".format(
                                    fragmentedSPToReassambleIT.locusTag,
                                    ListSPs.GetListProtIdsFromListSP(SPToReassambleUsingFragmentedSPsIT.listSiblingFragmentedSP),
                                    ListSPs.GetListProtIdsFromListSP(tmpListOFSPToReassambleUsingFragmentedSPsIT)
                                    )
                                # tmpListOFSPToReassambleUsingFragmentedSPsIT.index(fragmentedSPToReassambleIT, indexInMainSPList), not considering merging fragmented SP for this segment
                                icescreen_OO.addCommentToLocusTag2Comment(fragmentedSPToReassambleIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                                # return list with just self.list
                                listOfListOfColocalizedSPsToReturn.clear()
                                listOfListOfColocalizedSPsToReturn.append(self)
                                return listOfListOfColocalizedSPsToReturn
                            except Exception as caughtExceptionIT:
                                raise RuntimeError("Error in splitListOrderedSPs_colocalizeByFragmentedSPs: list SPs = {} ; initial self.list = {} ; try catch return an error : {}".format(
                                    ListSPs.GetListProtIdsFromListSP(tmpListOFSPToReassambleUsingFragmentedSPsIT),
                                    ListSPs.GetListProtIdsFromListSP(self.list),
                                    traceback.format_exc()
                                    ))
                else :
                    builtUpListOfColocalizedSPsToReturn.list.append(SPToReassambleUsingFragmentedSPsIT)
                #increment the loop
                indexInMainSPList += 1
            if len(builtUpListOfColocalizedSPsToReturn.list) > 0:
                listOfListOfColocalizedSPsToReturn.append(builtUpListOfColocalizedSPsToReturn)
            index_tmpListOFSPToReassambleUsingFragmentedSPsIT += 1
            if index_tmpListOFSPToReassambleUsingFragmentedSPsIT > 10000:
                raise RuntimeError("Error in splitListOrderedSPs_colocalizeByFragmentedSPs: index_tmpListOFSPToReassambleUsingFragmentedSPsIT > 10000 ; initial self.list = {}".format(
                    ListSPs.GetListProtIdsFromListSP(self.list)
                    ))
        return listOfListOfColocalizedSPsToReturn


    
    def splitListOrderedSPs_colocalizeByMaxNumberCDS(self):
        #commonMethods.ConfigParams.maxNumberCDSForSplitSPsByColocalizion
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
                                "Error in splitListOrderedSPs_colocalizeByMaxNumberCDS: distanceWithNextSp <= 0: {} for SP {} and {}".format(
                                        str(distanceWithNextSp),
                                        currSP.locusTag,
                                        otherSP.locusTag))
                    if (distanceWithNextSp > commonMethods.ConfigParams.maxNumberCDSForSplitSPsByColocalizion):
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


    def searchForFragmentedSPs(
        self,
        locusTagMerge2Comment,
        locusTagIntegrase2Comment,
        ):

        printInDebugIT = False

        SPTypeIT2setComplementarySPTypeIT = {}
        for (SPTypeIT1IT, SPTypeIT2IT) in list(itertools.combinations(self.list, 2)) :
            if printInDebugIT:
                print("searchForFragmentedSPs: {} and {}".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag))

            if SPTypeIT1IT.SPType != SPTypeIT2IT.SPType:
                if printInDebugIT:
                    print("SPTypeIT1IT.SPType != SPTypeIT2IT.SPType: {} and {}".format(SPTypeIT1IT.SPType, SPTypeIT2IT.SPType))
                continue

            if SPTypeIT1IT.strand != SPTypeIT2IT.strand:
                if printInDebugIT:
                    print("SPTypeIT1IT.strand != SPTypeIT2IT.strand: {} and {}".format(SPTypeIT1IT.strand, SPTypeIT2IT.strand))
                continue

            if SPTypeIT1IT.Blast_ali_start_Query_blast == -1 or SPTypeIT1IT.Blast_ali_end_Query_blast == -1 or SPTypeIT2IT.Blast_ali_start_Query_blast == -1 or SPTypeIT2IT.Blast_ali_end_Query_blast == -1:
                if printInDebugIT:
                    print("SPTypeIT1IT.Blast_ali_start_Query_blast == -1 ({}) or SPTypeIT1IT.Blast_ali_end_Query_blast == -1 ({}) or SPTypeIT2IT.Blast_ali_start_Query_blast == -1 ({}) or SPTypeIT2IT.Blast_ali_end_Query_blast == -1 ({})".format(str(SPTypeIT1IT.Blast_ali_start_Query_blast), str(SPTypeIT1IT.Blast_ali_end_Query_blast), str(SPTypeIT2IT.Blast_ali_start_Query_blast), str(SPTypeIT2IT.Blast_ali_end_Query_blast)))
                continue


            overlapIT = commonMethods.overlap_bounded_intervals(
                SPTypeIT1IT.Blast_ali_start_Query_blast,
                SPTypeIT1IT.Blast_ali_end_Query_blast,
                SPTypeIT2IT.Blast_ali_start_Query_blast,
                SPTypeIT2IT.Blast_ali_end_Query_blast,
                SPTypeIT1IT.genomeAccessionRank,
                SPTypeIT2IT.genomeAccessionRank
                )
            if overlapIT >= commonMethods.ConfigParams.maxOverlappingAliLenghtFragmentedSPs :
                #Those 2 SPs do not overlap
                if printInDebugIT:
                    print("overlapIT ({}) >= maxOverlappingAliLenghtFragmentedSPs ({})".format(str(overlapIT), str(commonMethods.ConfigParams.maxOverlappingAliLenghtFragmentedSPs)))
                continue

            # do not merge integrase of different types : tyr, ser, DDE
            if SPTypeIT1IT.SPType in icescreen_OO.setIntegraseNames and SPTypeIT2IT.SPType in icescreen_OO.setIntegraseNames:
                if SPTypeIT1IT.SPType in icescreen_OO.setIntegraseTyrNames and SPTypeIT2IT.SPType not in icescreen_OO.setIntegraseTyrNames :
                    commentITToAdd = "{} could complement the fragment {} to form a single integrase but they are not of the same type ({} and {}). ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT2IT.SPType)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    if printInDebugIT:
                        print(commentITToAdd)
                    continue
                elif SPTypeIT1IT.SPType in icescreen_OO.setIntegraseSerNames and SPTypeIT2IT.SPType not in icescreen_OO.setIntegraseSerNames :
                    commentITToAdd = "{} could complement the fragment {} to form a single integrase but they are not of the same type ({} and {}). ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT2IT.SPType)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    if printInDebugIT:
                        print(commentITToAdd)
                    continue
                elif SPTypeIT1IT.SPType in icescreen_OO.setIntegraseDDENames and SPTypeIT2IT.SPType not in icescreen_OO.setIntegraseDDENames :
                    commentITToAdd = "{} could complement the fragment {} to form a single integrase but they are not of the same type ({} and {}). ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT2IT.SPType)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    if printInDebugIT:
                        print(commentITToAdd)
                    continue



            # do not merge fragmented intergrase if triple serine in between
            if SPTypeIT1IT.SPType in icescreen_OO.setIntegraseSerNames and SPTypeIT2IT.SPType in icescreen_OO.setIntegraseSerNames:
                if SPTypeIT1IT.CDSPositionInGenome < SPTypeIT2IT.CDSPositionInGenome:
                    # downstream SP of SPTypeIT1IT
                    if self.list[SPTypeIT1IT.idxInListSP+1].SPType in icescreen_OO.setIntegraseSerNames :
                        if (abs(self.list[SPTypeIT1IT.idxInListSP+1].CDSPositionInGenome - SPTypeIT1IT.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            commentITToAdd = "{} could complement the fragment {} to form a single {} but there is a downstream Integrases {} adjacent in the genome in-between, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, self.list[SPTypeIT1IT.idxInListSP+1].locusTag)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            if printInDebugIT:
                                print(commentITToAdd)
                            continue
                    # upstream SP of SPTypeIT2IT
                    if self.list[SPTypeIT2IT.idxInListSP-1].SPType in icescreen_OO.setIntegraseSerNames :
                        if (abs(self.list[SPTypeIT2IT.idxInListSP-1].CDSPositionInGenome - SPTypeIT2IT.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            commentITToAdd = "{} could complement the fragment {} to form a single {} but there is a upstream Integrases {} adjacent in the genome in-between, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, self.list[SPTypeIT2IT.idxInListSP-1].locusTag)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            if printInDebugIT:
                                print(commentITToAdd)
                            continue
                else:
                    # downstream SP of SPTypeIT2IT
                    if self.list[SPTypeIT2IT.idxInListSP+1].SPType in icescreen_OO.setIntegraseSerNames :
                        if (abs(self.list[SPTypeIT2IT.idxInListSP+1].CDSPositionInGenome - SPTypeIT2IT.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            commentITToAdd = "{} could complement the fragment {} to form a single {} but there is a downstream Integrases {} adjacent in the genome in-between, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, self.list[SPTypeIT2IT.idxInListSP+1].locusTag)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            if printInDebugIT:
                                print(commentITToAdd)
                            continue
                    # upstream SP of SPTypeIT1IT
                    if self.list[SPTypeIT1IT.idxInListSP-1].SPType in icescreen_OO.setIntegraseSerNames :
                        if (abs(self.list[SPTypeIT1IT.idxInListSP-1].CDSPositionInGenome - SPTypeIT1IT.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            commentITToAdd = "{} could complement the fragment {} to form a single {} but there is a upstream Integrases {} adjacent in the genome in-between, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, self.list[SPTypeIT1IT.idxInListSP-1].locusTag)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            if printInDebugIT:
                                print(commentITToAdd)
                            continue

            # check setSPICEFamilyFromBlast compatible
            if SPTypeIT1IT.SPType not in icescreen_OO.setIntegraseNames:
                if len(SPTypeIT1IT.setSPICESuperFamilyFromBlast) > 0 and len(SPTypeIT2IT.setSPICESuperFamilyFromBlast) > 0:
                    if len(SPTypeIT1IT.setSPICESuperFamilyFromBlast.intersection(SPTypeIT2IT.setSPICESuperFamilyFromBlast)) == 0:
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but they have different ICE superfamily, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("len(SPTypeIT1IT.setSPICESuperFamilyFromBlast.intersection(SPTypeIT2IT.setSPICESuperFamilyFromBlast)) == 0 ; len(SPTypeIT1IT.setSPICESuperFamilyFromBlast = {} ; SPTypeIT2IT.setSPICESuperFamilyFromBlast = {} ".format(str(SPTypeIT1IT.setSPICESuperFamilyFromBlast), str(SPTypeIT2IT.setSPICESuperFamilyFromBlast)))
                        continue
                if SPTypeIT1IT.SPType == "Relaxase" :
                    if SPTypeIT1IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast != SPTypeIT2IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast :
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but their most similar reference SP have different relaxase family domain, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("SPTypeIT1IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast ({}) != SPTypeIT2IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast ({})".format(str(SPTypeIT1IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast), str(SPTypeIT2IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast)))
                        continue
            
            #Those 2 SPs do not overlap much, they may be complementary
            #check that the merge orientation is correct : i.e. can not merge if frgagment are the wrong way


            if SPTypeIT1IT.strand == "+" and SPTypeIT2IT.strand == "+" :
                # SP1   4.0	324.0
                # SP2   328.0	579.0
                if SPTypeIT1IT.CDSPositionInGenome < SPTypeIT2IT.CDSPositionInGenome:
                    if SPTypeIT1IT.Blast_ali_start_Query_blast > SPTypeIT2IT.Blast_ali_start_Query_blast or SPTypeIT1IT.Blast_ali_start_Query_blast > SPTypeIT1IT.Blast_ali_end_Query_blast or SPTypeIT2IT.Blast_ali_start_Query_blast > SPTypeIT2IT.Blast_ali_end_Query_blast :
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but the orientation of {} and {} is not coherent, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("SPTypeIT1IT.Blast_ali_start_Query_blast ({}) > SPTypeIT2IT.Blast_ali_start_Query_blast ({})".format(str(SPTypeIT1IT.Blast_ali_start_Query_blast), str(SPTypeIT2IT.Blast_ali_start_Query_blast)))
                        continue
                else :
                    if SPTypeIT2IT.Blast_ali_start_Query_blast > SPTypeIT1IT.Blast_ali_start_Query_blast or SPTypeIT1IT.Blast_ali_start_Query_blast > SPTypeIT1IT.Blast_ali_end_Query_blast or SPTypeIT2IT.Blast_ali_start_Query_blast > SPTypeIT2IT.Blast_ali_end_Query_blast :
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but the orientation of {} and {} is not coherent, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("SPTypeIT2IT.Blast_ali_start_Query_blast ({}) > SPTypeIT1IT.Blast_ali_start_Query_blast ({})".format(str(SPTypeIT2IT.Blast_ali_start_Query_blast), str(SPTypeIT1IT.Blast_ali_start_Query_blast)))
                        continue
            elif SPTypeIT1IT.strand == "-" and SPTypeIT2IT.strand == "-" :
                # SP1   328.0	579.0
                # SP2   4.0	324.0
                if SPTypeIT1IT.CDSPositionInGenome < SPTypeIT2IT.CDSPositionInGenome:
                    if SPTypeIT1IT.Blast_ali_start_Query_blast < SPTypeIT2IT.Blast_ali_start_Query_blast or SPTypeIT1IT.Blast_ali_start_Query_blast > SPTypeIT1IT.Blast_ali_end_Query_blast or SPTypeIT2IT.Blast_ali_start_Query_blast > SPTypeIT2IT.Blast_ali_end_Query_blast :
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but the orientation of {} and {} is not coherent, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("SPTypeIT1IT.Blast_ali_start_Query_blast ({}) > SPTypeIT2IT.Blast_ali_start_Query_blast ({})".format(str(SPTypeIT1IT.Blast_ali_start_Query_blast), str(SPTypeIT2IT.Blast_ali_start_Query_blast)))
                        continue
                else :
                    if SPTypeIT2IT.Blast_ali_start_Query_blast < SPTypeIT1IT.Blast_ali_start_Query_blast or SPTypeIT1IT.Blast_ali_start_Query_blast > SPTypeIT1IT.Blast_ali_end_Query_blast or SPTypeIT2IT.Blast_ali_start_Query_blast > SPTypeIT2IT.Blast_ali_end_Query_blast :
                        commentITToAdd = "{} could complement the fragment {} to form a single {} but the orientation of {} and {} is not coherent, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, SPTypeIT1IT.SPType, SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if printInDebugIT:
                            print("SPTypeIT2IT.Blast_ali_start_Query_blast ({}) > SPTypeIT1IT.Blast_ali_start_Query_blast ({})".format(str(SPTypeIT2IT.Blast_ali_start_Query_blast), str(SPTypeIT1IT.Blast_ali_start_Query_blast)))
                        continue


            else :
                raise RuntimeError(
                    "Error in searchForFragmentedSPs: SPTypeIT1IT.strand == {} and SPTypeIT2IT.strand == {}".format(
                            str(SPTypeIT1IT.strand),
                            str(SPTypeIT2IT.strand)))


            if SPTypeIT1IT in SPTypeIT2setComplementarySPTypeIT :
                setComplementarySPTypeIT = SPTypeIT2setComplementarySPTypeIT[SPTypeIT1IT]
                setComplementarySPTypeIT.add(SPTypeIT2IT)
            else :
                setComplementarySPTypeIT = set()
                setComplementarySPTypeIT.add(SPTypeIT2IT)
                SPTypeIT2setComplementarySPTypeIT[SPTypeIT1IT] = setComplementarySPTypeIT


        for SPTypeIT1IT, setComplementarySPTypeIT1IT in SPTypeIT2setComplementarySPTypeIT.items():
            # check that the other setComplementary are coherent
            booleanAllComplementarySPTypeITAreCoherent = True
            for SPTypeIT2IT in setComplementarySPTypeIT1IT :
                if SPTypeIT2IT in SPTypeIT2setComplementarySPTypeIT :
                    setComplementarySPTypeIT2IT = SPTypeIT2setComplementarySPTypeIT[SPTypeIT2IT]
                    for SPTypeITToAssertForCoherence in setComplementarySPTypeIT2IT :
                        if SPTypeITToAssertForCoherence not in setComplementarySPTypeIT1IT :
                            booleanAllComplementarySPTypeITAreCoherent = False
                            # SPTypeIT1IT could complement the fragments setComplementarySPTypeIT1IT to form a single SPTypeIT SP but SPTypeIT2IT was found to complement the fragment SPTypeITToAssertForCoherence as well, please manually review.
                            commentITToAdd = "{} could complement the fragments {} to form a single {} but {} was found to complement the fragment {} as well, please manually review. ".format(SPTypeIT1IT.locusTag, ListSPs.GetListProtIdsFromSetSP(setComplementarySPTypeIT1IT), SPTypeIT1IT.SPType, SPTypeIT2IT.locusTag, SPTypeITToAssertForCoherence.locusTag)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)


            if booleanAllComplementarySPTypeITAreCoherent :
                # if more than SP complementary, check that the combination is complementary
                listComplementarySPsToCheckFurther = []
                listComplementarySPsToCheckFurther.append(SPTypeIT1IT)
                listComplementarySPsToCheckFurther.extend(setComplementarySPTypeIT1IT)
                atLeastOneOverlapITGreaterEqualThanMaxAllowed = False
                for (SPTypeIT1ToCheckFurtherIT, SPTypeIT2ToCheckFurtherIT) in list(itertools.combinations(listComplementarySPsToCheckFurther, 2)) :
                    if SPTypeIT1ToCheckFurtherIT.Blast_ali_start_Query_blast == -1 or SPTypeIT1ToCheckFurtherIT.Blast_ali_end_Query_blast == -1 or SPTypeIT2ToCheckFurtherIT.Blast_ali_start_Query_blast == -1 or SPTypeIT2ToCheckFurtherIT.Blast_ali_end_Query_blast == -1:
                        raise RuntimeError(
                            "Error in searchForFragmentedSPs: SPTypeIT1ToCheckFurtherIT {} or SPTypeIT2ToCheckFurtherIT have Blast_ali__Query_blast  == -1".format(SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag))
                    overlapIT = commonMethods.overlap_bounded_intervals(
                        SPTypeIT1ToCheckFurtherIT.Blast_ali_start_Query_blast,
                        SPTypeIT1ToCheckFurtherIT.Blast_ali_end_Query_blast,
                        SPTypeIT2ToCheckFurtherIT.Blast_ali_start_Query_blast,
                        SPTypeIT2ToCheckFurtherIT.Blast_ali_end_Query_blast,
                        SPTypeIT1ToCheckFurtherIT.genomeAccessionRank,
                        SPTypeIT2ToCheckFurtherIT.genomeAccessionRank
                        )
                    if overlapIT >= commonMethods.ConfigParams.maxOverlappingAliLenghtFragmentedSPs :
                        atLeastOneOverlapITGreaterEqualThanMaxAllowed = True
                        # listComplementarySPsToCheckFurther were checked to complement each other as fragments to form a single SPTypeIT SP, but SPTypeIT1ToCheckFurtherIT was found to overlap with SPTypeIT2ToCheckFurtherIT, please manually review.
                        commentITToAdd = "{} were checked to complement each other as fragments to form a single {}, but {} was found to overlap with {}, please manually review. ".format(ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), SPTypeIT1IT.SPType, SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)

                if not atLeastOneOverlapITGreaterEqualThanMaxAllowed :
                    #OK those listComplementarySPsToCheckFurther can be checked to be merged together
                    # checking for maxCDSDistanceToMergeFragmentedSPs = 150
                    atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs = False
                    for (SPTypeIT1ToCheckFurtherIT, SPTypeIT2ToCheckFurtherIT) in list(itertools.combinations(listComplementarySPsToCheckFurther, 2)) :
                        if abs(SPTypeIT2ToCheckFurtherIT.CDSPositionInGenome - SPTypeIT1ToCheckFurtherIT.CDSPositionInGenome) > commonMethods.ConfigParams.maxCDSDistanceToMergeFragmentedSPs :
                            atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs = True
                            # listComplementarySPsToCheckFurther were checked to complement each other as fragments to form a single SPTypeIT SP, but SPTypeIT1ToCheckFurtherIT was found to be further away to SPTypeIT2ToCheckFurtherIT than cutoff, please manually review.
                            commentITToAdd = "{} were checked to complement each other as fragments to form a single {}, but {} was found to be further away to {} than cutoff {}, please manually review. ".format(ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), SPTypeIT1IT.SPType, SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag, str(commonMethods.ConfigParams.maxCDSDistanceToMergeFragmentedSPs))
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)

                    if not atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs :
                        #OK those listComplementarySPsToCheckFurther can be merged together
                        # add comment that listComplementarySPsToCheckFurther can be merged together in structure and SP
                        commentITToAdd = "{} complement each other as fragments to form a single {}. ".format(ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), SPTypeIT1IT.SPType)
                        setSPsFragmentsIT = set(listComplementarySPsToCheckFurther)
                        for SPsFragmentsIT in listComplementarySPsToCheckFurther:
                            icescreen_OO.addCommentToLocusTag2Comment(SPsFragmentsIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            SPsFragmentsIT.listSiblingFragmentedSP = sorted(listComplementarySPsToCheckFurther, key=lambda x: x.CDSPositionInGenome, reverse=False)



    def seedICEsIMEsStructure(
            self
            , locusTagIntegrase2Comment
            , SPsInSameFamilyMergeStructures2SameFamilyMergeStructure
            , listSameFamilyMergeStructures
            ):
        
        debug_seedICEsIMEsStructure = False

        def registerCurrentICEsIMEsStructureIfNeeded(
                currentICEsIMEsStructure
                , listICEsIMEsStructures
                , idxCurrentSP
                ):
            if debug_seedICEsIMEsStructure:
                print(" - registerCurrentICEsIMEsStructureIfNeeded")
            if currentICEsIMEsStructure.listOrderedSPs:
                if debug_seedICEsIMEsStructure:
                    print("\tcurrentICEsIMEsStructure.listOrderedSPs")
                currentICEsIMEsStructure.idxInSeedList = len(listICEsIMEsStructures)

                listICEsIMEsStructures.append(currentICEsIMEsStructure) # was after addPotentialUpstreamSPInConflict initially

                # EMStructure.BasicEMStructure.countInternalIdentifier += 1
                currentICEsIMEsStructure = EMStructure.ICEsIMEsStructure(True)
                return currentICEsIMEsStructure
            else:
                if debug_seedICEsIMEsStructure:
                    print("\tNOT currentICEsIMEsStructure.listOrderedSPs")
                return currentICEsIMEsStructure
    
        listICEsIMEsStructures = []
        currentICEsIMEsStructure = EMStructure.ICEsIMEsStructure(False)
        sameFamilyMergeStructureToCheck = None
        for idxCurrentSP, currentSP in enumerate(self.list):

            if debug_seedICEsIMEsStructure:
                print("** Dealing with idxCurrentSP: {}, currentSP: {}".format(str(idxCurrentSP), currentSP.locusTag))
                print(" - first step")

            greenLightAddSPConjugaisonModule = False
            transfertCommentFromSameFamilyMergeStructure = False
            if currentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:

                if debug_seedICEsIMEsStructure:
                    print("\tcurrentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure")

                sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP]
                greenLightAddSPConjugaisonModule = currentICEsIMEsStructure.listSPsIsContainedWithinOtherStructure(
                    currentSP
                    , sameFamilyMergeStructureToCheck
                    , True
                    )
                if greenLightAddSPConjugaisonModule:

                    if debug_seedICEsIMEsStructure:
                        print("\tgreenLightAddSPConjugaisonModule 1")
                    setAllowCheckingForMultipleDistantSameSPType = set()
                    greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currentICEsIMEsStructure,
                        currentSP,
                        setAllowCheckingForMultipleDistantSameSPType
                        )
                    if debug_seedICEsIMEsStructure:
                        if greenLightAddSPConjugaisonModule:
                            print("\tgreenLightAddSPConjugaisonModule 2")
                        else:
                            print("\tNOT greenLightAddSPConjugaisonModule 2")
                else:
                    if debug_seedICEsIMEsStructure:
                        print("\tNOT greenLightAddSPConjugaisonModule")

                transfertCommentFromSameFamilyMergeStructure = True
            else:

                if debug_seedICEsIMEsStructure:
                    print("\tcurrentSP NOT in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure")

                if sameFamilyMergeStructureToCheck is None:
                    if debug_seedICEsIMEsStructure:
                        print("\tsameFamilyMergeStructureToCheck is None")
                    setAllowCheckingForMultipleDistantSameSPType = set()
                    greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currentICEsIMEsStructure,
                        currentSP,
                        setAllowCheckingForMultipleDistantSameSPType
                        )
                    if debug_seedICEsIMEsStructure:
                        if greenLightAddSPConjugaisonModule:
                            print("\tgreenLightAddSPConjugaisonModule 3")
                        else:
                            print("\tNOT greenLightAddSPConjugaisonModule 3")
                else:
                    if debug_seedICEsIMEsStructure:
                        print("\tNOT sameFamilyMergeStructureToCheck is None")
                    greenLightAddSPConjugaisonModule = currentICEsIMEsStructure.listSPsIsContainedWithinOtherStructure(
                            currentSP
                            , sameFamilyMergeStructureToCheck
                            , True
                            )
                    if greenLightAddSPConjugaisonModule:
                        if debug_seedICEsIMEsStructure:
                            print("\tgreenLightAddSPConjugaisonModule 4")
                        setAllowCheckingForMultipleDistantSameSPType = set()
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                            currentICEsIMEsStructure,
                            currentSP,
                            setAllowCheckingForMultipleDistantSameSPType
                            )
                        if debug_seedICEsIMEsStructure:
                            if greenLightAddSPConjugaisonModule:
                                print("\tgreenLightAddSPConjugaisonModule 5")
                            else:
                                print("\tNOT greenLightAddSPConjugaisonModule 5")
                    else:
                        if debug_seedICEsIMEsStructure:
                            print("\tNOT greenLightAddSPConjugaisonModule 4")

            if debug_seedICEsIMEsStructure:
                print(" - second step")

            if greenLightAddSPConjugaisonModule:
                if debug_seedICEsIMEsStructure:
                    print("\tgreenLightAddSPConjugaisonModule")
                currentICEsIMEsStructure.addSPToConjugaisonModule(currentSP)
                if transfertCommentFromSameFamilyMergeStructure:
                    # transfer comment from merge structure to structure being built if yes
                    if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment not in currentICEsIMEsStructure.comment:
                        currentICEsIMEsStructure.comment += SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment
            else:
                if debug_seedICEsIMEsStructure:
                    print("\tNOT greenLightAddSPConjugaisonModule")
                # start new ICEsIMEsStructure and continue main loop
                # print("start new ICEsIMEsStructure and continue main loop")
                currentICEsIMEsStructure = registerCurrentICEsIMEsStructureIfNeeded(
                        currentICEsIMEsStructure
                        , listICEsIMEsStructures
                        , idxCurrentSP - 1
                        )
                # try adding again with new EM struct
                
                setAllowCheckingForMultipleDistantSameSPType = set()
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                    currentICEsIMEsStructure,
                    currentSP,
                    setAllowCheckingForMultipleDistantSameSPType
                    )
                if greenLightAddSPConjugaisonModule:
                    if debug_seedICEsIMEsStructure:
                        print("\tgreenLightAddSPConjugaisonModule 6")
                    currentICEsIMEsStructure.addSPToConjugaisonModule(currentSP)
                    if currentSP in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                        sameFamilyMergeStructureToCheck = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP]
                    else:
                        sameFamilyMergeStructureToCheck = None
                    if transfertCommentFromSameFamilyMergeStructure:
                        # transfer comment from merge structure to structure being built if yes
                        if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment not in currentICEsIMEsStructure.comment:
                            currentICEsIMEsStructure.comment += SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currentSP].comment
                else:
                    if debug_seedICEsIMEsStructure:
                        print("\tNOT greenLightAddSPConjugaisonModule 6")

        # last iteration
        registerCurrentICEsIMEsStructureIfNeeded(
                currentICEsIMEsStructure
                , listICEsIMEsStructures
                , len(self.list) - 1
                )

        if debug_seedICEsIMEsStructure:
            print(" - do backtracking adding SP toward upstream only after completion of the first pass toward downstream")
        # do backtracking adding SP toward upstream only after completion of the first pass toward downstream, not every time after a structure has just been stopped toward downstream
        for ICEsIMEsStructuresToCheckForConflictSPUpstream in reversed(listICEsIMEsStructures):
                ICEsIMEsStructuresToCheckForConflictSPUpstream.addPotentialUpstreamSPInConflict(
                        self.list
                        , listICEsIMEsStructures
                        , SPsInSameFamilyMergeStructures2SameFamilyMergeStructure
                        )

        if debug_seedICEsIMEsStructure:
            print(" - final listICEsIMEsStructures")
            for idx_ICEsIMEsStructuresDebug, ICEsIMEsStructuresDebug in enumerate(listICEsIMEsStructures):
                print("idx_ICEsIMEsStructuresDebug = {}: ICEsIMEsStructuresDebug = {} (list genes = {})".format(str(idx_ICEsIMEsStructuresDebug), ICEsIMEsStructuresDebug.internalIdentifier, ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructuresDebug.listOrderedSPs)))

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
                    # also consider fragment as integrase that follow up
                    if (abs(listUpstreamIntegraseToReturn[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2) or (listUpstreamIntegraseToReturn[-1] in currSp.listSiblingFragmentedSP and currSp in listUpstreamIntegraseToReturn[-1].listSiblingFragmentedSP):  # Integrases are adjacent in genome or separated by a CDS
                        if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, True)
                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                # only add comment for adjacent intergase
                                # rulesAddIntegrases.addCommentOnIntegraseNotCorrectlyOrientedForICEsIMEsStructure(currSp, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment, "upstream", "-", valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure)
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
                            rulesAddIntegrases.addCommentOnIntegraseNotCorrectlyOrientedForICEsIMEsStructure(currSp, downstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment, "upstream", "-", valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure, "getListIntegraseGroupJustUpstreamOfThisSP")
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
                    # also consider fragment as integrase that follow up
                    if (abs(listDownstreamIntegraseToReturn[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2) or (listDownstreamIntegraseToReturn[-1] in currSp.listSiblingFragmentedSP and currSp in listDownstreamIntegraseToReturn[-1].listSiblingFragmentedSP):  # Integrases are adjacent in genome or separated by a CDS
                        if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                            if upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented is not None :
                                valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = rulesAddIntegrases.isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, False)
                                if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                    # only add comment for adjacent intergase
                                    # rulesAddIntegrases.addCommentOnIntegraseNotCorrectlyOrientedForICEsIMEsStructure(currSp, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment, "downstream", "+", valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure)
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
                                rulesAddIntegrases.addCommentOnIntegraseNotCorrectlyOrientedForICEsIMEsStructure(currSp, upstreamICEsIMEsStructureToTestForIntegraseCorrectlyOriented, locusTagIntegrase2Comment, "downstream", "+", valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure, "getListIntegraseGroupJustUpstreamOfThisSP")
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

                if (currSP.SPType == "DDE transposase"):
                    strCurrSP = "D"
                elif (currSP.SPType in icescreen_OO.setIntegraseNames):
                    strCurrSP = "I"
                elif (currSP.SPType in icescreen_OO.listTypeSPConjModule):
                    strCurrSP = currSP.SPType[0]
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
                elif (currSP.SPType in icescreen_OO.listTypeSPConjModule):
                    strCurrSP = currSP.SPType[0]
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
