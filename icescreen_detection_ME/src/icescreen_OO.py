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

import argparse
import configparser
import os
# import sys
import time
import datetime
from datetime import timedelta
import csv
# import re
# import specific class OO for this script
import hit
import EMStructure
import rulesAddIntegrases
from EMTypeStructure import printOverallStatsToSummaryFile
import rulesMergeICEIMEStructures
import commonMethods

########
# GLOBAL VARS #
########
#scriptVersion = "1.4"
# setIntegraseTyrNames = {"IntTyr", "TyrInt"}
setIntegraseTyrNames = {"Tyrosine integrase"}
# setIntegraseSerNames = {"IntSer", "SerInt"}
setIntegraseSerNames = {"Serine integrase"}
setIntegraseDDENames = {"DDE transposase"}
setIntegraseNames = set()
setIntegraseNames.update(setIntegraseTyrNames, setIntegraseSerNames, setIntegraseDDENames)
# regexSplitBlastFamilyIntoTokens = '-|\||;|; |\n|, |,' # not multi familly in file detectedSP yet
segmentIdx2startGenomicRegion = {}
segmentIdx2stopGenomicRegion = {}

########
# GLOBAL METHODS #
########


# locusTag2Comment is a dictionary used to store the comments that will be visible in the output files for each SPs.
def addCommentToLocusTag2Comment(locusTagSent, commentSent, locusTag2Comment):
    if locusTagSent in locusTag2Comment:  # key already there
        commentIT = locusTag2Comment[locusTagSent]
        if commentSent not in commentIT:
            commentIT += commentSent
            locusTag2Comment[locusTagSent] = commentIT
    else:  # key not there
        locusTag2Comment[locusTagSent] = commentSent


def removeCommentToLocusTag2Comment(locusTagSent, commentSent, locusTag2Comment):
    if locusTagSent in locusTag2Comment:  # key already there
        commentIT = locusTag2Comment[locusTagSent]
        if commentSent in commentIT:
            locusTag2Comment[locusTagSent] = commentIT.replace(commentSent, "")
    else:  # key not there
        pass


# This method generates the content of the output csv file: input cvs file + information on the ICE or IME structure.
def printAllICEsIMEsStructureToInputFile(
        listOfListAllICEsIMEsStructure,
        listOfListSPsLonelyIntegrases,
        pathInputFile,
        modifiedInputFile,
        locusTagIntegrase2Comment,
        locusTagFinalize2Comment,
        locusTagMerge2Comment,
        dictIMEICEID2humanReadableIMEICEIIdentifier,
        hasMultipleGenomeAccesion
        ):

    printSegmentNumber = True

    # totalNumberSP = 0
    # totalNumberIntegrase = 0
    # totalNumberUnaffectedIntegrase = 0
    # totalNumberRelaxase = 0
    # totalNumberUnaffectedRelaxase = 0
    # totalNumberCoupling = 0
    # totalNumberUnaffectedCoupling = 0
    # totalNumberVirb4 = 0
    # totalNumberUnaffectedVirb4 = 0
    dictGenomeAccnum2totalNumberSP = {}
    dictGenomeAccnum2totalNumberIntegrase = {}
    dictGenomeAccnum2totalNumberUnaffectedIntegrase = {}
    dictGenomeAccnum2totalNumberRelaxase = {}
    dictGenomeAccnum2totalNumberUnaffectedRelaxase = {}
    dictGenomeAccnum2totalNumberCoupling = {}
    dictGenomeAccnum2totalNumberUnaffectedCoupling = {}
    dictGenomeAccnum2totalNumberVirb4 = {}
    dictGenomeAccnum2totalNumberUnaffectedVirb4 = {}
    dictGenomeAccnum2totalsetFragmentedSP = {}
    dictLocusTagSure2ICEIMEInternalId = {}
    dictLocusTagNotSure2ICEIMEInternalId = {}
    dictLocusTag2segmentNumber = {}

    for idx, currListAllICEsIMEsStructure in enumerate(listOfListAllICEsIMEsStructure):
        for currICEIMEStructure in currListAllICEsIMEsStructure:
            for currSP in currICEIMEStructure.listOrderedSPs:
                dictLocusTag2segmentNumber[currSP.locusTag] = idx + 1
                if currSP.locusTag in dictLocusTagSure2ICEIMEInternalId:
                    setInternalId = dictLocusTagSure2ICEIMEInternalId[currSP.locusTag]
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                else:
                    setInternalId = set()
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                    dictLocusTagSure2ICEIMEInternalId[currSP.locusTag] = setInternalId
                if currSP.genomeAccession in dictGenomeAccnum2totalsetFragmentedSP:
                    setFragmentedSPIT = dictGenomeAccnum2totalsetFragmentedSP[currSP.genomeAccession]
                    if len(currSP.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSP)
                else :
                    setFragmentedSPIT = set()
                    if len(currSP.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSP)
                    dictGenomeAccnum2totalsetFragmentedSP[currSP.genomeAccession] = setFragmentedSPIT

            for currSPConjModuleToManuallyCheckIT in currICEIMEStructure.setSPConjModuleToManuallyCheck:
                dictLocusTag2segmentNumber[currSPConjModuleToManuallyCheckIT.locusTag] = idx + 1
                if currSPConjModuleToManuallyCheckIT.locusTag in dictLocusTagNotSure2ICEIMEInternalId:
                    setInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPConjModuleToManuallyCheckIT.locusTag]
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                else:
                    setInternalId = set()
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                    dictLocusTagNotSure2ICEIMEInternalId[currSPConjModuleToManuallyCheckIT.locusTag] = setInternalId
                if currSPConjModuleToManuallyCheckIT.genomeAccession in dictGenomeAccnum2totalsetFragmentedSP:
                    setFragmentedSPIT = dictGenomeAccnum2totalsetFragmentedSP[currSPConjModuleToManuallyCheckIT.genomeAccession]
                    if len(currSPConjModuleToManuallyCheckIT.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSPConjModuleToManuallyCheckIT)
                else :
                    setFragmentedSPIT = set()
                    if len(currSPConjModuleToManuallyCheckIT.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSPConjModuleToManuallyCheckIT)
                    dictGenomeAccnum2totalsetFragmentedSP[currSPConjModuleToManuallyCheckIT.genomeAccession] = setFragmentedSPIT

                    
            for currSPIntegraseToManuallyCheckIT in currICEIMEStructure.setIntegraseToManuallyCheck:
                dictLocusTag2segmentNumber[currSPIntegraseToManuallyCheckIT.locusTag] = idx + 1
                if currSPIntegraseToManuallyCheckIT.locusTag in dictLocusTagNotSure2ICEIMEInternalId:
                    setInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPIntegraseToManuallyCheckIT.locusTag]
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                else:
                    setInternalId = set()
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                    dictLocusTagNotSure2ICEIMEInternalId[currSPIntegraseToManuallyCheckIT.locusTag] = setInternalId
                if currSPIntegraseToManuallyCheckIT.genomeAccession in dictGenomeAccnum2totalsetFragmentedSP:
                    setFragmentedSPIT = dictGenomeAccnum2totalsetFragmentedSP[currSPIntegraseToManuallyCheckIT.genomeAccession]
                    if len(currSPIntegraseToManuallyCheckIT.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSPIntegraseToManuallyCheckIT)
                else :
                    setFragmentedSPIT = set()
                    if len(currSPIntegraseToManuallyCheckIT.listSiblingFragmentedSP) > 0:
                        setFragmentedSPIT.add(currSPIntegraseToManuallyCheckIT)
                    dictGenomeAccnum2totalsetFragmentedSP[currSPIntegraseToManuallyCheckIT.genomeAccession] = setFragmentedSPIT



    for idx, currListSPsLonelyIntegrases in enumerate(listOfListSPsLonelyIntegrases):
        for currSPLonelyIntegrases in currListSPsLonelyIntegrases:
            dictLocusTag2segmentNumber[currSPLonelyIntegrases.locusTag] = idx + 1

    with open(pathInputFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        countIterRow = 0
        #idxColCDS = -1
        idxColCDS_locus_tag = -1
        idxColCDS_protein_id = -1
        idxColGenome_accession = -1
        idxColCDS_start = -1
        idxColCDS_Protein_Type = -1
        for row in reader:
            if countIterRow == 0:
                # check header
                countIterCol = 0
                for column in row:
                    # if column == "CDS":
                    #     idxColCDS = countIterCol
                    if column == "CDS_locus_tag":
                        idxColCDS_locus_tag = countIterCol
                    elif column == "CDS_protein_id":
                        idxColCDS_protein_id = countIterCol
                    elif column == "Genome_accession":
                        idxColGenome_accession = countIterCol
                    elif column == "CDS_start":
                        idxColCDS_start = countIterCol
                    elif column == "CDS_Protein_type":
                        idxColCDS_Protein_Type = countIterCol
                    countIterCol += 1
                if printSegmentNumber:
                    print("ICE_IME_id\tICE_IME_id_need_manual_curation\tSegment_number\tComments_ICE_IME_structure\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)
                else:
                    print("ICE_IME_id\tICE_IME_id_need_manual_curation\tComments_ICE_IME_structure\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)

            else:
                currSPLocusTag = ""
                currSPProteinId = ""
                currSPGenomeAccession = ""
                currSPstart = ""
                currSPType = ""  # Coupling, Relaxase, Virb4, or integrase
                #geneHasRealLocusTag = False
                # if idxColCDS >= 0:
                #     parsedCell = row[idxColCDS]
                #     currSPProteinId = parsedCell
                # else:
                #     raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS\"')
                if idxColCDS_locus_tag >= 0:
                    parsedCell = row[idxColCDS_locus_tag]
                    if len(parsedCell) > 0:
                        currSPLocusTag = parsedCell
                        #geneHasRealLocusTag = True
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS_locus_tag\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_protein_id >= 0:
                    parsedCell = row[idxColCDS_protein_id]
                    if len(parsedCell) > 0:
                        currSPProteinId = parsedCell
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS_protein_id\" in file {}'.format(str(pathInputFile)))
                if idxColGenome_accession >= 0:
                    parsedCell = row[idxColGenome_accession]
                    if len(parsedCell) > 0:
                        currSPGenomeAccession = parsedCell
                    else:
                        raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory information on Genome_accession for row {} in file {}'.format(str(row), str(pathInputFile)))
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"Genome_accession\" in file {}'.format(str(pathInputFile)))

                if idxColCDS_start >= 0:
                    parsedCell = row[idxColCDS_start]
                    if len(parsedCell) > 0:
                        parsedCell = int(parsedCell)
                        currSPstart = parsedCell
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS_start\"')
                
                #currSPLocusTag = currSPProteinId + "-" + str(currSPstart)
                #if not geneHasRealLocusTag:
                currSPLocusTag = commonMethods.makeCompositeUniqLocusTag(hasMultipleGenomeAccesion, currSPLocusTag, currSPProteinId, currSPGenomeAccession, currSPstart)

                if idxColCDS_Protein_Type >= 0:
                    parsedCell = row[idxColCDS_Protein_Type]
                    if (parsedCell == "Coupling protein" or parsedCell == "Relaxase" or parsedCell == "VirB4" or parsedCell in setIntegraseNames):  # parsedCell == "IntTyr" or parsedCell == "IntSer" or parsedCell == "DDE"
                        currSPType = parsedCell
                    else:
                        raise RuntimeError(
                                'printAllICEsIMEsStructureToInputFile error: could not parse value \"{}\" of column \"CDS_Protein_type\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS_Protein_type\"')

                strSure = ""
                if currSPLocusTag in dictLocusTagSure2ICEIMEInternalId:
                    currSetInternalId = dictLocusTagSure2ICEIMEInternalId[currSPLocusTag]
                    strSure = " ".join(str(i) for i in sorted(currSetInternalId))
                for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
                    strSure = strSure.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)


                strNotSure = ""
                if currSPLocusTag in dictLocusTagNotSure2ICEIMEInternalId:
                    currSetInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTag]
                    strNotSure = " ".join(str(i) for i in sorted(currSetInternalId))
                for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
                    strNotSure = strNotSure.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)

                strLocusTag2CommentIT = ""
                if currSPLocusTag in locusTagMerge2Comment:
                    if locusTagMerge2Comment[currSPLocusTag] not in strLocusTag2CommentIT:
                        strLocusTag2CommentIT += locusTagMerge2Comment[currSPLocusTag]
                if currSPLocusTag in locusTagIntegrase2Comment:
                    if locusTagIntegrase2Comment[currSPLocusTag] not in strLocusTag2CommentIT:
                        strLocusTag2CommentIT += locusTagIntegrase2Comment[currSPLocusTag]
                if currSPLocusTag in locusTagFinalize2Comment:
                    if locusTagFinalize2Comment[currSPLocusTag] not in strLocusTag2CommentIT:
                        strLocusTag2CommentIT += locusTagFinalize2Comment[currSPLocusTag]
                for IMEICEIDIT, humanReadableIMEICEIIdentifierIT in dictIMEICEID2humanReadableIMEICEIIdentifier.items():
                    strLocusTag2CommentIT = strLocusTag2CommentIT.replace(IMEICEIDIT, humanReadableIMEICEIIdentifierIT)

                if printSegmentNumber:
                    strSegmentNumber = ""
                    if currSPLocusTag in dictLocusTag2segmentNumber:
                        strSegmentNumber = dictLocusTag2segmentNumber[currSPLocusTag]
                    print(strSure + "\t" + strNotSure + "\t" + str(strSegmentNumber) + "\t" + strLocusTag2CommentIT + "\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)
                else:
                    print(strSure + "\t" + strNotSure + "\t" + strLocusTag2CommentIT + "\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)


                #totalNumberSP += 1
                if currSPGenomeAccession in dictGenomeAccnum2totalNumberSP:
                    dictGenomeAccnum2totalNumberSP[currSPGenomeAccession] += 1
                else :
                    dictGenomeAccnum2totalNumberSP[currSPGenomeAccession] = 1
                if currSPType == "Coupling protein":
                    #totalNumberCoupling += 1
                    if currSPGenomeAccession in dictGenomeAccnum2totalNumberCoupling:
                        dictGenomeAccnum2totalNumberCoupling[currSPGenomeAccession] += 1
                    else :
                        dictGenomeAccnum2totalNumberCoupling[currSPGenomeAccession] = 1
                    if not strSure:
                        #totalNumberUnaffectedCoupling += 1
                        if currSPGenomeAccession in dictGenomeAccnum2totalNumberUnaffectedCoupling:
                            dictGenomeAccnum2totalNumberUnaffectedCoupling[currSPGenomeAccession] += 1
                        else :
                            dictGenomeAccnum2totalNumberUnaffectedCoupling[currSPGenomeAccession] = 1
                elif currSPType == "Relaxase":
                    #totalNumberRelaxase += 1
                    if currSPGenomeAccession in dictGenomeAccnum2totalNumberRelaxase:
                        dictGenomeAccnum2totalNumberRelaxase[currSPGenomeAccession] += 1
                    else:
                        dictGenomeAccnum2totalNumberRelaxase[currSPGenomeAccession] = 1
                    if not strSure:
                        #totalNumberUnaffectedRelaxase += 1
                        if currSPGenomeAccession in dictGenomeAccnum2totalNumberUnaffectedRelaxase:
                            dictGenomeAccnum2totalNumberUnaffectedRelaxase[currSPGenomeAccession] += 1
                        else:
                            dictGenomeAccnum2totalNumberUnaffectedRelaxase[currSPGenomeAccession] = 1
                elif currSPType == "VirB4":
                    #totalNumberVirb4 += 1
                    if currSPGenomeAccession in dictGenomeAccnum2totalNumberVirb4:
                        dictGenomeAccnum2totalNumberVirb4[currSPGenomeAccession] += 1
                    else:
                        dictGenomeAccnum2totalNumberVirb4[currSPGenomeAccession] = 1
                    if not strSure:
                        #totalNumberUnaffectedVirb4 += 1
                        if currSPGenomeAccession in dictGenomeAccnum2totalNumberUnaffectedVirb4:
                            dictGenomeAccnum2totalNumberUnaffectedVirb4[currSPGenomeAccession] += 1
                        else:
                            dictGenomeAccnum2totalNumberUnaffectedVirb4[currSPGenomeAccession] = 1
                elif currSPType in setIntegraseNames:
                    #totalNumberIntegrase += 1
                    if currSPGenomeAccession in dictGenomeAccnum2totalNumberIntegrase:
                        dictGenomeAccnum2totalNumberIntegrase[currSPGenomeAccession] += 1
                    else:
                        dictGenomeAccnum2totalNumberIntegrase[currSPGenomeAccession] = 1
                    if not strSure:
                        #totalNumberUnaffectedIntegrase += 1
                        if currSPGenomeAccession in dictGenomeAccnum2totalNumberUnaffectedIntegrase:
                            dictGenomeAccnum2totalNumberUnaffectedIntegrase[currSPGenomeAccession] += 1
                        else:
                            dictGenomeAccnum2totalNumberUnaffectedIntegrase[currSPGenomeAccession] = 1
                else:
                    raise RuntimeError(
                            'printAllICEsIMEsStructureToInputFile error: unrecognized SPType {} for LocusTag {}'.format(
                                    currSPType, currSPLocusTag))
                
            countIterRow += 1

    csvfile.close()

    # return (totalNumberSP,
    #         totalNumberIntegrase,
    #         totalNumberUnaffectedIntegrase,
    #         totalNumberRelaxase,
    #         totalNumberUnaffectedRelaxase,
    #         totalNumberCoupling,
    #         totalNumberUnaffectedCoupling,
    #         totalNumberVirb4,
    #         totalNumberUnaffectedVirb4)
    return (dictGenomeAccnum2totalNumberSP,
        dictGenomeAccnum2totalNumberIntegrase,
        dictGenomeAccnum2totalNumberUnaffectedIntegrase,
        dictGenomeAccnum2totalNumberRelaxase,
        dictGenomeAccnum2totalNumberUnaffectedRelaxase,
        dictGenomeAccnum2totalNumberCoupling,
        dictGenomeAccnum2totalNumberUnaffectedCoupling,
        dictGenomeAccnum2totalNumberVirb4,
        dictGenomeAccnum2totalNumberUnaffectedVirb4,
        dictGenomeAccnum2totalsetFragmentedSP
        )




# This method parse the input file into objects used internally in the program.
#def parse_csv(outputPrintCDSNumberInsteadOfProteinIdAndStart, pathInputFile):
def parse_csv(pathInputFile, hasMultipleGenomeAccesion):
    listSPsParsed = hit.ListSPs()
    with open(pathInputFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        countIterRow = 0
        #idxColCDS = -1
        idxColCDS_locus_tag = -1
        idxColCDS_protein_id = -1
        idxColGenome_accession = -1
        idxColGenome_accession_rank = -1
        idxColCDS_start = -1
        idxColCDS_end = -1
        idxColCDS_strand = -1
        idxColCDS_num = -1
        idxColCDS_Protein_Type = -1
        idxColIs_hit_blast = -1
        idxColIs_hit_HMM = -1
        idxColIs_pseudo = -1
        idxColLength_of_blast_most_similar_ref_SP = -1
        idxColBlast_ali_length = -1
        idxColBlast_ali_start_CDS = -1
        idxColBlast_ali_end_CDS = -1
        idxColBlast_ali_start_Query_blast = -1
        idxColBlast_ali_end_Query_blast = -1
        idxColBlast_ali_identity_perc = -1
        idxColE_value_blast = -1
        idxColBlast_ali_bitscore = -1
        idxColCDS_coverage_blast = -1
        idxColBlast_ali_coverage_most_similar_ref_SP = -1
        # idxColIs_hit_HMM_CC = -1
        # idxColElement_family = -1
        idxColICE_superfamily_of_most_similar_ref_SP = -1
        idxColICE_family_of_most_similar_ref_SP = -1
        idxColIME_superfamily_of_most_similar_ref_SP = -1
        idxColRelaxase_family_domain_of_most_similar_ref_SP = -1
        idxColRelaxase_family_MOB_of_most_similar_ref_SP = -1
        idxColCoupling_type_of_most_similar_ref_SP = -1
        idxColFalse_positives = -1
        idxColSP_blast_validation = -1
        idxColUse_annotation = -1
        idxColBest_hmmprofile = -1
        # idxColBest_hmmprofile_CC = -1
        dictCheckUniquLocusTag = {}
        for row in reader:
            if countIterRow == 0:
                # check header
                countIterCol = 0
                for column in row:
                    # if column == "CDS":
                    #     idxColCDS = countIterCol
                    if column == "CDS_locus_tag":
                        idxColCDS_locus_tag = countIterCol
                    elif column == "CDS_protein_id":
                        idxColCDS_protein_id = countIterCol
                    elif column == "Genome_accession":
                        idxColGenome_accession = countIterCol
                    elif column == "Genome_accession_rank":
                        idxColGenome_accession_rank = countIterCol
                    elif column == "CDS_start":
                        idxColCDS_start = countIterCol
                    elif column == "CDS_end":
                        idxColCDS_end = countIterCol
                    elif column == "CDS_strand":
                        idxColCDS_strand = countIterCol
                    elif column == "CDS_num":
                        idxColCDS_num = countIterCol
                    elif column == "CDS_Protein_type":
                        idxColCDS_Protein_Type = countIterCol
                    elif column == "Is_hit_blast":
                        idxColIs_hit_blast = countIterCol
                    elif column == "Is_hit_HMM":
                        idxColIs_hit_HMM = countIterCol
                    elif column == "Is_pseudo":
                        idxColIs_pseudo = countIterCol
                    elif column == "Length_of_blast_most_similar_ref_SP":
                        idxColLength_of_blast_most_similar_ref_SP = countIterCol
                    elif column == "Blast_ali_length":
                        idxColBlast_ali_length = countIterCol
                    elif column == "Blast_ali_start_CDS":
                        idxColBlast_ali_start_CDS = countIterCol
                    elif column == "Blast_ali_end_CDS":
                        idxColBlast_ali_end_CDS = countIterCol
                    elif column == "Blast_ali_start_Query_blast":
                        idxColBlast_ali_start_Query_blast = countIterCol
                    elif column == "Blast_ali_end_Query_blast":
                        idxColBlast_ali_end_Query_blast = countIterCol
                    elif column == "Blast_ali_identity_perc":
                        idxColBlast_ali_identity_perc = countIterCol
                    elif column == "Blast_ali_E-value":
                        idxColE_value_blast = countIterCol
                    elif column == "Blast_ali_bitscore":
                        idxColBlast_ali_bitscore = countIterCol
                    elif column == "CDS_coverage_blast":
                        idxColCDS_coverage_blast = countIterCol
                    elif column == "Blast_ali_coverage_most_similar_ref_SP":
                        idxColBlast_ali_coverage_most_similar_ref_SP = countIterCol
                    # elif column=="Is_hit_HMM_CC":
                    #    idxColIs_hit_HMM_CC = countIterCol
                    # elif column=="Element_family":
                    #    idxColElement_family = countIterCol
                    elif column == "ICE_superfamily_of_most_similar_ref_SP":
                        idxColICE_superfamily_of_most_similar_ref_SP = countIterCol
                    elif column == "ICE_family_of_most_similar_ref_SP":
                        idxColICE_family_of_most_similar_ref_SP = countIterCol
                    elif column == "IME_superfamily_of_most_similar_ref_SP":
                        idxColIME_superfamily_of_most_similar_ref_SP = countIterCol
                    elif column == "Relaxase_family_domain_of_most_similar_ref_SP":
                        idxColRelaxase_family_domain_of_most_similar_ref_SP = countIterCol
                    elif column == "Relaxase_family_MOB_of_most_similar_ref_SP":
                        idxColRelaxase_family_MOB_of_most_similar_ref_SP = countIterCol
                    elif column == "Coupling_type_of_most_similar_ref_SP":
                        idxColCoupling_type_of_most_similar_ref_SP = countIterCol
                    elif column == "False_positives":
                        idxColFalse_positives = countIterCol
                    elif column == "SP_blast_validation":
                        idxColSP_blast_validation = countIterCol
                    elif column == "Use_annotation":
                        idxColUse_annotation = countIterCol
                    elif column == "Description_of_matching_HMM_profile":
                        # was elif column=="Best_hmmprofile":
                        idxColBest_hmmprofile = countIterCol
                    # elif column=="Best_hmmprofile_CC":
                    #    idxColBest_hmmprofile_CC = countIterCol
                    countIterCol += 1
                if ((idxColIs_hit_blast + idxColIs_hit_HMM) == -2):  # + idxColIs_hit_HMM_CC
                    raise RuntimeError('Input file error: absence of at least one of the mandatory columns \"Is_hit_blast\" or \"Is_hit_HMM\"')
            else:
                currSP = hit.SP()
                #geneHasRealLocusTag = False
                # if idxColCDS >= 0:
                #     parsedCell = row[idxColCDS]
                #     currSP.proteinId = parsedCell
                #     if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                #         currSP.locusTag = currSP.proteinId + "-" + str(currSP.start)
                #     elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                #         pass
                #     else:
                #         raise RuntimeError(
                #                 "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                #                         outputPrintCDSNumberInsteadOfProteinIdAndStart))
                # else:
                #     raise RuntimeError('Input file error: missing mandatory column \"CDS\"')
                locusTagParsed = ""
                if idxColCDS_locus_tag >= 0:
                    parsedCell = row[idxColCDS_locus_tag]
                    if len(parsedCell) > 0:
                        locusTagParsed = parsedCell
                        #geneHasRealLocusTag = True
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_locus_tag\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_protein_id >= 0:
                    parsedCell = row[idxColCDS_protein_id]
                    if len(parsedCell) > 0:
                        currSP.proteinId = parsedCell
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_protein_id\" in file {}'.format(str(pathInputFile)))
                if idxColGenome_accession >= 0:
                    parsedCell = row[idxColGenome_accession]
                    if len(parsedCell) > 0:
                        currSP.genomeAccession = parsedCell
                    else:
                        raise RuntimeError('Input file error: missing mandatory information on Genome_accession for row {} in file {}'.format(str(row), str(pathInputFile)))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Genome_accession\" in file {}'.format(str(pathInputFile)))

                if idxColGenome_accession_rank >= 0:
                    parsedCell = row[idxColGenome_accession_rank]
                    if len(parsedCell) > 0:
                        parsedCell = int(parsedCell)
                        currSP.genomeAccessionRank = parsedCell
                    else:
                        raise RuntimeError('Input file error: missing mandatory information on Genome_accession_rank for row {} in file {}'.format(str(row), str(pathInputFile)))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Genome_accession_rank\" in file {}'.format(str(pathInputFile)))

                if idxColCDS_start >= 0:
                    parsedCell = row[idxColCDS_start]
                    if len(parsedCell) > 0:
                        parsedCell = int(parsedCell)
                        currSP.start = parsedCell
                    # if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                    #     currSP.locusTag = currSP.proteinId + "-" + str(currSP.start)
                    # elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                    #     pass
                    # else:
                    #     raise RuntimeError(
                    #             "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                    #                     outputPrintCDSNumberInsteadOfProteinIdAndStart))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_start\" in file {}'.format(str(pathInputFile)))
                
                #if not geneHasRealLocusTag:
                currSP.locusTag = commonMethods.makeCompositeUniqLocusTag(hasMultipleGenomeAccesion, locusTagParsed, currSP.proteinId, currSP.genomeAccession, currSP.start)
                if len(currSP.locusTag) == 0:
                    raise RuntimeError('Input file error: could not determine locusTag for row\n{}\n in file {}'.format(str(row), str(pathInputFile)))
                if currSP.locusTag in dictCheckUniquLocusTag:
                    raise RuntimeError('Input file error: duplicate locusTag {}'.format(currSP.locusTag))
                else:
                    dictCheckUniquLocusTag[currSP.locusTag] = 1

                if idxColCDS_end >= 0:
                    parsedCell = row[idxColCDS_end]
                    parsedCell = int(parsedCell)
                    currSP.stop = parsedCell
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_end\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_strand >= 0:
                    parsedCell = row[idxColCDS_strand]
                    if (parsedCell == "+" or parsedCell == "-"):
                        currSP.strand = parsedCell
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"CDS_strand\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_strand\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_num >= 0:
                    parsedCell = row[idxColCDS_num]
                    parsedCell = int(parsedCell)
                    currSP.CDSPositionInGenome = parsedCell
                    # if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                    #     pass
                    # elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                    #     currSP.locusTag = str(currSP.CDSPositionInGenome)
                    # else:
                    #     raise RuntimeError(
                    #             "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                    #                     outputPrintCDSNumberInsteadOfProteinIdAndStart))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_num\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_Protein_Type >= 0:
                    parsedCell = row[idxColCDS_Protein_Type]
                    if (parsedCell == "Coupling protein" or parsedCell == "Relaxase" or parsedCell == "VirB4" or parsedCell in setIntegraseNames):  # parsedCell == "IntTyr" or parsedCell == "IntSer" or parsedCell == "DDE"
                        currSP.SPType = parsedCell
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"CDS_Protein_type\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_Protein_type\" in file {}'.format(str(pathInputFile)))
                if idxColIs_hit_blast >= 0:
                    parsedCell = row[idxColIs_hit_blast]
                    parsedCell = int(parsedCell[0])  # [0] to get only the first character and avoid issues with number formatting according to different locale
                    if (parsedCell == 0 or parsedCell == 1):
                        currSP.SPDetectedByBlast = parsedCell
                        # if ( parsedCell == 1 and idxColElement_family >= 0 ):
                        #    parsedCell = row[idxColElement_family]
                        #    listFamilyTokens = re.split(regexSplitBlastFamilyIntoTokens, parsedCell)
                        #    for currFamilyToken in listFamilyTokens:
                        #        if currFamilyToken:# add only if not empty or null
                        #            currSP.setSPFamilyFromBlast.add(currFamilyToken)
                        if parsedCell == 1:
                            if idxColICE_superfamily_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColICE_superfamily_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPICESuperFamilyFromBlast.add(parsedCell)
                            if idxColICE_family_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColICE_family_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPICEFamilyFromBlast.add(parsedCell)
                            if idxColIME_superfamily_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColIME_superfamily_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPIMESuperFamilyFromBlast.add(parsedCell)
                            if idxColRelaxase_family_domain_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColRelaxase_family_domain_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    if currSP.SPType == "Relaxase":  # col CDS_Protein_Type ; Values: Coupling, Relaxase, Virb4
                                        currSP.Relaxase_family_domain_of_most_similar_ref_SPFromBlast = parsedCell
                                        #print("{} have Relaxase_family_domain_of_most_similar_ref_SPFromBlast \"{}\"".format(currSP.locusTag, currSP.Relaxase_family_domain_of_most_similar_ref_SPFromBlast))
                                    else:
                                        raise RuntimeError('Input file error: adding Relaxase_family_domain_of_most_similar_ref_SP to SP {} ({}-{}) of Type {}'.format(currSP.locusTag, currSP.start, currSP.stop, currSP.SPType))
                            if idxColRelaxase_family_MOB_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColRelaxase_family_MOB_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    if currSP.SPType == "Relaxase":  # col CDS_Protein_Type ; Values: Coupling, Relaxase, Virb4
                                        currSP.Relaxase_family_MOB_of_most_similar_ref_SPFromBlast = parsedCell
                                    else:
                                        raise RuntimeError('Input file error: adding Relaxase_family_MOB_of_most_similar_ref_SP to SP {} ({}-{}) of Type {}'.format(currSP.locusTag, currSP.start, currSP.stop, currSP.SPType))
                            if idxColCoupling_type_of_most_similar_ref_SP >= 0:
                                parsedCell = row[idxColCoupling_type_of_most_similar_ref_SP]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    if currSP.SPType == "Coupling protein":  # col CDS_Protein_Type ; Values: Coupling, Relaxase, Virb4
                                        currSP.Coupling_type_of_most_similar_ref_SPFromBlast = parsedCell
                                    else:
                                        raise RuntimeError('Input file error: adding Coupling_type_of_most_similar_ref_SP to SP {} ({}-{}) of Type {}'.format(currSP.locusTag, currSP.start, currSP.stop, currSP.SPType))
                            if idxColFalse_positives >= 0:
                                parsedCell = row[idxColFalse_positives]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.False_positivesFromBlast = parsedCell
                            if idxColSP_blast_validation >= 0:
                                parsedCell = row[idxColSP_blast_validation]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.blast_validation = parsedCell
                            if idxColUse_annotation >= 0:
                                parsedCell = row[idxColUse_annotation]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.Use_annotationFromBlast = parsedCell
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"Is_hit_blast\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                if idxColIs_hit_HMM >= 0:
                    parsedCell = row[idxColIs_hit_HMM]
                    parsedCell = int(parsedCell[0])  # [0] to get only the first character and avoid issues with number formatting according to different locale
                    if (parsedCell == 0 or parsedCell == 1):
                        currSP.SPDetectedByHMM = parsedCell
                        if (parsedCell == 1 and idxColBest_hmmprofile >= 0):
                            parsedCell = row[idxColBest_hmmprofile]
                            # listFamilyTokens = re.split(regexSplitBlastFamilyIntoTokens, parsedCell)
                            # for currFamilyToken in listFamilyTokens:
                            #    if currFamilyToken:# add only if not empty or null
                            #        currSP.setSPFamilyFromHMM.add(currFamilyToken)
                            if parsedCell:  # add only if not empty or null
                                currSP.setSPFamilyFromHMM.add(parsedCell)
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"Is_hit_HMM\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                if idxColIs_pseudo >= 0:
                    parsedCell = row[idxColIs_pseudo]
                    if parsedCell == "Pseudo" or parsedCell == "pseudo" or parsedCell == "True" or parsedCell == "true" :
                        currSP.pseudo = True
                if idxColLength_of_blast_most_similar_ref_SP >= 0:
                    parsedCell = row[idxColLength_of_blast_most_similar_ref_SP]
                    try:
                        currSP.Length_of_blast_most_similar_ref_SP = int(float(parsedCell)) # 796
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Length_of_blast_most_similar_ref_SP\" in file {}'.format(str(pathInputFile)))
                if idxColBlast_ali_length >= 0:
                    parsedCell = row[idxColBlast_ali_length]
                    try:
                        currSP.Blast_ali_length = int(float(parsedCell)) # 786
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Blast_ali_length\" in file {}'.format(str(pathInputFile)))
                if idxColBlast_ali_start_CDS >= 0:
                    parsedCell = row[idxColBlast_ali_start_CDS]
                    try:
                        currSP.Blast_ali_start_CDS = int(float(parsedCell))
                    except:
                        pass
                if idxColBlast_ali_end_CDS >= 0:
                    parsedCell = row[idxColBlast_ali_end_CDS]
                    try:
                        currSP.Blast_ali_end_CDS = int(float(parsedCell))
                    except:
                        pass
                if idxColBlast_ali_start_Query_blast >= 0:
                    parsedCell = row[idxColBlast_ali_start_Query_blast]
                    try:
                        currSP.Blast_ali_start_Query_blast = int(float(parsedCell))
                    except:
                        pass
                if idxColBlast_ali_end_Query_blast >= 0:
                    parsedCell = row[idxColBlast_ali_end_Query_blast]
                    try:
                        currSP.Blast_ali_end_Query_blast = int(float(parsedCell))
                    except:
                        pass
                if idxColBlast_ali_identity_perc >= 0:
                    parsedCell = row[idxColBlast_ali_identity_perc]
                    try:
                        currSP.Blast_ali_identity_perc = float(parsedCell) # 100.0
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Blast_ali_identity_perc\" in file {}'.format(str(pathInputFile)))
                if idxColE_value_blast >= 0:
                    parsedCell = row[idxColE_value_blast]
                    try:
                        currSP.E_value_blast = float(parsedCell)
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Blast_ali_E-value\" in file {}'.format(str(pathInputFile)))
                if idxColBlast_ali_bitscore >= 0:
                    parsedCell = row[idxColBlast_ali_bitscore]
                    try:
                        currSP.Blast_ali_bitscore = float(parsedCell)
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Blast_ali_bitscore\" in file {}'.format(str(pathInputFile)))
                if idxColCDS_coverage_blast >= 0:
                    parsedCell = row[idxColCDS_coverage_blast]
                    try:
                        currSP.CDS_coverage_blast = float(parsedCell)
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_coverage_blast\" in file {}'.format(str(pathInputFile)))
                if idxColBlast_ali_coverage_most_similar_ref_SP >= 0:
                    parsedCell = row[idxColBlast_ali_coverage_most_similar_ref_SP]
                    try:
                        currSP.Blast_ali_coverage_most_similar_ref_SP = float(parsedCell)
                    except:
                        pass
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"Blast_ali_coverage_most_similar_ref_SP\" in file {}'.format(str(pathInputFile)))

#                if idxColIs_hit_HMM_CC >= 0:
#                    if (currSP.SPDetectedByHMM == 1):
#                        pass # already detected by Is_hit_HMM_JL
#                    else:
#                        parsedCell = row[idxColIs_hit_HMM_CC]
#                        parsedCell = int(parsedCell[0]) # [0] to get only the first character and avoid issues with number formatting according to different locale
#                        if (parsedCell == 0 or parsedCell == 1):
#                            currSP.SPDetectedByHMM = parsedCell
#                            if ( parsedCell == 1 and idxColBest_hmmprofile_CC >= 0 ):
#                                parsedCell = row[idxColBest_hmmprofile_CC]
#                                listFamilyTokens = re.split(regexSplitBlastFamilyIntoTokens, parsedCell)
#                                for currFamilyToken in listFamilyTokens:
#                                    if currFamilyToken:# add only if not empty or null
#                                        currSP.setSPFamilyFromHMM.add(currFamilyToken)
#                        else:
#                            raise RuntimeError('Input file error: could not parse value \"{}\" of column \"Is_hit_HMM_CC\" in row number {}'.format(parsedCell, countIterRow+1))
                listSPsParsed.list.append(currSP)
            countIterRow += 1
    csvfile.close()
    return listSPsParsed



def createDictIMEICEID2humanReadableIMEICEIIdentifier(listOfListAllICEsIMEsStructure):

    dictIMEICEID2humanReadableIMEICEIIdentifier = {}
    idxCountIdStructure = 1
    for currListAllICEsIMEsStructure in listOfListAllICEsIMEsStructure:
        for currICEIMEStructure in currListAllICEsIMEsStructure:
            #print ("internalIdentifier = {}".format(currICEIMEStructure.internalIdentifier))
            if currICEIMEStructure.internalIdentifier in dictIMEICEID2humanReadableIMEICEIIdentifier:
                raise RuntimeError('Error createDictIMEICEID2humanReadableIMEICEIIdentifier: currICEIMEStructure.internalIdentifier {} already in dictIMEICEID2humanReadableIMEICEIIdentifier'.format(currICEIMEStructure.internalIdentifier))
            if currICEIMEStructure.catStructConjModule == "unsure" and currICEIMEStructure.categoryRegardingIntegrase == "no integrase assigned":
                pass
            else :
                dictIMEICEID2humanReadableIMEICEIIdentifier[currICEIMEStructure.internalIdentifier] = "ID_"+str(idxCountIdStructure)
                idxCountIdStructure += 1
    return dictIMEICEID2humanReadableIMEICEIIdentifier

    #carefull, changed internalIdentifier to a non digit, can cause the following problems:
    # - src/EMStructure.py : problems maybe ??
        # DONE !!! if self.internalIdentifier >= otherICEsIMEsStructureToMerge.internalIdentifier:
        # DONE !!! if self.internalIdentifier == lowerBoundMergedInternalIdentifierIT and otherICEsIMEsStructureToMerge.internalIdentifier == higherBoundMergedInternalIdentifierIT:
        # DONE !!! and rest of the method
        # DONE !!! ICEsIMEsStructure.listLowerBoundHigherBoundStructureMergedInternalIdentifier.append(
        # DONE !!! if self.internalIdentifier < currICEIMEInConflict.internalIdentifier:



# This method prints the output file (tsv file) that details the list of ICEs and IMEs including their signature proteins.
def printAllICEsIMEsStructureToOutputFile(
        listOfListAllICEsIMEsStructure,
        # outputPrintCDSNumberInsteadOfProteinIdAndStart,
        outputFile,
        dictIMEICEID2humanReadableIMEICEIIdentifier):


    # print("\n** Detailed ICE / IME structures:", file=outputFile)
    listStHeaderToPrint = EMStructure.ICEsIMEsStructure.GetSummaryObjectHeaderAsTsv()
    stHeaderToPrintPrefix = listStHeaderToPrint[0]
    stHeaderToPrintPostfix = listStHeaderToPrint[1]
    print(stHeaderToPrintPrefix + stHeaderToPrintPostfix, file=outputFile)
    for idx, currListAllICEsIMEsStructure in enumerate(listOfListAllICEsIMEsStructure):
        # print("--------------------- Genomic region {}: {} - {} -------------------".format(idx+1, segmentIdx2startGenomicRegion[idx], segmentIdx2stopGenomicRegion[idx]), file=outputFile)
        for currICEIMEStructure in currListAllICEsIMEsStructure:
            if currICEIMEStructure.catStructConjModule == "unsure" and currICEIMEStructure.categoryRegardingIntegrase == "no integrase assigned":
                pass
            else :
                listStToPrint = currICEIMEStructure.GetSummaryObjectAsTsv(
                        # outputPrintCDSNumberInsteadOfProteinIdAndStart,
                        idx + 1,
                        dictIMEICEID2humanReadableIMEICEIIdentifier)
                stToPrintPrefix = listStToPrint[0]
                stToPrintPostfix = listStToPrint[1]
                print(stToPrintPrefix + stToPrintPostfix, file=outputFile)
                # print("{}".format(currICEIMEStructure.GetObjectAsJson(True, "")))

########
# MAIN #
########


def main():

    # Get folder of icescreen.py (folder src)
    # src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    # Change work directory to src
    # os.chdir(src_path)

    # Parse script arguments
    parser = argparse.ArgumentParser(description="ICEscreen (Searching for ICEs and IMEs)")
    # parser.add_argument("res", metavar="RES", help="CSV file (Results of clustering of hits of ICE Finder)", nargs=1, type=str)
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
            '-i',
            '--input',
            help="path to input file containsing signature protein (SP) data from the first part of the pipeline (.tsv)",
            required=True)
    required.add_argument(
            '-c',
            '--config',
            help="path to config file that contains the major settings for this script",
            required=True)
    required.add_argument(
            '-o',
            '--output',
            help="path to output file (tsv file)",
            required=True)
    required.add_argument(
            '-m',
            '--modified_input',
            help="path to output csv file: input cvs file + information on the ICE or IME structure",
            required=True)
    required.add_argument(
            '-l',
            '--log',
            help="path to the log file",
            required=True)
    required.add_argument(
            '-verbose',
            help="Turn on verbose mode ; default=False",
            action='store_true')

    args = parser.parse_args()
    pathInputFile = args.input
    pathOutputFile = args.output
    pathModifiedInputFile = args.modified_input
    pathLogFile = args.log
    verbose = args.verbose

    start_time = time.time()
    if verbose:
        print("\n\nStarting icescreen_OO.py at {}".format(datetime.datetime.now()))

    # Get configuration file
    config = configparser.ConfigParser()
    config.read(args.config)

    maxNumberCDSForSplitSPsByColocalizion = int(config["PARAMS"]["maxNumberCDSForSplitSPsByColocalizion"])
    maxNumberCDSForFilterIMESize = int(config["PARAMS"]["maxNumberCDSForFilterIMESize"])
    groupListSPintoICEsIMEsUsingFamilyInfo = config["PARAMS"]["groupListSPintoICEsIMEsUsingFamilyInfo"]
    useFamilyInfoToTryToResolveSPModuleConjConflict = config["PARAMS"]["useFamilyInfoToTryToResolveSPModuleConjConflict"]
    useDistanceCDSInfoToTryToResolveSPModuleConjConflict = config["PARAMS"]["useDistanceCDSInfoToTryToResolveSPModuleConjConflict"]
    moveSingleSPToCheck = config["PARAMS"]["moveSingleSPToCheck"]
    # outputPrintCDSNumberInsteadOfProteinIdAndStart = config["PARAMS"]["outputPrintCDSNumberInsteadOfProteinIdAndStart"]
    allowAdjacentIntegraseOnlyForSer = config["PARAMS"]["allowAdjacentIntegraseOnlyForSer"]
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance = int(config["PARAMS"]["useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance"])
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance = int(config["PARAMS"]["useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance"])
    maxOverlappingAliLenghtFragmentedSPs = int(config["PARAMS"]["maxOverlappingAliLenghtFragmentedSPs"])
    maxCDSDistanceToMergeFragmentedSPs = int(config["PARAMS"]["maxCDSDistanceToMergeFragmentedSPs"])
    EMStructure.BasicEMStructure.threshold_blast_ali_identity_perc_transfert_ICEFamily_to_structure_module_conj = int(config["PARAMS"]["threshold_blast_ali_identity_perc_transfert_ICEFamily_to_structure_module_conj"])

    # open and truncate outputFile and logFile
    outputFile = open(pathOutputFile, "w")
    modifiedInputFile = open(pathModifiedInputFile, "w")
    logFile = open(pathLogFile, "w")

    # print path to file
    # print("Input file: \"{}\"".format(pathInputFile))
    if verbose:
        print("Output file: \"{}\"".format(pathOutputFile))
        print("Modified input file: \"{}\"".format(pathModifiedInputFile))
        print("Log file: \"{}\"".format(pathLogFile))

    print("** Path to file:", file=logFile)
    print("Input file: \"{}\"".format(
            pathInputFile), file=logFile)
    print("Output file: \"{}\"".format(
            pathOutputFile), file=logFile)
    print("Modified input file: \"{}\"".format(
            pathModifiedInputFile), file=logFile)

    #print("\n** Script version: \"{}\"".format( scriptVersion), file=logFile)
    print("\n** Parameters:", file=logFile)
    print("maxNumberCDSForSplitSPsByColocalizion: \"{}\"".format(
            maxNumberCDSForSplitSPsByColocalizion), file=logFile)
    print("maxNumberCDSForFilterIMESize: \"{}\"".format(
            maxNumberCDSForFilterIMESize), file=logFile)
    print("groupListSPintoICEsIMEsUsingFamilyInfo: \"{}\"".format(
            groupListSPintoICEsIMEsUsingFamilyInfo), file=logFile)
    print("useFamilyInfoToTryToResolveSPModuleConjConflict: \"{}\"".format(
            useFamilyInfoToTryToResolveSPModuleConjConflict), file=logFile)
    print("useDistanceCDSInfoToTryToResolveSPModuleConjConflict: \"{}\"".format(
            useDistanceCDSInfoToTryToResolveSPModuleConjConflict), file=logFile)
    print("moveSingleSPToCheck: \"{}\"".format(
            moveSingleSPToCheck), file=logFile)
    # print("outputPrintCDSNumberInsteadOfProteinIdAndStart: \"{}\"".format(
    #         outputPrintCDSNumberInsteadOfProteinIdAndStart), file=logFile)
    print("allowAdjacentIntegraseOnlyForSer: \"{}\"".format(
            allowAdjacentIntegraseOnlyForSer), file=logFile)
    print("useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance: \"{}\"".format(
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance), file=logFile)
    print("useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance: \"{}\"".format(
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance), file=logFile)

    
    hasMultipleGenomeAccesion = commonMethods.determineIfResultSPFileHasMultipleGenomeAccesion(pathInputFile)


    SPsWholeGenome = parse_csv(
            # outputPrintCDSNumberInsteadOfProteinIdAndStart,
            pathInputFile,
            hasMultipleGenomeAccesion
            )
    
    SPsWholeGenome.sortListSPsByStart()

    # for testing
    # for currSP in SPsWholeGenome.list:
    #    print("{}".format(currSP.GetObjectAsJson()))

    listOfListOfColocalizedSPs = SPsWholeGenome.splitListOrderedSPs_colocalizeByMaxNumberCDS(maxNumberCDSForSplitSPsByColocalizion)
    listOfListAllICEsIMEsStructure = [] # [ICEsIMEsStructure] the outer broader list correspond to segments, while each inner sub-list correspond to ICEsIMEsStructure within the segments
    listOfListSPsLonelyIntegrases = [] # [SP] the outer broader list correspond to segments, while each inner sub-list correspond to integrase that are not assigned to any ICE or IME structures
    locusTagIntegrase2Comment = {} # locusTagIntegrase2Comment is a dictionary used to store the comments that will be visible in the output files for each integrase
    locusTagFinalize2Comment = {} # locusTagFinalize2Comment is a dictionary generated at the last step of the algorithm and  used to store the comments that will be visible in the output files for each SP
    locusTagMerge2Comment = {} # locusTagMerge2Comment is a dictionary generated at the merge step to find nested structures and used to store the comments that will be visible in the output files for each SP

    # make ListOfSubsegmentListColocalizedSPs accordingly to FragmentedSPs
    listOfSegmentListOfSubsegmentListColocalizedSPs_afterSecondSplitByFragmentedSPs = []
    for currListSPs in listOfListOfColocalizedSPs:
        currListSPs.fillUpIdxOfSPsInListSP()
        currListSPs.searchForFragmentedSPs(
            #currListSPs.list,
            locusTagMerge2Comment,
            locusTagIntegrase2Comment,
            groupListSPintoICEsIMEsUsingFamilyInfo,
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
            maxOverlappingAliLenghtFragmentedSPs,
            maxCDSDistanceToMergeFragmentedSPs
            )
        listOfListOfColocalizedSPs_newTakingIntoConsiderationSplitByFragmentedSPs = currListSPs.splitListOrderedSPs_colocalizeByFragmentedSPs(locusTagMerge2Comment)
        listOfSegmentListOfSubsegmentListColocalizedSPs_afterSecondSplitByFragmentedSPs.append(listOfListOfColocalizedSPs_newTakingIntoConsiderationSplitByFragmentedSPs)


    # The following algorithm generates listOfListAllICEsIMEsStructure and listOfListSPsLonelyIntegrases from listOfListOfColocalizedSPs
    #OLD for currListSPs in listOfListOfColocalizedSPs:
    for segmentListOfSubsegmentListColocalizedSPs_afterSecondSplitByFragmentedSPsIT in listOfSegmentListOfSubsegmentListColocalizedSPs_afterSecondSplitByFragmentedSPs: # search ICEIMEStrucutre in subsegments

        tmpBuildUp_listICEsIMEsStructures = []
        tmpBuildUp_currListSPs = hit.ListSPs()

        for currListSPs in segmentListOfSubsegmentListColocalizedSPs_afterSecondSplitByFragmentedSPsIT:

            currListSPs.fillUpIdxOfSPsInListSP()

            listSameFamilyMergeStructures = []  # listSameFamilyMergeStructures contains information about the SPs that belong to the same family and that should be grouped in priority.
            SPsInSameFamilyMergeStructures2SameFamilyMergeStructure = {}  # key = SP, value = EMStructureMerged
            if useFamilyInfoToTryToResolveSPModuleConjConflict == "YES":
                (listSameFamilyMergeStructures,
                SPsInSameFamilyMergeStructures2SameFamilyMergeStructure) = rulesMergeICEIMEStructures.buildSameFamilyMergeStructures(
                        currListSPs.list,
                        locusTagMerge2Comment,
                        groupListSPintoICEsIMEsUsingFamilyInfo)


            # seedICEsIMEsStructure is a method that initiate anchors from subsequent SPs of the conjugation module. Anchors will eventually be finalized as ICEs and IME structures.
            listICEsIMEsStructures = currListSPs.seedICEsIMEsStructure(
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    useFamilyInfoToTryToResolveSPModuleConjConflict,
                    useDistanceCDSInfoToTryToResolveSPModuleConjConflict,
                    allowAdjacentIntegraseOnlyForSer,
                    locusTagIntegrase2Comment,
                    SPsInSameFamilyMergeStructures2SameFamilyMergeStructure,
                    listSameFamilyMergeStructures,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)

            # integrase can now be SP.listSiblingFragmentedSP too
            # checkForObviousIntegraseUpstreamAndDownstreamToAdd
            rulesAddIntegrases.addObviousIntegraseUpstreamAndDownstream_priorMerging(
                listICEsIMEsStructures,
                currListSPs.list,
                locusTagIntegrase2Comment,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                allowAdjacentIntegraseOnlyForSer
            )

            # tryMergeSameFamilyStructures merge distant anchors of SPs of the conjugation module that are of the same family
            rulesMergeICEIMEStructures.tryMergeSameFamilyStructures(
                    listICEsIMEsStructures,
                    listSameFamilyMergeStructures,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    locusTagMerge2Comment,
                    locusTagIntegrase2Comment,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)

            # tryMergeSameFamilyStructures merge distant anchors of SPs of the conjugation module that are not of the same family but compatible
            rulesMergeICEIMEStructures.tryMergeNestedICEsIMEsStructures(
                    listICEsIMEsStructures,
                    maxNumberCDSForFilterIMESize,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    locusTagMerge2Comment,
                    locusTagIntegrase2Comment,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)

            for currICEsIMEsStructure in listICEsIMEsStructures:
                # tryResolveSPsInConflictAfterMergeEvents solves SPs that were attributed to multiple anchors initially and whose at least 1 anchor was involved in a merge with another anchor.
                currICEsIMEsStructure.tryResolveSPsInConflictAfterMergeEvents()

            # integrase can now be SP.listSiblingFragmentedSP too
            rulesAddIntegrases.addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure(
                    currListSPs.list,
                    listICEsIMEsStructures,
                    maxNumberCDSForFilterIMESize,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    set(),
                    0,
                    allowAdjacentIntegraseOnlyForSer,
                    locusTagIntegrase2Comment,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
            
            tmpBuildUp_listICEsIMEsStructures.extend(listICEsIMEsStructures)
            tmpBuildUp_currListSPs.list.extend(currListSPs.list)


        # reattach all subsegment into one listICEsIMEsStructures as well as all SP in new all containing currListSPs
        tmpBuildUp_listICEsIMEsStructures = sorted(tmpBuildUp_listICEsIMEsStructures, key=lambda x: x.listOrderedSPs[0].start, reverse=False)
        tmpBuildUp_currListSPs.sortListSPsByStart()
        tmpBuildUp_currListSPs.fillUpIdxOfSPsInListSP()

        dictIntegraseLocusTagFoundInStructure = {}  # dictIntegraseLocusTagFoundInStructure is used so to differenciate integrase that are assigned to structures and those who are not.
        for currICEsIMEsStructure in tmpBuildUp_listICEsIMEsStructures:
            if currICEsIMEsStructure.delMerging_idxListUpstreamStructure == -1:
                # finalizeICEIMEStruct assigns the type of the ICE or IME and check for potential errors in the structures
                currICEsIMEsStructure.finalizeICEIMEStruct(
                        tmpBuildUp_listICEsIMEsStructures,
                        maxNumberCDSForFilterIMESize,
                        moveSingleSPToCheck,
                        locusTagFinalize2Comment)
                # List integrases found in structure, so that we can find integrases not found in structure latter
                for currSPIntegrase in currICEsIMEsStructure.listIntegraseUpstream:
                    dictIntegraseLocusTagFoundInStructure[currSPIntegrase.locusTag] = 1
                for currSPIntegrase in currICEsIMEsStructure.listIntegraseDownstream:
                    dictIntegraseLocusTagFoundInStructure[currSPIntegrase.locusTag] = 1
                # print("{}".format(currICEsIMEsStructure.GetObjectAsJson(True, "")))

        # del if delMerging_idxListUpstreamStructure
        for i in range(len(tmpBuildUp_listICEsIMEsStructures) - 1, -1, -1):
            currICEsIMEsStructure = tmpBuildUp_listICEsIMEsStructures[i]
            if currICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
                del tmpBuildUp_listICEsIMEsStructures[i]
        # add to listOfListAllICEsIMEsStructure
        if len(tmpBuildUp_listICEsIMEsStructures) >= 1:
            segmentIdx2startGenomicRegion[len(listOfListAllICEsIMEsStructure)] = tmpBuildUp_currListSPs.list[0].start
            segmentIdx2stopGenomicRegion[len(listOfListAllICEsIMEsStructure)] = tmpBuildUp_currListSPs.list[-1].stop
            listOfListAllICEsIMEsStructure.append(tmpBuildUp_listICEsIMEsStructures)
            # find integrases not found in structure
            listIntegraseNotFoundInStructureSegment = []
            for currSP in tmpBuildUp_currListSPs.list:
                if (currSP.SPType in setIntegraseNames):  # currSP.SPType == "IntTyr" or currSP.SPType == "IntSer" or currSP.SPType == "DDE"
                    if currSP.locusTag in dictIntegraseLocusTagFoundInStructure:
                        pass
                    else:
                        listIntegraseNotFoundInStructureSegment.append(currSP)
                        # dictIntegraseNotFoundInStructure[currSP] = 1
            listOfListSPsLonelyIntegrases.append(listIntegraseNotFoundInStructureSegment)
        else:
            listOfListAllICEsIMEsStructure.append(tmpBuildUp_listICEsIMEsStructures)
            listOfListSPsLonelyIntegrases.append(tmpBuildUp_currListSPs.list)

    dictIMEICEID2humanReadableIMEICEIIdentifier = createDictIMEICEID2humanReadableIMEICEIIdentifier(listOfListAllICEsIMEsStructure)

    # print in output files
    printAllICEsIMEsStructureToOutputFile(
            listOfListAllICEsIMEsStructure,
            # outputPrintCDSNumberInsteadOfProteinIdAndStart,
            outputFile,
            dictIMEICEID2humanReadableIMEICEIIdentifier)
    
    
    # (totalNumberSP,
    #  totalNumberIntegrase,
    #  totalNumberUnaffectedIntegrase,
    #  totalNumberRelaxase,
    #  totalNumberUnaffectedRelaxase,
    #  totalNumberCoupling,
    #  totalNumberUnaffectedCoupling,
    #  totalNumberVirb4,
    #  totalNumberUnaffectedVirb4) = printAllICEsIMEsStructureToInputFile(
    (dictGenomeAccnum2totalNumberSP,
    dictGenomeAccnum2totalNumberIntegrase,
    dictGenomeAccnum2totalNumberUnaffectedIntegrase,
    dictGenomeAccnum2totalNumberRelaxase,
    dictGenomeAccnum2totalNumberUnaffectedRelaxase,
    dictGenomeAccnum2totalNumberCoupling,
    dictGenomeAccnum2totalNumberUnaffectedCoupling,
    dictGenomeAccnum2totalNumberVirb4,
    dictGenomeAccnum2totalNumberUnaffectedVirb4,
    dictGenomeAccnum2totalsetFragmentedSP
    ) = printAllICEsIMEsStructureToInputFile(
        listOfListAllICEsIMEsStructure,
        listOfListSPsLonelyIntegrases,
        pathInputFile,
        modifiedInputFile,
        locusTagIntegrase2Comment,
        locusTagFinalize2Comment,
        locusTagMerge2Comment,
        dictIMEICEID2humanReadableIMEICEIIdentifier,
        hasMultipleGenomeAccesion
        )
    pathSummaryFile = os.path.splitext(pathOutputFile)[0] + '.summary'
    summaryFile = open(pathSummaryFile, "w")
    # print summary file
    # printOverallStatsToSummaryFile(
    #         totalNumberSP,
    #         totalNumberIntegrase,
    #         totalNumberUnaffectedIntegrase,
    #         totalNumberRelaxase,
    #         totalNumberUnaffectedRelaxase,
    #         totalNumberCoupling,
    #         totalNumberUnaffectedCoupling,
    #         totalNumberVirb4,
    #         totalNumberUnaffectedVirb4,
    #         listOfListAllICEsIMEsStructure,
    #         maxNumberCDSForSplitSPsByColocalizion,
    #         maxNumberCDSForFilterIMESize,
    #         summaryFile)
    printOverallStatsToSummaryFile(
        dictGenomeAccnum2totalNumberSP,
        dictGenomeAccnum2totalNumberIntegrase,
        dictGenomeAccnum2totalNumberUnaffectedIntegrase,
        dictGenomeAccnum2totalNumberRelaxase,
        dictGenomeAccnum2totalNumberUnaffectedRelaxase,
        dictGenomeAccnum2totalNumberCoupling,
        dictGenomeAccnum2totalNumberUnaffectedCoupling,
        dictGenomeAccnum2totalNumberVirb4,
        dictGenomeAccnum2totalNumberUnaffectedVirb4,
        dictGenomeAccnum2totalsetFragmentedSP,
        listOfListAllICEsIMEsStructure,
        maxNumberCDSForSplitSPsByColocalizion,
        maxNumberCDSForFilterIMESize,
        summaryFile)

    summaryFile.close()

    outputFile.close()
    modifiedInputFile.close()
    logFile.close()

    # end of script
    elapsed_time_secs = time.time() - start_time
    if verbose:
        print("Execution took: %s secs (Wall clock time)\n\n" % timedelta(seconds=round(elapsed_time_secs)))


if __name__ == '__main__':
    main()
