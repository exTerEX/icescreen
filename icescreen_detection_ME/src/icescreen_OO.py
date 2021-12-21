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


########
# GLOBAL VARS #
########
#scriptVersion = "1.4"
# setIntegraseTyrNames = {"IntTyr", "TyrInt"}
setIntegraseTyrNames = {"Tyrosine integrase"}
# setIntegraseSerNames = {"IntSer", "SerInt"}
setIntegraseSerNames = {"Serine integrase"}
setIntegraseNames = set()
setIntegraseNames.update(setIntegraseTyrNames, setIntegraseSerNames)
# setIntegraseNames.add("DDE")
setIntegraseNames.add("DDE transposase")
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


# This method generates the content of the output file (csv file) that is the input file with columns added that reflect the output file (i.e. \"ICE IME Number\" and \"ICE IME Number (To review)\").
def printAllICEsIMEsStructureToInputFile(
        listOfListAllICEsIMEsStructure,
        listOfListSPsLonelyIntegrases,
        pathInputFile,
        modifiedInputFile,
        locusTagIntegrase2Comment,
        locusTagFinalize2Comment,
        locusTagMerge2Comment):

    printSegmentNumber = True

    totalNumberSP = 0
    totalNumberIntegrase = 0
    totalNumberUnaffectedIntegrase = 0
    totalNumberRelaxase = 0
    totalNumberUnaffectedRelaxase = 0
    totalNumberCoupling = 0
    totalNumberUnaffectedCoupling = 0
    totalNumberVirb4 = 0
    totalNumberUnaffectedVirb4 = 0
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

            for currSPLocusTags in currICEIMEStructure.setSPConjModuleLocusTagsToManuallyCheck:
                dictLocusTag2segmentNumber[currSPLocusTags] = idx + 1
                if currSPLocusTags in dictLocusTagNotSure2ICEIMEInternalId:
                    setInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTags]
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                else:
                    setInternalId = set()
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                    dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTags] = setInternalId

            for currSPLocusTags in currICEIMEStructure.setIntegraseLocusTagsToManuallyCheck:
                dictLocusTag2segmentNumber[currSPLocusTags] = idx + 1
                if currSPLocusTags in dictLocusTagNotSure2ICEIMEInternalId:
                    setInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTags]
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                else:
                    setInternalId = set()
                    setInternalId.add(currICEIMEStructure.internalIdentifier)
                    dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTags] = setInternalId

    for idx, currListSPsLonelyIntegrases in enumerate(listOfListSPsLonelyIntegrases):
        for currSPLonelyIntegrases in currListSPsLonelyIntegrases:
            dictLocusTag2segmentNumber[currSPLonelyIntegrases.locusTag] = idx + 1

    with open(pathInputFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        countIterRow = 0
        idxColCDS = -1
        idxColCDS_start = -1
        idxColCDS_Protein_Type = -1
        for row in reader:
            if countIterRow == 0:
                # check header
                countIterCol = 0
                for column in row:
                    if column == "CDS":
                        idxColCDS = countIterCol
                    elif column == "CDS_start":
                        idxColCDS_start = countIterCol
                    elif column == "CDS_Protein_type":
                        idxColCDS_Protein_Type = countIterCol
                    countIterCol += 1
                if printSegmentNumber:
                    print("ICE IME Number\tICE IME Number (To review)\tSegment Number\tComment regarding ICE IME structure\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)
                else:
                    print("ICE IME Number\tICE IME Number (To review)\tComment regarding ICE IME structure\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)

            else:
                currSPProteinId = ""
                currSPstart = ""
                currSPType = ""  # Coupling, Relaxase, Virb4, or integrase

                if idxColCDS >= 0:
                    parsedCell = row[idxColCDS]
                    currSPProteinId = parsedCell
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS\"')
                if idxColCDS_start >= 0:
                    parsedCell = row[idxColCDS_start]
                    parsedCell = int(parsedCell)
                    currSPstart = parsedCell
                else:
                    raise RuntimeError('printAllICEsIMEsStructureToInputFile error: missing mandatory column \"CDS_start\"')
                currSPLocusTag = currSPProteinId + "-" + str(currSPstart)
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

                strNotSure = ""
                if currSPLocusTag in dictLocusTagNotSure2ICEIMEInternalId:
                    currSetInternalId = dictLocusTagNotSure2ICEIMEInternalId[currSPLocusTag]
                    strNotSure = " ".join(str(i) for i in sorted(currSetInternalId))

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

                if printSegmentNumber:
                    strSegmentNumber = ""
                    if currSPLocusTag in dictLocusTag2segmentNumber:
                        strSegmentNumber = dictLocusTag2segmentNumber[currSPLocusTag]
                    print(strSure + "\t" + strNotSure + "\t" + str(strSegmentNumber) + "\t" + strLocusTag2CommentIT + "\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)
                else:
                    print(strSure + "\t" + strNotSure + "\t" + strLocusTag2CommentIT + "\t" + "\t".join(str(i) for i in row), file=modifiedInputFile)

                totalNumberSP += 1
                if currSPType == "Coupling protein":
                    totalNumberCoupling += 1
                    if not strSure:
                        totalNumberUnaffectedCoupling += 1
                elif currSPType == "Relaxase":
                    totalNumberRelaxase += 1
                    if not strSure:
                        totalNumberUnaffectedRelaxase += 1
                elif currSPType == "VirB4":
                    totalNumberVirb4 += 1
                    if not strSure:
                        totalNumberUnaffectedVirb4 += 1
                elif currSPType in setIntegraseNames:
                    totalNumberIntegrase += 1
                    if not strSure:
                        totalNumberUnaffectedIntegrase += 1
                else:
                    raise RuntimeError(
                            'printAllICEsIMEsStructureToInputFile error: unrecognized SPType {} for LocusTag {}'.format(
                                    currSPType, currSPLocusTag))

            countIterRow += 1

    csvfile.close()

    return (totalNumberSP,
            totalNumberIntegrase,
            totalNumberUnaffectedIntegrase,
            totalNumberRelaxase,
            totalNumberUnaffectedRelaxase,
            totalNumberCoupling,
            totalNumberUnaffectedCoupling,
            totalNumberVirb4,
            totalNumberUnaffectedVirb4)


# This method parse the input file into objects used internally in the program.
def parse_csv(outputPrintCDSNumberInsteadOfProteinIdAndStart, pathInputFile):
    listSPsParsed = hit.ListSPs()
    with open(pathInputFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        countIterRow = 0
        idxColCDS = -1
        idxColCDS_start = -1
        idxColCDS_end = -1
        idxColCDS_strand = -1
        idxColCDS_num = -1
        idxColCDS_Protein_Type = -1
        idxColhit_blast = -1
        idxColhit_HMM = -1
        # idxColhit_HMM_CC = -1
        # idxColElement_family = -1
        idxColICE_superfamily = -1
        idxColICE_family = -1
        idxColIME_family = -1

        idxColBest_hmmprofile = -1
        # idxColBest_hmmprofile_CC = -1
        for row in reader:
            if countIterRow == 0:
                # check header
                countIterCol = 0
                for column in row:
                    if column == "CDS":
                        idxColCDS = countIterCol
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
                    elif column == "hit_blast":
                        idxColhit_blast = countIterCol
                    elif column == "hit_HMM":
                        idxColhit_HMM = countIterCol
                    # elif column=="hit_HMM_CC":
                    #    idxColhit_HMM_CC = countIterCol
                    # elif column=="Element_family":
                    #    idxColElement_family = countIterCol
                    elif column == "ICE_superfamily":
                        idxColICE_superfamily = countIterCol
                    elif column == "ICE_family":
                        idxColICE_family = countIterCol
                    elif column == "IME_family":
                        idxColIME_family = countIterCol
                    elif column == "Profile_description":
                        # was elif column=="Best_hmmprofile":
                        idxColBest_hmmprofile = countIterCol
                    # elif column=="Best_hmmprofile_CC":
                    #    idxColBest_hmmprofile_CC = countIterCol
                    countIterCol += 1
                if ((idxColhit_blast + idxColhit_HMM) == -2):  # + idxColhit_HMM_CC
                    raise RuntimeError('Input file error: absence of at least one of the mandatory columns \"hit_blast\" or \"hit_HMM\"')
            else:
                currSP = hit.SP()
                if idxColCDS >= 0:
                    parsedCell = row[idxColCDS]
                    currSP.proteinId = parsedCell
                    if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                        currSP.locusTag = currSP.proteinId + "-" + str(currSP.start)
                    elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                        pass
                    else:
                        raise RuntimeError(
                                "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                                        outputPrintCDSNumberInsteadOfProteinIdAndStart))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS\"')
                if idxColCDS_start >= 0:
                    parsedCell = row[idxColCDS_start]
                    parsedCell = int(parsedCell)
                    currSP.start = parsedCell
                    if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                        currSP.locusTag = currSP.proteinId + "-" + str(currSP.start)
                    elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                        pass
                    else:
                        raise RuntimeError(
                                "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                                        outputPrintCDSNumberInsteadOfProteinIdAndStart))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_start\"')
                if idxColCDS_end >= 0:
                    parsedCell = row[idxColCDS_end]
                    parsedCell = int(parsedCell)
                    currSP.stop = parsedCell
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_end\"')
                if idxColCDS_strand >= 0:
                    parsedCell = row[idxColCDS_strand]
                    if (parsedCell == "+" or parsedCell == "-"):
                        currSP.strand = parsedCell
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"CDS_strand\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_strand\"')
                if idxColCDS_num >= 0:
                    parsedCell = row[idxColCDS_num]
                    parsedCell = int(parsedCell)
                    currSP.CDSPositionInGenome = parsedCell
                    if outputPrintCDSNumberInsteadOfProteinIdAndStart == "NO":
                        pass
                    elif outputPrintCDSNumberInsteadOfProteinIdAndStart == "YES":
                        currSP.locusTag = str(currSP.CDSPositionInGenome)
                    else:
                        raise RuntimeError(
                                "Error in parse_csv: unrecognized outputPrintCDSNumberInsteadOfProteinIdAndStart = {}".format(
                                        outputPrintCDSNumberInsteadOfProteinIdAndStart))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_num\"')
                if idxColCDS_Protein_Type >= 0:
                    parsedCell = row[idxColCDS_Protein_Type]
                    if (parsedCell == "Coupling protein" or parsedCell == "Relaxase" or parsedCell == "VirB4" or parsedCell in setIntegraseNames):  # parsedCell == "IntTyr" or parsedCell == "IntSer" or parsedCell == "DDE"
                        currSP.SPType = parsedCell
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"CDS_Protein_type\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                else:
                    raise RuntimeError('Input file error: missing mandatory column \"CDS_Protein_type\"')
                if idxColhit_blast >= 0:
                    parsedCell = row[idxColhit_blast]
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
                            if idxColICE_superfamily >= 0:
                                parsedCell = row[idxColICE_superfamily]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPICESuperFamilyFromBlast.add(parsedCell)
                            if idxColICE_family >= 0:
                                parsedCell = row[idxColICE_family]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPICEFamilyFromBlast.add(parsedCell)
                            if idxColIME_family >= 0:
                                parsedCell = row[idxColIME_family]
                                if parsedCell and parsedCell != "-":  # add only if not empty or null
                                    currSP.setSPIMEFamilyFromBlast.add(parsedCell)
                    else:
                        raise RuntimeError(
                                'Input file error: could not parse value \"{}\" of column \"hit_blast\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
                if idxColhit_HMM >= 0:
                    parsedCell = row[idxColhit_HMM]
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
                                'Input file error: could not parse value \"{}\" of column \"hit_HMM\" in row number {}'.format(
                                        parsedCell, countIterRow + 1))
#                if idxColhit_HMM_CC >= 0:
#                    if (currSP.SPDetectedByHMM == 1):
#                        pass # already detected by hit_HMM_JL
#                    else:
#                        parsedCell = row[idxColhit_HMM_CC]
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
#                            raise RuntimeError('Input file error: could not parse value \"{}\" of column \"hit_HMM_CC\" in row number {}'.format(parsedCell, countIterRow+1))
                listSPsParsed.list.append(currSP)
            countIterRow += 1
    csvfile.close()
    return listSPsParsed


# This method prints the output file (tsv file) that details the list of ICEs and IMEs including their signature proteins.
def printAllICEsIMEsStructureToOutputFile(
        listOfListAllICEsIMEsStructure,
        outputPrintCDSNumberInsteadOfProteinIdAndStart,
        outputFile):

    # print("\n** Detailed ICE / IME structures:", file=outputFile)
    listStHeaderToPrint = EMStructure.ICEsIMEsStructure.GetSummaryObjectHeaderAsTsv()
    stHeaderToPrintPrefix = listStHeaderToPrint[0]
    stHeaderToPrintPostfix = listStHeaderToPrint[1]
    print(stHeaderToPrintPrefix + stHeaderToPrintPostfix, file=outputFile)
    for idx, currListAllICEsIMEsStructure in enumerate(listOfListAllICEsIMEsStructure):
        # print("--------------------- Genomic region {}: {} - {} -------------------".format(idx+1, segmentIdx2startGenomicRegion[idx], segmentIdx2stopGenomicRegion[idx]), file=outputFile)
        for currICEIMEStructure in currListAllICEsIMEsStructure:
            listStToPrint = currICEIMEStructure.GetSummaryObjectAsTsv(
                    outputPrintCDSNumberInsteadOfProteinIdAndStart,
                    idx + 1)
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
            help="path to output file (tsv file) that is the input file with columns added that reflect the output file (i.e. \"ICE IME Number\" and \"ICE IME Number (To review)\")",
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
    outputPrintCDSNumberInsteadOfProteinIdAndStart = config["PARAMS"]["outputPrintCDSNumberInsteadOfProteinIdAndStart"]
    allowAdjacentIntegraseOnlyForSer = config["PARAMS"]["allowAdjacentIntegraseOnlyForSer"]
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance = int(
            config["PARAMS"]["useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance"])
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance = int(
            config["PARAMS"]["useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance"])

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
    print("outputPrintCDSNumberInsteadOfProteinIdAndStart: \"{}\"".format(
            outputPrintCDSNumberInsteadOfProteinIdAndStart), file=logFile)
    print("allowAdjacentIntegraseOnlyForSer: \"{}\"".format(
            allowAdjacentIntegraseOnlyForSer), file=logFile)
    print("useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance: \"{}\"".format(
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance), file=logFile)
    print("useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance: \"{}\"".format(
            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance), file=logFile)

    SPsWholeGenome = parse_csv(
            outputPrintCDSNumberInsteadOfProteinIdAndStart,
            pathInputFile)
    SPsWholeGenome.sortListSPsByStart()

    # for testing
    # for currSP in SPsWholeGenome.list:
    #    print("{}".format(currSP.GetObjectAsJson()))

    listOfListOfColocalizedSPs = SPsWholeGenome.splitListOrderedSPsByColocalizion(maxNumberCDSForSplitSPsByColocalizion)
    listOfListAllICEsIMEsStructure = []  # [ICEsIMEsStructure] the outer broader list correspond to segments, while each inner sub-list correspond to ICEsIMEsStructure within the segments
    listOfListSPsLonelyIntegrases = []  # [SP] the outer broader list correspond to segments, while each inner sub-list correspond to integrase that are not assigned to any ICE or IME structures
    locusTagIntegrase2Comment = {}  # locusTagIntegrase2Comment is a dictionary used to store the comments that will be visible in the output files for each integrase
    locusTagFinalize2Comment = {}  # locusTagFinalize2Comment is a dictionary generated at the last step of the algorithm and  used to store the comments that will be visible in the output files for each SP
    locusTagMerge2Comment = {}  # locusTagMerge2Comment is a dictionary generated at the merge step to find nested structures and used to store the comments that will be visible in the output files for each SP

    # The following algorithm generates listOfListAllICEsIMEsStructure and listOfListSPsLonelyIntegrases from listOfListOfColocalizedSPs
    for currListSPs in listOfListOfColocalizedSPs:
        # for testing
        # print("--------------------- Colocalized: -------------------")
        # for currSP in currListSPs.list:
        #    print("{}".format(currSP.GetObjectAsJson()))
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

        # addSPIntegraseUpstreamAndDownstream add integrase to anchor to create ICE and IME structures
        rulesAddIntegrases.addSPIntegraseUpstreamAndDownstream(
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

        dictIntegraseLocusTagFoundInStructure = {}  # dictIntegraseLocusTagFoundInStructure is used so to differenciate integrase that are assigned to structures and those who are not.
        for currICEsIMEsStructure in listICEsIMEsStructures:
            if currICEsIMEsStructure.delMerging_idxListUpstreamStructure == -1:
                # finalizeICEIMEStruct assigns the type of the ICE or IME and check for potential errors in the structures
                currICEsIMEsStructure.finalizeICEIMEStruct(
                        listICEsIMEsStructures,
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
        for i in range(len(listICEsIMEsStructures) - 1, -1, -1):
            currICEsIMEsStructure = listICEsIMEsStructures[i]
            if currICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
                del listICEsIMEsStructures[i]
        # add to listOfListAllICEsIMEsStructure
        if len(listICEsIMEsStructures) >= 1:
            segmentIdx2startGenomicRegion[len(listOfListAllICEsIMEsStructure)] = currListSPs.list[0].start
            segmentIdx2stopGenomicRegion[len(listOfListAllICEsIMEsStructure)] = currListSPs.list[-1].stop
            listOfListAllICEsIMEsStructure.append(listICEsIMEsStructures)
            # find integrases not found in structure
            listIntegraseNotFoundInStructureSegment = []
            for currSP in currListSPs.list:
                if (currSP.SPType in setIntegraseNames):  # currSP.SPType == "IntTyr" or currSP.SPType == "IntSer" or currSP.SPType == "DDE"
                    if currSP.locusTag in dictIntegraseLocusTagFoundInStructure:
                        pass
                    else:
                        listIntegraseNotFoundInStructureSegment.append(currSP)
                        # dictIntegraseNotFoundInStructure[currSP] = 1
            listOfListSPsLonelyIntegrases.append(listIntegraseNotFoundInStructureSegment)

        else:
            listOfListAllICEsIMEsStructure.append(listICEsIMEsStructures)
            listOfListSPsLonelyIntegrases.append(currListSPs.list)

    # print in output files
    printAllICEsIMEsStructureToOutputFile(
            listOfListAllICEsIMEsStructure,
            outputPrintCDSNumberInsteadOfProteinIdAndStart,
            outputFile)
    (totalNumberSP,
     totalNumberIntegrase,
     totalNumberUnaffectedIntegrase,
     totalNumberRelaxase,
     totalNumberUnaffectedRelaxase,
     totalNumberCoupling,
     totalNumberUnaffectedCoupling,
     totalNumberVirb4,
     totalNumberUnaffectedVirb4) = printAllICEsIMEsStructureToInputFile(
             listOfListAllICEsIMEsStructure,
             listOfListSPsLonelyIntegrases,
             pathInputFile,
             modifiedInputFile,
             locusTagIntegrase2Comment,
             locusTagFinalize2Comment,
             locusTagMerge2Comment)
    pathSummaryFile = os.path.splitext(pathOutputFile)[0] + '.summary'
    summaryFile = open(pathSummaryFile, "w")
    # print summary file
    printOverallStatsToSummaryFile(
            totalNumberSP,
            totalNumberIntegrase,
            totalNumberUnaffectedIntegrase,
            totalNumberRelaxase,
            totalNumberUnaffectedRelaxase,
            totalNumberCoupling,
            totalNumberUnaffectedCoupling,
            totalNumberVirb4,
            totalNumberUnaffectedVirb4,
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
