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
# Remark: copy similar file in icescreen_formatting/ and icescreen_detection_ME/src/ and icescreen_detection_SP/src/

import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# https://stackoverflow.com/questions/6821156/how-to-find-range-overlap-in-python
def overlap_bounded_intervals(start1, end1, start2, end2, genomeAccessionRank1, genomeAccessionRank2):
    if genomeAccessionRank1 != genomeAccessionRank2:
        return 0
    else:
        """how much does the range (start1, end1) overlap with (start2, end2)"""
        return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)
# def range_intersect(r1, r2):
#     return range(max(r1.start,r2.start), min(r1.stop,r2.stop)) or None


def makeCompositeUniqLocusTag(isMultiGenbankFile, locusTagSent, proteinIdSent, genomeAccessionSent, startSent):
    strLocusTagToReturn = ""

    if len(locusTagSent) == 0 or locusTagSent is None or locusTagSent == "-":
        if isMultiGenbankFile: 
            strLocusTagToReturn = genomeAccessionSent + "-" + proteinIdSent + "-" + str(startSent)
        else:
            strLocusTagToReturn = proteinIdSent + "-" + str(startSent)
    else :
        strLocusTagToReturn = locusTagSent
        #return locusTagSent
    
    if len(strLocusTagToReturn) == 0:
        raise RuntimeError('Error makeCompositeUniqLocusTag: len(strLocusTagToReturn) == 0 for proteinId {} genomeAccession {} start {} '.format(str(proteinIdSent), genomeAccessionSent, str(startSent)))
    
    #print("makeCompositeUniqLocusTag {} from isMultiGenbankFile {}, locusTagSent {}, proteinIdSent {}, genomeAccessionSent {}, startSent {}".format(strLocusTagToReturn, str(isMultiGenbankFile), locusTagSent, proteinIdSent, genomeAccessionSent, str(startSent)))

    return strLocusTagToReturn


def parseListGenomeAccessionFromGenbankFile(pathGbInputFile):
    hasMultipleGenomeAccesion = True
    listGenomeAccessionFromFastaFile = []
    dictGenomeIds = {}
    if os.path.isfile(pathGbInputFile):
        pathGbInputFile_handle = open (pathGbInputFile, 'r', encoding='utf-8')        
        record_iterator = SeqIO.parse(pathGbInputFile_handle, "genbank")
        for record in record_iterator:
            genome_accession = record.id
            if len(genome_accession) == 0:
                raise RuntimeError("Error in parseListGenomeAccessionFromGenbankFile: empty genome_accession in file {}".format(
                    str(pathGbInputFile)))
            if genome_accession in dictGenomeIds :
                raise RuntimeError("Error in parseListGenomeAccessionFromGenbankFile: duplicate genome_accession {} in file {}".format(
                    str(genome_accession),
                    str(pathGbInputFile)))
            else :
                dictGenomeIds[genome_accession] = 1
        pathGbInputFile_handle.close()
    listGenomeAccessionFromFastaFile = list(dictGenomeIds.keys())

    # For fasta file
    # with open(pathGbInputFile, 'r', encoding='utf8') as fastafileIT:
    #     for line in fastafileIT:
    #         line = line.strip()
    #         if line.startswith(">"):
    #             # >1&locus_tag=YYYYY1&protein_id=XXXXX1&genome_accession=NZ_TEST1&genome_accession_rank=1|+|1..2
    #             matchLineGenomeAccessionIT = re.match(r"^.+&genome_accession=(.+)&genome_accession_rank=.+$", line)
    #             if matchLineGenomeAccessionIT:
    #                 genomeAccessionIT = matchLineGenomeAccessionIT.group(1)
    #                 if genomeAccessionIT in setGenomeAccessionProcessed :
    #                     pass
    #                 else :
    #                     listGenomeAccessionFromFastaFile.append(genomeAccessionIT)
    #                     setGenomeAccessionProcessed.add(genomeAccessionIT)
    #             else :
    #                 raise RuntimeError('error parse_genbank: line does not match the regex: {} '.format(str(line)))
                
    if len(listGenomeAccessionFromFastaFile) == 0:
        raise RuntimeError('error parseListGenomeAccessionFromGenbankFile: len(listGenomeAccessionFromFastaFile) == 0 for file: {} '.format(str(pathGbInputFile)))
    elif len(listGenomeAccessionFromFastaFile) == 1:
        hasMultipleGenomeAccesion = False
    else:
        hasMultipleGenomeAccesion = True
    return (hasMultipleGenomeAccesion, listGenomeAccessionFromFastaFile)


#return boolean hasMultipleGenomeAccesion
def determineIfResultSPFileHasMultipleGenomeAccesion(pathInputFile):
    dictGenomeAccesion = {}
    with open(pathInputFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        countIterRow = 0
        idxColGenome_accession = -1
        for row in reader:
            if countIterRow == 0:
                # check header
                countIterCol = 0
                for column in row:
                    if column == "Genome_accession":
                        idxColGenome_accession = countIterCol
                    countIterCol += 1
            else:
                if idxColGenome_accession >= 0:
                    parsedCell = row[idxColGenome_accession]
                    if len(parsedCell) > 0:
                        if parsedCell not in dictGenomeAccesion:
                            dictGenomeAccesion[parsedCell] = 1
                    else:
                        raise RuntimeError('determineIfResultSPFileHasMultipleGenomeAccesion error: missing mandatory information on Genome_accession for row {} in file {}. idxColGenome_accession = {}'.format(str(row), str(pathInputFile), str(idxColGenome_accession)))
                else:
                    raise RuntimeError('determineIfResultSPFileHasMultipleGenomeAccesion error: missing mandatory column \"Genome_accession\" in file {}'.format(str(pathInputFile)))
            countIterRow += 1

    if len(dictGenomeAccesion) == 0:
        raise RuntimeError('determineIfResultSPFileHasMultipleGenomeAccesion error: len(dictGenomeAccesion) == 0 for file {}'.format(str(pathInputFile)))
    elif len(dictGenomeAccesion) == 1:
        return False
    else:
        return True

def search_feature(item, feature):
    """Search qualifier of a feature (:class:`Bio.SeqFeature.SeqFeature`) and
    return it formatted:
        - Remove comma: ',' -> ''
        - Replace space by underscore: ' ' -> '_'
        - Replace vertical bar by underscore: '|' -> '_'

    :param item: Name of the qualifier searched
    :param feature: The feature
    :type item: :class:`str`
    :type feature: Object of :class:`Bio.SeqFeature.SeqFeature`
    :return: The value of the qualifier
    :rtype: :class:`str`
    """
    res = feature.qualifiers[item][0]
    res = formatQualifier(res)
    return(res)

def formatQualifier(strQualifierToFormat):
    strQualifierToFormat = strQualifierToFormat.replace(',', '')
    strQualifierToFormat = strQualifierToFormat.replace(' ', '_')
    strQualifierToFormat = strQualifierToFormat.replace('|', '_')
    return(strQualifierToFormat)