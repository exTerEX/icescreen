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
import rulesSeedSPExtension
import hit
import icescreen_OO
import commonMethods
import itertools

# SP of similar family are preferrably grouped in an anchor of conj module
def buildSameFamilyMergeStructures(currListSPs, locusTagMerge2Comment,
                                   groupListSPintoICEsIMEsUsingFamilyInfo):
    listSameFamilyMergeStructures = []  # EMStructures
    SPsInSameFamilyMergeStructures2SameFamilyMergeStructure = {}
    # key = SP, value = EMStructureMerged

    family2ListSPs = {}
    for currSP in currListSPs:
        if (currSP.SPType in icescreen_OO.setIntegraseNames):
            continue
        if currSP.SPDetectedByBlast == 1:

            for familyFromBlastIT in currSP.setSPICEFamilyFromBlast:
                # print("setSPICEFamilyFromBlast: {} for {}".format(
                # familyFromBlastIT, currSP.locusTag))
                if familyFromBlastIT in family2ListSPs:
                    currListSPsFromFamily2ListSPs = family2ListSPs[familyFromBlastIT]
                    currListSPsFromFamily2ListSPs.append(currSP)
                else:
                    currListSPsFromFamily2ListSPs = []
                    currListSPsFromFamily2ListSPs.append(currSP)
                    family2ListSPs[familyFromBlastIT] = currListSPsFromFamily2ListSPs

            for familyFromBlastIT in currSP.setSPIMESuperFamilyFromBlast:
                # print("setSPIMESuperFamilyFromBlast: {} for {}".format(
                # familyFromBlastIT, currSP.locusTag))

                if familyFromBlastIT in family2ListSPs:
                    currListSPsFromFamily2ListSPs = family2ListSPs[familyFromBlastIT]
                    currListSPsFromFamily2ListSPs.append(currSP)
                else:
                    currListSPsFromFamily2ListSPs = []
                    currListSPsFromFamily2ListSPs.append(currSP)
                    family2ListSPs[familyFromBlastIT] = currListSPsFromFamily2ListSPs

    # calculate number CDSs between registred SP of the same family
    distBetweenSPs2duoSPToMerge = {}
    for keyFamily, valueListSPs in family2ListSPs.items():
        # print(" ** distBetweenSPs2duoSPToMerge: keyFamily = {}"
        # .format(keyFamily))
        if len(valueListSPs) <= 1:
            continue
        for currIdxSP, currValueSP in enumerate(valueListSPs):
            # print("\tcurrValueSP = {}".format(currValueSP.locusTag))
            if currIdxSP != 0:  # not the first element
                # was CDSPositionInGenome  instead of idxInListSP ?
                distBetweenSPs = abs(currValueSP.idxInListSP -
                                     valueListSPs[currIdxSP - 1].idxInListSP)
                tupleSPsIT = (valueListSPs[currIdxSP - 1], currValueSP)
                # print("\tdistBetweenSPs = {}".format(str(distBetweenSPs)))
                if distBetweenSPs in distBetweenSPs2duoSPToMerge:
                    currListSPsForDistBetweenSPs2duoSPToMerge = distBetweenSPs2duoSPToMerge[distBetweenSPs]
                    currListSPsForDistBetweenSPs2duoSPToMerge.append(
                            tupleSPsIT)
                else:
                    currListSPsForDistBetweenSPs2duoSPToMerge = []
                    currListSPsForDistBetweenSPs2duoSPToMerge.append(
                            tupleSPsIT)
                    distBetweenSPs2duoSPToMerge[distBetweenSPs] = currListSPsForDistBetweenSPs2duoSPToMerge

    # start by merging the closer ones
    for sortedKeyDistBetweenSPs2duoSPToMerge in sorted(
            distBetweenSPs2duoSPToMerge):
        listDuoSPToMerge = distBetweenSPs2duoSPToMerge[sortedKeyDistBetweenSPs2duoSPToMerge]
        for (currOneInDuoSPToMerge, currTwoInDuoSPToMerge) in listDuoSPToMerge:
            doAddAMergeComment = True
            distBetweenSPs = abs(
                    currOneInDuoSPToMerge.idxInListSP -
                    currTwoInDuoSPToMerge.idxInListSP)
            if distBetweenSPs == 1:
                doAddAMergeComment = False

            if currOneInDuoSPToMerge in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and currTwoInDuoSPToMerge in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # both SPs have already be merged
                if doAddAMergeComment:
                    commentITToAdd = "The SPs {} and {} of same family already have both been merged with other SPs of the same family so are not merged together. ".format(
                            currOneInDuoSPToMerge.locusTag, currTwoInDuoSPToMerge.locusTag)
                    # if commentITToAdd not in keyICEIMEPrimary.comment:
                    #    keyICEIMEPrimary.comment += commentITToAdd
                    icescreen_OO.addCommentToLocusTag2Comment(currOneInDuoSPToMerge.locusTag, commentITToAdd, locusTagMerge2Comment)
                    icescreen_OO.addCommentToLocusTag2Comment(currTwoInDuoSPToMerge.locusTag, commentITToAdd, locusTagMerge2Comment)

            elif currOneInDuoSPToMerge in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and currTwoInDuoSPToMerge not in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # currOneInDuoSPToMerge already merged, not currTwoInDuoSPToMerge
                # try to add currTwoInDuoSPToMerge in structure of currOneInDuoSPToMerge
                currSameFamilyMergeStructure = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currOneInDuoSPToMerge]
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currSameFamilyMergeStructure,
                        currTwoInDuoSPToMerge,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        False,
                        False,
                        False)
                if greenLightAddSPConjugaisonModule:
                    if doAddAMergeComment:
                        commentITToAdd = "Conjugaison module SP {} is being associated with same family {}. ".format(
                                currTwoInDuoSPToMerge.locusTag,
                                hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))

                        if commentITToAdd not in currSameFamilyMergeStructure.comment:
                            currSameFamilyMergeStructure.comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currTwoInDuoSPToMerge.locusTag, commentITToAdd, locusTagMerge2Comment)
                    # print("HERE1: "+currTwoInDuoSPToMerge.locusTag+"="+str(currSameFamilyMergeStructure.internalIdentifier))
                    currSameFamilyMergeStructure.addSPToConjugaisonModule(currTwoInDuoSPToMerge)
                    SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currTwoInDuoSPToMerge] = currSameFamilyMergeStructure

            elif currOneInDuoSPToMerge not in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and currTwoInDuoSPToMerge in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # currTwoInDuoSPToMerge already merged, not currOneInDuoSPToMerge
                # try to add currOneInDuoSPToMerge in structure of currTwoInDuoSPToMerge
                currSameFamilyMergeStructure = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currTwoInDuoSPToMerge]
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currSameFamilyMergeStructure,
                        currOneInDuoSPToMerge,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        False,
                        False,
                        False)
                if greenLightAddSPConjugaisonModule:
                    if doAddAMergeComment:
                        commentITToAdd = "Conjugaison module SP {} is being associated with same family {}. ".format(
                                currOneInDuoSPToMerge.locusTag,
                                hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))

                        if commentITToAdd not in currSameFamilyMergeStructure.comment:
                            currSameFamilyMergeStructure.comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currOneInDuoSPToMerge.locusTag, commentITToAdd, locusTagMerge2Comment)
                    # print("HERE2.2: "+currTwoInDuoSPToMerge.locusTag+"="+str(currSameFamilyMergeStructure.internalIdentifier))
                    currSameFamilyMergeStructure.addSPToConjugaisonModule(currOneInDuoSPToMerge)
                    SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currOneInDuoSPToMerge] = currSameFamilyMergeStructure

            else:
                # currTwoInDuoSPToMerge and currOneInDuoSPToMerge not already merged
                currSameFamilyMergeStructure = EMStructure.ICEsIMEsStructure(True)
                currSameFamilyMergeStructure.addSPToConjugaisonModule(currOneInDuoSPToMerge)
                # print("HERE3.1: "+currOneInDuoSPToMerge.locusTag+"="+str(currSameFamilyMergeStructure.internalIdentifier))
                SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currOneInDuoSPToMerge] = currSameFamilyMergeStructure
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currSameFamilyMergeStructure,
                        currTwoInDuoSPToMerge,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        False,
                        False,
                        False)
                if greenLightAddSPConjugaisonModule:
                    if doAddAMergeComment:
                        commentITToAdd = "Conjugaison module SP {} is being associated with same family {}. ".format(
                                currTwoInDuoSPToMerge.locusTag,
                                hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))

                        if commentITToAdd not in currSameFamilyMergeStructure.comment:
                            currSameFamilyMergeStructure.comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(currTwoInDuoSPToMerge.locusTag, commentITToAdd, locusTagMerge2Comment)
                    # print("HERE3.2: "+currTwoInDuoSPToMerge.locusTag+"="+str(currSameFamilyMergeStructure.internalIdentifier))
                    currSameFamilyMergeStructure.addSPToConjugaisonModule(currTwoInDuoSPToMerge)
                    SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currTwoInDuoSPToMerge] = currSameFamilyMergeStructure
                listSameFamilyMergeStructures.append(currSameFamilyMergeStructure)

    # loop through all SameFamilyMergeStructure and remove those len(.listOrderedSPs) <= 1
    for i in range(len(listSameFamilyMergeStructures) - 1, -1, -1):
        if len(listSameFamilyMergeStructures[i].listOrderedSPs) <= 1:
            for currSPRemovedFromMergeStruct in listSameFamilyMergeStructures[i].listOrderedSPs:
                del SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currSPRemovedFromMergeStruct]
            del listSameFamilyMergeStructures[i]

    # for SPIT, EMStructureMergedIT  in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure.items():  # key = SP, value = EMStructureMerged
    #    print(SPIT.locusTag+" | "+str(EMStructureMergedIT.internalIdentifier))
    # print("----------------------")

    # look around SPs registred for merge for SP greenLightAddSPConjugaisonModule without family
    for currSameFamilyMergeStructure in listSameFamilyMergeStructures:
        # print("Look around greenLightAddSPConjugaisonModule "+str(currSameFamilyMergeStructure.internalIdentifier))
        setGreenLightedSPConjugaisonModuleWithoutFamily = set()
        for currSPameFamilyMergeStructure in currSameFamilyMergeStructure.listOrderedSPs:
            # deal with upstream SP
            idxUpstreamSPIT = currSPameFamilyMergeStructure.idxInListSP - 1
            if idxUpstreamSPIT >= 0:
                currUpstreamSPIT = currListSPs[idxUpstreamSPIT]
                # preconditionsCompatibleNeigbhorFamilyStatus = False
                # if len(currUpstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currUpstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currUpstreamSPIT.setSPIMESuperFamilyFromBlast) == 0:
                #    preconditionsCompatibleNeigbhorFamilyStatus = True
                preconditionsCompatibleNeigbhorFamilyStatus = True
                if len(currSameFamilyMergeStructure.listVirB4) > 0 and len(currUpstreamSPIT.setSPIMESuperFamilyFromBlast) > 0:
                    preconditionsCompatibleNeigbhorFamilyStatus = False
                if currUpstreamSPIT in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and currSPameFamilyMergeStructure in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currUpstreamSPIT] == SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currSPameFamilyMergeStructure]:
                        preconditionsCompatibleNeigbhorFamilyStatus = False
                if preconditionsCompatibleNeigbhorFamilyStatus is True and currUpstreamSPIT in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    preconditionsCompatibleNeigbhorFamilyStatus = False
                    # commentITToAdd = "Conjugaison module SP {} could be associated with downstream neighbor SPs {} but was not because it is already associated with another same family structure. "\
                    # .format(currUpstreamSPIT.locusTag, hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))
                    # if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    #    currSameFamilyMergeStructure.comment += commentITToAdd
                    # icescreen_OO.addCommentToLocusTag2Comment(currUpstreamSPIT.locusTag, commentITToAdd, locusTagMerge2Comment)

                if preconditionsCompatibleNeigbhorFamilyStatus is True:
                    if currUpstreamSPIT not in currSameFamilyMergeStructure.listOrderedSPs:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                                currSameFamilyMergeStructure,
                                currUpstreamSPIT,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                False,
                                False,
                                False)
                        if greenLightAddSPConjugaisonModule:
                            setGreenLightedSPConjugaisonModuleWithoutFamily.add(currUpstreamSPIT)
                            # print("added currUpstreamSPIT "+currUpstreamSPIT.locusTag+"to setGreenLightedSPConjugaisonModuleWithoutFamily")

            # deal with downstream SP
            idxDownstreamSPIT = currSPameFamilyMergeStructure.idxInListSP + 1
            if idxDownstreamSPIT < len(currListSPs):
                currDownstreamSPIT = currListSPs[idxDownstreamSPIT]

                # preconditionsCompatibleNeigbhorFamilyStatus = False
                # if len(currDownstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPIMESuperFamilyFromBlast) == 0:
                #    preconditionsCompatibleNeigbhorFamilyStatus = True
                preconditionsCompatibleNeigbhorFamilyStatus = True
                if len(currSameFamilyMergeStructure.listVirB4) > 0 and len(currDownstreamSPIT.setSPIMESuperFamilyFromBlast) > 0:
                    preconditionsCompatibleNeigbhorFamilyStatus = False
                if currDownstreamSPIT in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and currSPameFamilyMergeStructure in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    if SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currDownstreamSPIT] == SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currSPameFamilyMergeStructure]:
                        preconditionsCompatibleNeigbhorFamilyStatus = False
                if preconditionsCompatibleNeigbhorFamilyStatus is True and currDownstreamSPIT in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                    preconditionsCompatibleNeigbhorFamilyStatus = False
                    # commentITToAdd = "Conjugaison module SP {} could be associated with upstream neighbor SPs {} but was not because it is already associated with another same family structure. "\
                    # .format(currDownstreamSPIT.locusTag, hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))
                    # if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    #    currSameFamilyMergeStructure.comment += commentITToAdd
                    # icescreen_OO.addCommentToLocusTag2Comment(currDownstreamSPIT.locusTag, commentITToAdd, locusTagMerge2Comment)

                if preconditionsCompatibleNeigbhorFamilyStatus is True:
                    # if len(currDownstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPIMESuperFamilyFromBlast) == 0:
                    if currDownstreamSPIT not in currSameFamilyMergeStructure.listOrderedSPs:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                                currSameFamilyMergeStructure,
                                currDownstreamSPIT,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                False,
                                False,
                                False)
                        if greenLightAddSPConjugaisonModule:
                            setGreenLightedSPConjugaisonModuleWithoutFamily.add(currDownstreamSPIT)
                            # print("added currDownstreamSPIT "+currDownstreamSPIT.locusTag+"to setGreenLightedSPConjugaisonModuleWithoutFamily")

        if len(setGreenLightedSPConjugaisonModuleWithoutFamily) == 0:
            pass
        elif len(setGreenLightedSPConjugaisonModuleWithoutFamily) == 1:
            # just one choice possible, add it

            currSPCompatibleNeigbhor = next(iter(setGreenLightedSPConjugaisonModuleWithoutFamily))
            # for currSPCompatibleNeigbhor in setGreenLightedSPConjugaisonModuleWithoutFamily:

            # for SPIT, EMStructureMergedIT  in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure.items():  # key = SP, value = EMStructureMerged
            #    print(SPIT.locusTag+" | "+str(EMStructureMergedIT.internalIdentifier))

            # if add a SP that way, check if present in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and give priority to proximity ?
            if currSPCompatibleNeigbhor in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                currSameFamilyMergeStructureToAlter = SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currSPCompatibleNeigbhor]
                currSameFamilyMergeStructureToAlter.removeSPConjugaisonModule(currSPCompatibleNeigbhor, None, False, True)
                commentITToAdd = "Conjugaison module SP {} could be associated with same family SPs {} but was not because it is a compatible neigbhor with another same family structure. ".format(
                        currSPCompatibleNeigbhor.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructureToAlter.listOrderedSPs))
                # print(commentITToAdd)
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                if commentITToAdd not in currSameFamilyMergeStructureToAlter.comment:
                    currSameFamilyMergeStructureToAlter.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(currSPCompatibleNeigbhor.locusTag, commentITToAdd, locusTagMerge2Comment)
                commentITToDel = "Conjugaison module SP {} is being associated with same family {}. ".format(
                        currSPCompatibleNeigbhor.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructureToAlter.listOrderedSPs))
                currSameFamilyMergeStructureToAlter.comment = currSameFamilyMergeStructureToAlter.comment.replace(commentITToDel, "")
                icescreen_OO.removeCommentToLocusTag2Comment(currSPCompatibleNeigbhor.locusTag, commentITToDel, locusTagMerge2Comment)
                del SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[currSPCompatibleNeigbhor]

            commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and is being associated with same family SPs {}. ".format(
                    currSPCompatibleNeigbhor.locusTag, hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))
            # print(commentITToAdd)

            if commentITToAdd not in currSameFamilyMergeStructure.comment:
                currSameFamilyMergeStructure.comment += commentITToAdd
            icescreen_OO.addCommentToLocusTag2Comment(currSPCompatibleNeigbhor.locusTag, commentITToAdd, locusTagMerge2Comment)
            currSameFamilyMergeStructure.addSPToConjugaisonModule(currSPCompatibleNeigbhor)

        elif len(setGreenLightedSPConjugaisonModuleWithoutFamily) == 2:
            ObjectListGreenLightedSPConjugaisonModuleWithoutFamily = hit.ListSPs()
            ObjectListGreenLightedSPConjugaisonModuleWithoutFamily.list.extend(setGreenLightedSPConjugaisonModuleWithoutFamily)
            # listGreenLightedSPConjugaisonModuleWithoutFamily = list(setGreenLightedSPConjugaisonModuleWithoutFamily)
            # ObjectListGreenLightedSPConjugaisonModuleWithoutFamily.list.sort(key=lambda x: x.locusTag, reverse=True)
            # sort by closer to more distant to a EM structure
            ObjectListGreenLightedSPConjugaisonModuleWithoutFamily.sortListSPsByProximityToMEStructure(currSameFamilyMergeStructure)
            listGreenLightedSPConjugaisonModuleWithoutFamily = ObjectListGreenLightedSPConjugaisonModuleWithoutFamily.list
            firstSPConjugaisonModuleWithoutFamilyInSet = listGreenLightedSPConjugaisonModuleWithoutFamily[0]
            secondSPConjugaisonModuleWithoutFamilyInSet = listGreenLightedSPConjugaisonModuleWithoutFamily[1]

            # like I did for if len(setGreenLightedSPConjugaisonModuleWithoutFamily) == 1:
            if firstSPConjugaisonModuleWithoutFamilyInSet in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and secondSPConjugaisonModuleWithoutFamilyInSet in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # both SP could be match to another same family structure
                commentITToAdd = "Conjugaison module SP {} and {} are compatible neigbhors and could have been associated with SPs {} but were not as they both have other distant structures with same family. ".format(
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs))
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(firstSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)
                icescreen_OO.addCommentToLocusTag2Comment(secondSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)

            elif firstSPConjugaisonModuleWithoutFamilyInSet in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and secondSPConjugaisonModuleWithoutFamilyInSet not in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # firstSPConjugaisonModuleWithoutFamilyInSet could be match to another same family structure, not secondSPConjugaisonModuleWithoutFamilyInSet
                commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and is being associated with same family SPs {}. Please note that there is another compatible neigbhors that could have been associated with same family SPs {}: {} but was not because it can be associated with another same family structure: {}. ".format(
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[firstSPConjugaisonModuleWithoutFamilyInSet].listOrderedSPs))
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(secondSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)
                currSameFamilyMergeStructure.addSPToConjugaisonModule(secondSPConjugaisonModuleWithoutFamilyInSet)

                commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and could have been associated with SPs {} but was not ({} was) because it can also be associated with same family structure {}. ".format(
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[firstSPConjugaisonModuleWithoutFamilyInSet].listOrderedSPs))
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(firstSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)

            elif firstSPConjugaisonModuleWithoutFamilyInSet not in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure and secondSPConjugaisonModuleWithoutFamilyInSet in SPsInSameFamilyMergeStructures2SameFamilyMergeStructure:
                # secondSPConjugaisonModuleWithoutFamilyInSet could be match to another same family structure, not firstSPConjugaisonModuleWithoutFamilyInSet
                commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and is being associated with same family SPs {}. Please note that there is another compatible neigbhors that could have been associated with same family SPs {}: {} but was not because it can be associated with another same family structure: {}. ".format(
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag, hit.ListSPs.GetListProtIdsFromListSP(SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[secondSPConjugaisonModuleWithoutFamilyInSet].listOrderedSPs))
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(firstSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)
                currSameFamilyMergeStructure.addSPToConjugaisonModule(firstSPConjugaisonModuleWithoutFamilyInSet)

                commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and could have been associated with SPs {} but was not ({} was) because it can also be associated with same family structure {}. ".format(
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(SPsInSameFamilyMergeStructures2SameFamilyMergeStructure[secondSPConjugaisonModuleWithoutFamilyInSet].listOrderedSPs))
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(secondSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)

            else:
                # associate the closer one which is firstSPConjugaisonModuleWithoutFamilyInSet
                # both SP could NOT be match to another same family structure, choose firstSPConjugaisonModuleWithoutFamilyInSet over secondSPConjugaisonModuleWithoutFamilyInSet
                commentITToAdd = "Conjugaison module SP {} is the closer compatible neigbhor and is being associated with same family SPs {}. Please note that there is another compatible neigbhors that could have been associated with same family SPs {}: {}. ".format(
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag)
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(firstSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)
                currSameFamilyMergeStructure.addSPToConjugaisonModule(firstSPConjugaisonModuleWithoutFamilyInSet)

                commentITToAdd = "Conjugaison module SP {} is a compatible neigbhor and could have been associated with same family SPs {} but was not ({} was choosen based on proximity). ".format(
                        secondSPConjugaisonModuleWithoutFamilyInSet.locusTag,
                        hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs),
                        firstSPConjugaisonModuleWithoutFamilyInSet.locusTag)
                if commentITToAdd not in currSameFamilyMergeStructure.comment:
                    currSameFamilyMergeStructure.comment += commentITToAdd
                icescreen_OO.addCommentToLocusTag2Comment(secondSPConjugaisonModuleWithoutFamilyInSet.locusTag, commentITToAdd, locusTagMerge2Comment)

        else:
            raise RuntimeError(
                    "Error in buildSameFamilyMergeStructures look around SPs registred for merge for SP greenLightAddSPConjugaisonModule without family: len(setGreenLightedSPConjugaisonModuleWithoutFamily) > 2 ({}) for ICEsIMEsStructures {} and {}".format(
                            str(len(setGreenLightedSPConjugaisonModuleWithoutFamily)),
                            currSameFamilyMergeStructure.internalIdentifier,
                            hit.ListSPs.GetListProtIdsFromListSP(currSameFamilyMergeStructure.listOrderedSPs)))

    # for familyMergeStructuresIT in listSameFamilyMergeStructures:
    #    print("familyMergeStructuresIT {}".format(familyMergeStructuresIT.GetObjectAsJson(True, "")))

    return (listSameFamilyMergeStructures,
            SPsInSameFamilyMergeStructures2SameFamilyMergeStructure)



def checkAndPerformMergeFragmentedForSpecificTypeSP(
    typeOfSPConsideredForFragmentsIT,
    spTypeIT2ICEsIMEsStructure,
    allowCheckingForMultipleDistantRelaxaseIT,
    allowCheckingForMultipleDistantCouplingIT,
    allowCheckingForMultipleDistantVirB4IT,
    listKeysSPTypeIT2ali_start_end_Query_blast,
    groupListSPintoICEsIMEsUsingFamilyInfo,
    maxOverlappingAliLenghtFragmentedSPs,
    locusTagMerge2Comment,
    maxCDSDistanceToMergeFragmentedSPs,
    locusTagIntegrase2Comment,
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance
    ):

    #print("* Starting checkAndPerformMergeFragmentedForSpecificTypeSP for typeOfSPConsideredForFragmentsIT = {}".format(typeOfSPConsideredForFragmentsIT))

    SPTypeIT2setComplementarySPTypeIT = {}
    for (SPTypeIT1IT, SPTypeIT2IT) in list(itertools.combinations(listKeysSPTypeIT2ali_start_end_Query_blast, 2)) :
        if SPTypeIT1IT.Blast_ali_start_Query_blast == -1 or SPTypeIT1IT.Blast_ali_end_Query_blast == -1 or SPTypeIT2IT.Blast_ali_start_Query_blast == -1 or SPTypeIT2IT.Blast_ali_end_Query_blast == -1:
            continue

        # compatible IME_superfamily_of_most_similar_ref_SP and ICE_superfamily_of_most_similar_ref_SP
        greenLightFamilyInfoTest1 = rulesSeedSPExtension.testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
            SPTypeIT1IT,
            spTypeIT2ICEsIMEsStructure[SPTypeIT2IT],
            groupListSPintoICEsIMEsUsingFamilyInfo)
        greenLightFamilyInfoTest2 = rulesSeedSPExtension.testSPFamilyInfoCompatibilityWithICEsIMEsStructure(
            SPTypeIT2IT,
            spTypeIT2ICEsIMEsStructure[SPTypeIT1IT],
            groupListSPintoICEsIMEsUsingFamilyInfo)
        if greenLightFamilyInfoTest1 and greenLightFamilyInfoTest2 :
            # if SPTypeIT : same SPTypeIT_family_domain_of_most_similar_ref_SPFromBlast
            ifSPTypeITSameSPTypeITFamilyDomain = True
            if SPTypeIT1IT.SPType == "Relaxase" and SPTypeIT2IT.SPType == "Relaxase" :
                if SPTypeIT1IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast != SPTypeIT2IT.Relaxase_family_domain_of_most_similar_ref_SPFromBlast :
                    ifSPTypeITSameSPTypeITFamilyDomain = False
            if ifSPTypeITSameSPTypeITFamilyDomain :
                overlapIT = commonMethods.overlap_bounded_intervals(
                    SPTypeIT1IT.Blast_ali_start_Query_blast,
                    SPTypeIT1IT.Blast_ali_end_Query_blast,
                    SPTypeIT2IT.Blast_ali_start_Query_blast,
                    SPTypeIT2IT.Blast_ali_end_Query_blast,
                    SPTypeIT1IT.genomeAccessionRank,
                    SPTypeIT2IT.genomeAccessionRank
                    )
                if overlapIT < maxOverlappingAliLenghtFragmentedSPs :
                    #Those 2 SPs do not overlap much, they are complementary

                    # check that there is no other SP attached to their structure in between the potential parts of the fragmented SP
                    if SPTypeIT1IT.CDSPositionInGenome < SPTypeIT2IT.CDSPositionInGenome :
                        if spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].listOrderedSPs[-1] is not SPTypeIT1IT or spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].listOrderedSPs[0] is not SPTypeIT2IT :
                            commentITToAdd = "{} could complement the other fragment {} to form a single {} but there are other SPs attributed to their structures in between which render their merging impossible, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, typeOfSPConsideredForFragmentsIT)
                            if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment:
                                spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment += commentITToAdd
                            if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment:
                                spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment += commentITToAdd
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            continue
                    else :
                        if spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].listOrderedSPs[0] is not SPTypeIT1IT or spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].listOrderedSPs[-1] is not SPTypeIT2IT :
                            commentITToAdd = "{} could complement the other fragment {} to form a single {} but there are other SPs attributed to their structures in between which render their merging impossible, please manually review. ".format(SPTypeIT1IT.locusTag, SPTypeIT2IT.locusTag, typeOfSPConsideredForFragmentsIT)
                            if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment:
                                spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment += commentITToAdd
                            if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment:
                                spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment += commentITToAdd
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            continue


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
                        commentITToAdd = "{} could complement the fragments {} to form a single {} but {} was found to complement the fragment {} as well, please manually review. ".format(SPTypeIT1IT.locusTag, hit.ListSPs.GetListProtIdsFromSetSP(setComplementarySPTypeIT1IT), typeOfSPConsideredForFragmentsIT, SPTypeIT2IT.locusTag, SPTypeITToAssertForCoherence.locusTag)
                        if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment:
                            spTypeIT2ICEsIMEsStructure[SPTypeIT1IT].comment += commentITToAdd
                        if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment:
                            spTypeIT2ICEsIMEsStructure[SPTypeIT2IT].comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1IT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2IT.locusTag, commentITToAdd, locusTagMerge2Comment)

        if booleanAllComplementarySPTypeITAreCoherent :
            # if more than SP complementary, check that the combination is complementary
            # maxOverlappingAliLenghtFragmentedSPs = 40
            listComplementarySPsToCheckFurther = []
            listComplementarySPsToCheckFurther.append(SPTypeIT1IT)
            listComplementarySPsToCheckFurther.extend(setComplementarySPTypeIT1IT)
            atLeastOneOverlapITGreaterEqualThanMaxAllowed = False
            for (SPTypeIT1ToCheckFurtherIT, SPTypeIT2ToCheckFurtherIT) in list(itertools.combinations(listComplementarySPsToCheckFurther, 2)) :
                if SPTypeIT1ToCheckFurtherIT.Blast_ali_start_Query_blast == -1 or SPTypeIT1ToCheckFurtherIT.Blast_ali_end_Query_blast == -1 or SPTypeIT2ToCheckFurtherIT.Blast_ali_start_Query_blast == -1 or SPTypeIT2ToCheckFurtherIT.Blast_ali_end_Query_blast == -1:
                    raise RuntimeError(
                        "Error in tryMergeFragmentedSPs: SPTypeIT1ToCheckFurtherIT {} or SPTypeIT2ToCheckFurtherIT have Blast_ali__Query_blast  == -1".format(SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag))
                overlapIT = commonMethods.overlap_bounded_intervals(
                    SPTypeIT1ToCheckFurtherIT.Blast_ali_start_Query_blast,
                    SPTypeIT1ToCheckFurtherIT.Blast_ali_end_Query_blast,
                    SPTypeIT2ToCheckFurtherIT.Blast_ali_start_Query_blast,
                    SPTypeIT2ToCheckFurtherIT.Blast_ali_end_Query_blast,
                    SPTypeIT1ToCheckFurtherIT.genomeAccessionRank,
                    SPTypeIT2ToCheckFurtherIT.genomeAccessionRank
                    )
                if overlapIT >= maxOverlappingAliLenghtFragmentedSPs :
                    atLeastOneOverlapITGreaterEqualThanMaxAllowed = True
                    # listComplementarySPsToCheckFurther were checked to complement each other as fragments to form a single SPTypeIT SP, but SPTypeIT1ToCheckFurtherIT was found to overlap with SPTypeIT2ToCheckFurtherIT, please manually review.
                    commentITToAdd = "{} were checked to complement each other as fragments to form a single {}, but {} was found to overlap with {}, please manually review. ".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), typeOfSPConsideredForFragmentsIT, SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag)
                    if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT1ToCheckFurtherIT].comment:
                        spTypeIT2ICEsIMEsStructure[SPTypeIT1ToCheckFurtherIT].comment += commentITToAdd
                    if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT2ToCheckFurtherIT].comment:
                        spTypeIT2ICEsIMEsStructure[SPTypeIT2ToCheckFurtherIT].comment += commentITToAdd
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                    icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)

            if not atLeastOneOverlapITGreaterEqualThanMaxAllowed :
                #OK those listComplementarySPsToCheckFurther can be checked to be merged together
                # checking for maxCDSDistanceToMergeFragmentedSPs = 150
                atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs = False
                for (SPTypeIT1ToCheckFurtherIT, SPTypeIT2ToCheckFurtherIT) in list(itertools.combinations(listComplementarySPsToCheckFurther, 2)) :
                    if abs(SPTypeIT2ToCheckFurtherIT.CDSPositionInGenome - SPTypeIT1ToCheckFurtherIT.CDSPositionInGenome) > maxCDSDistanceToMergeFragmentedSPs :
                        atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs = True
                        # listComplementarySPsToCheckFurther were checked to complement each other as fragments to form a single SPTypeIT SP, but SPTypeIT1ToCheckFurtherIT was found to be further away to SPTypeIT2ToCheckFurtherIT than cutoff, please manually review.
                        commentITToAdd = "{} were checked to complement each other as fragments to form a single {}, but {} was found to be further away to {} than cutoff {}, please manually review. ".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), typeOfSPConsideredForFragmentsIT, SPTypeIT1ToCheckFurtherIT.locusTag, SPTypeIT2ToCheckFurtherIT.locusTag, str(maxCDSDistanceToMergeFragmentedSPs))
                        if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT1ToCheckFurtherIT].comment:
                            spTypeIT2ICEsIMEsStructure[SPTypeIT1ToCheckFurtherIT].comment += commentITToAdd
                        if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPTypeIT2ToCheckFurtherIT].comment:
                            spTypeIT2ICEsIMEsStructure[SPTypeIT2ToCheckFurtherIT].comment += commentITToAdd
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT1ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        icescreen_OO.addCommentToLocusTag2Comment(SPTypeIT2ToCheckFurtherIT.locusTag, commentITToAdd, locusTagMerge2Comment)


                if not atLeastOneCombinationGreaterThenMaxCDSDistanceToMergeFragmentedSPs :
                    #OK those listComplementarySPsToCheckFurther can be merged together

                    firstStructureToMergeBcFragmented = -1
                    secondToLastSetStructuresToMergeBcFragmented = set()
                    for indexInListComplementarySPsToCheckFurther, complementarySPsToCheckFurtherIT in enumerate(listComplementarySPsToCheckFurther):
                        if indexInListComplementarySPsToCheckFurther == 0 :
                            firstStructureToMergeBcFragmented = spTypeIT2ICEsIMEsStructure[complementarySPsToCheckFurtherIT]
                        else :
                            if spTypeIT2ICEsIMEsStructure[complementarySPsToCheckFurtherIT] is not firstStructureToMergeBcFragmented :
                                secondToLastSetStructuresToMergeBcFragmented.add(spTypeIT2ICEsIMEsStructure[complementarySPsToCheckFurtherIT])                    
                    if firstStructureToMergeBcFragmented == -1 :
                        raise RuntimeError(
                            "Error in checkAndPerformMergeFragmentedForSpecificTypeSP: firstStructureToMergeBcFragmented == -1 for {}".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther)))
                    if len(secondToLastSetStructuresToMergeBcFragmented) == 0 :
                        raise RuntimeError(
                            "Error in checkAndPerformMergeFragmentedForSpecificTypeSP: len(secondToLastSetStructuresToMergeBcFragmented) == 0 for {}".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther)))

                    (mostUpstreamICEIMEThatWillHoldTheMerge, listICEIMECompatibleToMerge) = getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge(
                            firstStructureToMergeBcFragmented,
                            secondToLastSetStructuresToMergeBcFragmented,
                            groupListSPintoICEsIMEsUsingFamilyInfo,
                            allowCheckingForMultipleDistantRelaxaseIT,
                            allowCheckingForMultipleDistantCouplingIT,
                            allowCheckingForMultipleDistantVirB4IT)
                    #print("mostUpstreamICEIMEThatWillHoldTheMerge = {}".format(mostUpstreamICEIMEThatWillHoldTheMerge)) # , hit.ListSPs.GetListProtIdsFromListSP(mostUpstreamICEIMEThatWillHoldTheMerge.listOrderedSPs)
                    #print("listICEIMECompatibleToMerge = {}".format(listICEIMECompatibleToMerge))

                    if mostUpstreamICEIMEThatWillHoldTheMerge and len(listICEIMECompatibleToMerge) > 0:
                        # add comment that listComplementarySPsToCheckFurther can be merged together in structure and SP
                        commentITToAdd = "{} complement each other as fragments to form a single {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), typeOfSPConsideredForFragmentsIT)
                        for SPsTthatComplementAsFragmentIT in listComplementarySPsToCheckFurther:
                            icescreen_OO.addCommentToLocusTag2Comment(SPsTthatComplementAsFragmentIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                        if commentITToAdd not in mostUpstreamICEIMEThatWillHoldTheMerge.comment:
                            mostUpstreamICEIMEThatWillHoldTheMerge.comment += commentITToAdd
                        for currICEIMEStructToMerge in listICEIMECompatibleToMerge:
                            if commentITToAdd not in currICEIMEStructToMerge.comment:
                                currICEIMEStructToMerge.comment += commentITToAdd
                            performedAMergeEvent_partOfMultipleMergeEvents = mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                                    currICEIMEStructToMerge,
                                    groupListSPintoICEsIMEsUsingFamilyInfo, locusTagMerge2Comment,
                                    locusTagIntegrase2Comment,
                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                    allowCheckingForMultipleDistantRelaxaseIT,
                                    allowCheckingForMultipleDistantCouplingIT,
                                    allowCheckingForMultipleDistantVirB4IT)
                    else :
                        # add comment that listComplementarySPsToCheckFurther can be merged together in structure and SP
                        commentITToAdd = "{} complement each other as fragments to form a single {} but their structures couldn't be merged, please check. ".format(hit.ListSPs.GetListProtIdsFromListSP(listComplementarySPsToCheckFurther), typeOfSPConsideredForFragmentsIT)
                        for SPsTthatComplementAsFragmentIT in listComplementarySPsToCheckFurther:
                            icescreen_OO.addCommentToLocusTag2Comment(SPsTthatComplementAsFragmentIT.locusTag, commentITToAdd, locusTagMerge2Comment)
                            if commentITToAdd not in spTypeIT2ICEsIMEsStructure[SPsTthatComplementAsFragmentIT].comment:
                                spTypeIT2ICEsIMEsStructure[SPsTthatComplementAsFragmentIT].comment += commentITToAdd
                            
                                    
                                    

def tryMergeFragmentedSPs(
    listICEsIMEsStructures,
    locusTagMerge2Comment,
    locusTagIntegrase2Comment,
    groupListSPintoICEsIMEsUsingFamilyInfo,
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
    maxOverlappingAliLenghtFragmentedSPs,
    maxCDSDistanceToMergeFragmentedSPs
    ):

    #print("**** Starting tryMergeFragmentedSPs")

# - improvment to do : deteccte SP fragmentées. Test 1 segment 16 : WP_158396091.1 should be merged with WP_158396083.1
# Blast_ali_start_Query_blast	Blast_ali_end_Query_blast
# WP_158396083.1 : 1.0	279.0
# WP_158396091.1 : 251.0	432.0
# same protein type
# maxOverlappingAliLenghtFragmentedSPs = 40
# maxCDSDistanceToMergeFragmentedSPs = 150
# compatible IME_superfamily_of_most_similar_ref_SP and ICE_superfamily_of_most_similar_ref_SP
# if relaxase : same Relaxase_family_domain_of_most_similar_ref_SP

    # maxOverlappingAliLenghtFragmentedSPs = 40
    # maxCDSDistanceToMergeFragmentedSPs = 150

    relaxase2ali_start_end_Query_blast = {}
    couplingProtein2ali_start_end_Query_blast = {}
    virB42ali_start_end_Query_blast = {}
    relaxase2ICEsIMEsStructure = {}
    couplingProtein2ICEsIMEsStructure = {}
    virB42ICEsIMEsStructure = {}
    #firstStructureInListICEsIMEsStructures = -1
    #secondToLastStructuresAsSetFromListICEsIMEsStructures = set()

    # for currICEsIMEsStructure in listICEsIMEsStructures:
    for indexInListICEsIMEsStructuresIT, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):

        if currICEsIMEsStructure.delMerging_idxListUpstreamStructure != -1 :
            continue
        # if indexInListICEsIMEsStructuresIT == 0 :
        #     firstStructureInListICEsIMEsStructures = currICEsIMEsStructure
        # else :
        #     secondToLastStructuresAsSetFromListICEsIMEsStructures.add(currICEsIMEsStructure)

        for relaxaseIT in currICEsIMEsStructure.listRelaxase:
            if relaxaseIT in relaxase2ali_start_end_Query_blast :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: relaxaseIT {} in relaxase2ali_start_end_Query_blast".format(relaxaseIT.locusTag))
            relaxase2ali_start_end_Query_blast[relaxaseIT] = (relaxaseIT.Blast_ali_start_Query_blast, relaxaseIT.Blast_ali_end_Query_blast)
            if relaxaseIT in relaxase2ICEsIMEsStructure :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: relaxaseIT {} in relaxase2ICEsIMEsStructure".format(relaxaseIT.locusTag))
            relaxase2ICEsIMEsStructure[relaxaseIT] = currICEsIMEsStructure
        for couplingProteinIT in currICEsIMEsStructure.listCouplingProtein:
            if couplingProteinIT in couplingProtein2ali_start_end_Query_blast :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: couplingProteinIT {} in couplingProtein2ali_start_end_Query_blast".format(couplingProteinIT.locusTag))
            couplingProtein2ali_start_end_Query_blast[couplingProteinIT] = (couplingProteinIT.Blast_ali_start_Query_blast, couplingProteinIT.Blast_ali_end_Query_blast)
            if couplingProteinIT in couplingProtein2ICEsIMEsStructure :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: couplingProteinIT {} in couplingProtein2ICEsIMEsStructure".format(couplingProteinIT.locusTag))
            couplingProtein2ICEsIMEsStructure[couplingProteinIT] = currICEsIMEsStructure
        for virB4IT in currICEsIMEsStructure.listVirB4:
            if virB4IT in virB42ali_start_end_Query_blast :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: virB4IT {} in virB42ali_start_end_Query_blast".format(virB4IT.locusTag))
            virB42ali_start_end_Query_blast[virB4IT] = (virB4IT.Blast_ali_start_Query_blast, virB4IT.Blast_ali_end_Query_blast)
            if virB4IT in virB42ICEsIMEsStructure :
                raise RuntimeError(
                    "Error in tryMergeFragmentedSPs: virB4IT {} in virB42ICEsIMEsStructure".format(virB4IT.locusTag))
            virB42ICEsIMEsStructure[virB4IT] = currICEsIMEsStructure

    # "Relaxase"
    checkAndPerformMergeFragmentedForSpecificTypeSP(
        "Relaxase",
        relaxase2ICEsIMEsStructure,
        True,
        False,
        False,
        list(relaxase2ali_start_end_Query_blast.keys()),
        groupListSPintoICEsIMEsUsingFamilyInfo,
        maxOverlappingAliLenghtFragmentedSPs,
        locusTagMerge2Comment,
        maxCDSDistanceToMergeFragmentedSPs,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance
        )
    # "Coupling protein"
    checkAndPerformMergeFragmentedForSpecificTypeSP(
        "Coupling protein",
        couplingProtein2ICEsIMEsStructure,
        False,
        True,
        False,
        list(couplingProtein2ali_start_end_Query_blast.keys()),
        groupListSPintoICEsIMEsUsingFamilyInfo,
        maxOverlappingAliLenghtFragmentedSPs,
        locusTagMerge2Comment,
        maxCDSDistanceToMergeFragmentedSPs,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance
        )
    # "VirB4"
    checkAndPerformMergeFragmentedForSpecificTypeSP(
        "VirB4",
        virB42ICEsIMEsStructure,
        False,
        False,
        True,
        list(virB42ali_start_end_Query_blast.keys()),
        groupListSPintoICEsIMEsUsingFamilyInfo,
        maxOverlappingAliLenghtFragmentedSPs,
        locusTagMerge2Comment,
        maxCDSDistanceToMergeFragmentedSPs,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance
        )





# this method do not perform the actual merging but test all the merging possibilities between anchor of conj module of similar family within a segment
def tryMergeSameFamilyStructures(
        listICEsIMEsStructures,
        listSameFamilyMergeStructures,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        locusTagMerge2Comment,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):
    for currSameFamilyMergeStructure in listSameFamilyMergeStructures:

        # print(" ********* tryMergeSameFamilyStructures 0: {}".format(currSameFamilyMergeStructure.GetObjectAsJson(True, "")))

        setSPsInCurrSameFamilyMergeStructure = set(currSameFamilyMergeStructure.listOrderedSPs)
        mostUpstreamICEIMEThatWillHoldTheMerge = None
        for currICEsIMEsStructure in listICEsIMEsStructures:

            # setSPsInCurrICEsIMEsStructure = set(currICEsIMEsStructure.listOrderedSPs)  # not listOrderedSPs as we can get obvious integrase that have been added early on
            setSPsInCurrICEsIMEsStructure = set(currICEsIMEsStructure.listRelaxase)
            setSPsInCurrICEsIMEsStructure.update(currICEsIMEsStructure.listCouplingProtein)
            setSPsInCurrICEsIMEsStructure.update(currICEsIMEsStructure.listVirB4)

            # lenInterectIT = len(setSPsInCurrICEsIMEsStructure.intersection(setSPsInCurrSameFamilyMergeStructure))
            # print("\n* EVALUATING: {}".format(currICEsIMEsStructure.GetObjectAsJson(True, "")))
            # print("\n* lenInterectIT: {}".format(str(lenInterectIT)))
            # if lenInterectIT > 0 and lenInterectIT == len(currICEsIMEsStructure.listOrderedSPs):
            if setSPsInCurrICEsIMEsStructure.issubset(setSPsInCurrSameFamilyMergeStructure):

                # found one of the structure to merge
                # print("\n* tryMergeSameFamilyStructures 1 FOUND EQUIVALENT: {}".format(currICEsIMEsStructure.GetObjectAsJson(True, "")))

                if mostUpstreamICEIMEThatWillHoldTheMerge is None:
                    mostUpstreamICEIMEThatWillHoldTheMerge = currICEsIMEsStructure
                else:
                    # print("---> tryMergeSameFamilyStructures 2: {} and {}".format(mostUpstreamICEIMEThatWillHoldTheMerge.internalIdentifier, currICEsIMEsStructure.internalIdentifier))
                    mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                            currICEsIMEsStructure,
                            groupListSPintoICEsIMEsUsingFamilyInfo,
                            locusTagMerge2Comment,
                            locusTagIntegrase2Comment,
                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                            False,
                            False,
                            False)


# this method do not perform the actual merging but test all the merging possibilities between anchor of conj module of different but compatible family within a segment
def tryMergeNestedICEsIMEsStructures(
        listICEsIMEsStructures,
        maxNumberCDSForFilterIMESize,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        locusTagMerge2Comment,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):
    # for currICEsIMEsStructure in listICEsIMEsStructures:
    # assignPutativeTypeStructure(currICEsIMEsStructure, maxNumberCDSForFilterIMESize)
    # score all combination of pair of ICEsIMEsStructures to see if they could be merged
    # dictScoreMerge_2_listPairsICEsIMEsStructure1ToMergeWithICEsIMEsStructure2 = {}

    allowCheckingForMultipleDistantRelaxaseIT = False
    allowCheckingForMultipleDistantCouplingIT = False
    allowCheckingForMultipleDistantVirB4IT = False

    dictICEIME2otherICEIMETomerge2CDSInBetween = {}
    # for currIndexOUTER, currOuterLoopICEsIMEsStructure in enumerate(listICEsIMEsStructures):  # enumerate not working for nested indexes
    for currIndexOUTER in range(len(listICEsIMEsStructures)):  # LOOP_OUTER_COMBI_MERGE
        currOuterLoopICEsIMEsStructure = listICEsIMEsStructures[currIndexOUTER]

        # idxInSeedListMostUpstreamICEsIMEsStructureNested = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostUpstreamICEsIMEsStructureMerged()
        idxInSeedListMostDownstreamICEsIMEsStructureMerged = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostDownstreamICEsIMEsStructureMerged()

        # for currIndexINNER, currInnerLoopICEsIMEsStructure in enumerate(listICEsIMEsStructures):  # LOOP_INNER_COMBI_MERGE
        for currIndexINNER in range(idxInSeedListMostDownstreamICEsIMEsStructureMerged + 1, len(listICEsIMEsStructures)):  # LOOP_INNER_COMBI_MERGE
            currInnerLoopICEsIMEsStructure = listICEsIMEsStructures[currIndexINNER]

            numberOfICEsIMEsStructuresInBetween = currIndexINNER - idxInSeedListMostDownstreamICEsIMEsStructureMerged - 1
            if numberOfICEsIMEsStructuresInBetween < 0:
                raise RuntimeError("Error in tryMergeNestedICEsIMEsStructures: numberOfICEsIMEsStructuresInBetween <= 0 ({}) for ICEsIMEsStructures {} and {}".format(
                        str(numberOfICEsIMEsStructuresInBetween), currOuterLoopICEsIMEsStructure.internalIdentifier, currInnerLoopICEsIMEsStructure.internalIdentifier))
            scoreMerge = scoreMergeTwoICEsIMEsStructures(
                    currOuterLoopICEsIMEsStructure,
                    currInnerLoopICEsIMEsStructure,
                    numberOfICEsIMEsStructuresInBetween,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    allowCheckingForMultipleDistantRelaxaseIT,
                    allowCheckingForMultipleDistantCouplingIT,
                    allowCheckingForMultipleDistantVirB4IT)
            # print("HERE scoreMerge = {} for {} and {}".format(str(scoreMerge), str(currOuterLoopICEsIMEsStructure.internalIdentifier), str(currInnerLoopICEsIMEsStructure.internalIdentifier)))

            if scoreMerge > 0:
                if currOuterLoopICEsIMEsStructure in dictICEIME2otherICEIMETomerge2CDSInBetween:  # key already there
                    dictOtherICEIMETomerge2CDSInBetween = dictICEIME2otherICEIMETomerge2CDSInBetween[currOuterLoopICEsIMEsStructure]
                    dictOtherICEIMETomerge2CDSInBetween[currInnerLoopICEsIMEsStructure] = scoreMerge
                else:  # key not there
                    dictOtherICEIMETomerge2CDSInBetween = {}
                    dictOtherICEIMETomerge2CDSInBetween[currInnerLoopICEsIMEsStructure] = scoreMerge
                    dictICEIME2otherICEIMETomerge2CDSInBetween[currOuterLoopICEsIMEsStructure] = dictOtherICEIMETomerge2CDSInBetween

                if currInnerLoopICEsIMEsStructure in dictICEIME2otherICEIMETomerge2CDSInBetween:  # key already there
                    dictOtherICEIMETomerge2CDSInBetween = dictICEIME2otherICEIMETomerge2CDSInBetween[currInnerLoopICEsIMEsStructure]
                    dictOtherICEIMETomerge2CDSInBetween[currOuterLoopICEsIMEsStructure] = scoreMerge
                else:  # key not there
                    dictOtherICEIMETomerge2CDSInBetween = {}
                    dictOtherICEIMETomerge2CDSInBetween[currOuterLoopICEsIMEsStructure] = scoreMerge
                    dictICEIME2otherICEIMETomerge2CDSInBetween[currInnerLoopICEsIMEsStructure] = dictOtherICEIMETomerge2CDSInBetween

    performedAMergeEvent = False
    for keyICEIMEPrimary in dictICEIME2otherICEIMETomerge2CDSInBetween:

        if keyICEIMEPrimary.delMerging_idxListUpstreamStructure >= 0:
            continue

        # print("keyICEIMEPrimary {} ".format(str(keyICEIMEPrimary.internalIdentifier)))

        dictOtherICEIMETomerge2CDSInBetween = dictICEIME2otherICEIMETomerge2CDSInBetween[keyICEIMEPrimary]
        if len(dictOtherICEIMETomerge2CDSInBetween) == 1:
            for keyICEIMESecondary in dictOtherICEIMETomerge2CDSInBetween:
                if keyICEIMESecondary.delMerging_idxListUpstreamStructure >= 0:
                    continue

                # print("\tkeyICEIMEPrimary {} merge with single keyICEIMESecondary {}".format(str(keyICEIMEPrimary.internalIdentifier), str(keyICEIMESecondary.internalIdentifier)))

                if len(dictICEIME2otherICEIMETomerge2CDSInBetween[keyICEIMESecondary]) == 1:
                    # merge of 1 - 1
                    # print("merge of 1 - 1 with keyICEIMESecondary {}".format(keyICEIMESecondary.internalIdentifier))
                    if keyICEIMEPrimary.listOrderedSPs[0].CDSPositionInGenome < keyICEIMESecondary.listOrderedSPs[0].CDSPositionInGenome:
                        performedAMergeEvent = keyICEIMEPrimary.mergeWith(
                                keyICEIMESecondary,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                locusTagMerge2Comment,
                                locusTagIntegrase2Comment,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                allowCheckingForMultipleDistantRelaxaseIT,
                                allowCheckingForMultipleDistantCouplingIT,
                                allowCheckingForMultipleDistantVirB4IT)
                        # performedAMergeEvent = True
                    else:
                        performedAMergeEvent = keyICEIMESecondary.mergeWith(
                                keyICEIMEPrimary,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                locusTagMerge2Comment,
                                locusTagIntegrase2Comment,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                allowCheckingForMultipleDistantRelaxaseIT,
                                allowCheckingForMultipleDistantCouplingIT,
                                allowCheckingForMultipleDistantVirB4IT)
                        # performedAMergeEvent = True
                else:
                    # print("merge of n - 1 with keyICEIMESecondary {}".format(dictICEIME2otherICEIMETomerge2CDSInBetween[keyICEIMESecondary]))
                    pass
        else:
            # print("merge of 1 - n with dictOtherICEIMETomerge2CDSInBetween {}".format(dictOtherICEIMETomerge2CDSInBetween))
            # multiple putative ICE IME struct to merge, merge first with the more likely
            # 1/ find the merge with identical blast family
            # 2/ among those, find the merge with lower scoreMerge
            dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct = {}
            for keyICEIMESecondary, valueCDSInBetween in dictOtherICEIMETomerge2CDSInBetween.items():

                if keyICEIMESecondary.delMerging_idxListUpstreamStructure >= 0:
                    continue
                # print("assessing {}".format(keyICEIMESecondary.internalIdentifier))

                # print("\tkeyICEIMEPrimary {} merge with multiple keyICEIMESecondary {}".format(str(keyICEIMEPrimary.internalIdentifier), str(keyICEIMESecondary.internalIdentifier)))
                setInteresctFamilies = set()
                setInteresctFamilies.update(keyICEIMEPrimary.setICEFamilyFromBlastOfSPConjModule.intersection(keyICEIMESecondary.setICEFamilyFromBlastOfSPConjModule))
                setInteresctFamilies.update(keyICEIMEPrimary.setIMESuperFamilyFromBlastOfSPConjModule.intersection(keyICEIMESecondary.setIMESuperFamilyFromBlastOfSPConjModule))
                if len(setInteresctFamilies) in dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct:
                    dictCDSInBetween2setICEIMEStruct = dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct[len(setInteresctFamilies)]
                    if valueCDSInBetween in dictCDSInBetween2setICEIMEStruct:
                        setICEIMEStruct = dictCDSInBetween2setICEIMEStruct[valueCDSInBetween]
                        setICEIMEStruct.add(keyICEIMESecondary)
                    else:
                        setICEIMEStruct = set()
                        setICEIMEStruct.add(keyICEIMESecondary)
                        dictCDSInBetween2setICEIMEStruct[valueCDSInBetween] = setICEIMEStruct
                else:
                    dictCDSInBetween2setICEIMEStruct = {}
                    setICEIMEStruct = set()
                    setICEIMEStruct.add(keyICEIMESecondary)
                    dictCDSInBetween2setICEIMEStruct[valueCDSInBetween] = setICEIMEStruct
                    dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct[len(setInteresctFamilies)] = dictCDSInBetween2setICEIMEStruct

            for keyLenSetInteresctFamilies in sorted(dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct.keys(), reverse=True):
                dictCDSInBetween2setICEIMEStruct = dictLenInteresctSetFamilyFromBlastOfSPConjModule2dictCDSInBetween2setICEIMEStruct[keyLenSetInteresctFamilies]
                for keyCDSInBetween in sorted(dictCDSInBetween2setICEIMEStruct.keys(), reverse=False):

                    setICEIMEStruct = dictCDSInBetween2setICEIMEStruct[keyCDSInBetween]

                    if len(setICEIMEStruct) == 1:
                        # found the one to do in priority
                        for currICEIMEStructToMerge in setICEIMEStruct:

                            if len(dictICEIME2otherICEIMETomerge2CDSInBetween[currICEIMEStructToMerge]) == 1:
                                # merge of 1 - 1
                                # print("\tmerge of 1 - 1 with currICEIMEStructToMerge {}".format(currICEIMEStructToMerge.internalIdentifier))
                                if keyICEIMEPrimary.listOrderedSPs[0].CDSPositionInGenome < currICEIMEStructToMerge.listOrderedSPs[0].CDSPositionInGenome:
                                    performedAMergeEvent = keyICEIMEPrimary.mergeWith(
                                            currICEIMEStructToMerge,
                                            groupListSPintoICEsIMEsUsingFamilyInfo,
                                            locusTagMerge2Comment,
                                            locusTagIntegrase2Comment,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                            allowCheckingForMultipleDistantRelaxaseIT,
                                            allowCheckingForMultipleDistantCouplingIT,
                                            allowCheckingForMultipleDistantVirB4IT)
                                    # performedAMergeEvent = True
                                else:
                                    performedAMergeEvent = currICEIMEStructToMerge.mergeWith(
                                            keyICEIMEPrimary,
                                            groupListSPintoICEsIMEsUsingFamilyInfo,
                                            locusTagMerge2Comment, locusTagIntegrase2Comment,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                            allowCheckingForMultipleDistantRelaxaseIT,
                                            allowCheckingForMultipleDistantCouplingIT,
                                            allowCheckingForMultipleDistantVirB4IT)
                                    # performedAMergeEvent = True
                            else:

                                # print("merge of n - 1 with keyICEIMEPrimary {}, currICEIMEStructToMerge {}, len {}:"\
                                #      .format(keyICEIMEPrimary.internalIdentifier, currICEIMEStructToMerge.internalIdentifier, str(len(dictICEIME2otherICEIMETomerge2CDSInBetween[keyICEIMESecondary]))))
                                # for currOtherICEIMETomerge in dictICEIME2otherICEIMETomerge2CDSInBetween[currICEIMEStructToMerge]:
                                #    print("\t{}".format(currOtherICEIMETomerge.internalIdentifier))

                                (mostUpstreamICEIMEThatWillHoldTheMerge,
                                 listICEIMECompatibleToMerge) = getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge(
                                        currICEIMEStructToMerge,
                                        dictICEIME2otherICEIMETomerge2CDSInBetween[currICEIMEStructToMerge],
                                        groupListSPintoICEsIMEsUsingFamilyInfo,
                                        allowCheckingForMultipleDistantRelaxaseIT,
                                        allowCheckingForMultipleDistantCouplingIT,
                                        allowCheckingForMultipleDistantVirB4IT)
                                if mostUpstreamICEIMEThatWillHoldTheMerge and len(listICEIMECompatibleToMerge) > 0:
                                    for currICEIMEStructToMerge in listICEIMECompatibleToMerge:
                                        performedAMergeEvent_partOfMultipleMergeEvents = mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                                                currICEIMEStructToMerge,
                                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                                locusTagMerge2Comment,
                                                locusTagIntegrase2Comment,
                                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                                allowCheckingForMultipleDistantRelaxaseIT,
                                                allowCheckingForMultipleDistantCouplingIT,
                                                allowCheckingForMultipleDistantVirB4IT)
                                        if performedAMergeEvent_partOfMultipleMergeEvents:
                                            performedAMergeEvent = True
                                    # performedAMergeEvent = True
                                else:
                                    strAddToComment = ""
                                    for currOtherICEIMETomerge in dictICEIME2otherICEIMETomerge2CDSInBetween[currICEIMEStructToMerge]:
                                        if strAddToComment:
                                            strAddToComment += ", "
                                        strAddToComment += str(currOtherICEIMETomerge.internalIdentifier)

                                    # print("The merge of this ICE / IME structure {} ({}) with {} was discarded because the latter could be match to multiple incompatible ICE / IME structures: {}. Please manually check which one to merge. "\
                                    #    .format(keyICEIMEPrimary.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(keyICEIMEPrimary.listOrderedSPs), currICEIMEStructToMerge.internalIdentifier, strAddToComment))
                                    commentITToAdd = "The merge of the ICE / IME structure {} ({}) with ICE / IME structure {} was discarded because the latter could be match to multiple incompatible ICE / IME structures: {}. Please manually check which one to merge. ".format(
                                            keyICEIMEPrimary.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(keyICEIMEPrimary.listOrderedSPs), currICEIMEStructToMerge.internalIdentifier, strAddToComment)
                                    if commentITToAdd not in keyICEIMEPrimary.comment:
                                        keyICEIMEPrimary.comment += commentITToAdd
                                    for currSPKeyICEIMEPrimary in keyICEIMEPrimary.listOrderedSPs:
                                        icescreen_OO.addCommentToLocusTag2Comment(currSPKeyICEIMEPrimary.locusTag, commentITToAdd, locusTagMerge2Comment)
                            break

                    else:
                        # there are multiple ICE IME struct that could be merged
                        (mostUpstreamICEIMEThatWillHoldTheMerge, listICEIMECompatibleToMerge) = getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge(
                                keyICEIMEPrimary,
                                setICEIMEStruct,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                allowCheckingForMultipleDistantRelaxaseIT,
                                allowCheckingForMultipleDistantCouplingIT,
                                allowCheckingForMultipleDistantVirB4IT)
                        if mostUpstreamICEIMEThatWillHoldTheMerge and len(listICEIMECompatibleToMerge) > 0:
                            for currICEIMEStructToMerge in listICEIMECompatibleToMerge:
                                performedAMergeEvent_partOfMultipleMergeEvents = mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                                        currICEIMEStructToMerge,
                                        groupListSPintoICEsIMEsUsingFamilyInfo, locusTagMerge2Comment,
                                        locusTagIntegrase2Comment,
                                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                        allowCheckingForMultipleDistantRelaxaseIT,
                                        allowCheckingForMultipleDistantCouplingIT,
                                        allowCheckingForMultipleDistantVirB4IT)
                                if performedAMergeEvent_partOfMultipleMergeEvents:
                                    performedAMergeEvent = True
                            # performedAMergeEvent = True
                        else:
                            commentITToAdd = "The ICE / IME structure {} ({}) could be merged with multiple incompatible ICE / IME structures: {}. Please manually check which one to merge. ".format(
                                    keyICEIMEPrimary.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(keyICEIMEPrimary.listOrderedSPs), EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(setICEIMEStruct))
                            if commentITToAdd not in keyICEIMEPrimary.comment:
                                keyICEIMEPrimary.comment += commentITToAdd
                            for currSPKeyICEIMEPrimary in keyICEIMEPrimary.listOrderedSPs:
                                icescreen_OO.addCommentToLocusTag2Comment(currSPKeyICEIMEPrimary.locusTag, commentITToAdd, locusTagMerge2Comment)

                    break  # just for the lower CDS In Between
                break  # just for the bigger intersect

    if performedAMergeEvent:  # recursive to merge more than 2 ICE /IME
        tryMergeNestedICEsIMEsStructures(
                listICEsIMEsStructures,
                maxNumberCDSForFilterIMESize,
                groupListSPintoICEsIMEsUsingFamilyInfo,
                locusTagMerge2Comment,
                locusTagIntegrase2Comment,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)



# return scoreMerge Int. A score > 0 reflecting that a merge is possible. The score indicate the number of other ICEsIMEsStructures in between the 2 that are assessed. Priority will be given to merge with the lowest score.
def scoreMergeTwoICEsIMEsStructures(
        currOuterLoopICEsIMEsStructure,
        currInnerLoopICEsIMEsStructure,
        numberOfICEsIMEsStructuresInBetween,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        allowCheckingForMultipleDistantRelaxase,
        allowCheckingForMultipleDistantCoupling,
        allowCheckingForMultipleDistantVirB4):

    DEBug_scoreMergeTwoICEsIMEsStructures = False
    if DEBug_scoreMergeTwoICEsIMEsStructures :
        print(" ** starting scoreMergeTwoICEsIMEsStructures with currOuterLoopICEsIMEsStructure {} ({}) and currInnerLoopICEsIMEsStructure {} ({})".format(currOuterLoopICEsIMEsStructure.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(currOuterLoopICEsIMEsStructure.listOrderedSPs), currInnerLoopICEsIMEsStructure.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(currInnerLoopICEsIMEsStructure.listOrderedSPs)))

    if currOuterLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome > currInnerLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:
        raise RuntimeError("Error in scoreMergeTwoICEsIMEsStructures: currOuterLoopICEsIMEsStructure {} is not upstream ({}) of currInnerLoopICEsIMEsStructure {} ({})".format(
                str(currOuterLoopICEsIMEsStructure.internalIdentifier),
                str(currOuterLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome),
                str(currInnerLoopICEsIMEsStructure.internalIdentifier),
                str(currInnerLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome)))

    if currOuterLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : currOuterLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0")
        return -1
    if currInnerLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : currInnerLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0")
        return -1

    if currOuterLoopICEsIMEsStructure.listIntegraseDownstream:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : currOuterLoopICEsIMEsStructure.listIntegraseDownstream")
        return -1  # can not merge a upstream structure with downstream structure if registred downstream integrase
    if currInnerLoopICEsIMEsStructure.listIntegraseUpstream:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : currInnerLoopICEsIMEsStructure.listIntegraseUpstream")
        return -1  # can not merge a downstream structure with upstream structure if registred upstream integrase

    # disallow merge of anchors of ICE family with distant anchor of IME superfamily
    # currOuterLoopICEsIMEsStructure has IME superfamily and no ICE family
    # currInnerLoopICEsIMEsStructure has no IME superfamily and ICE family
    if len(currOuterLoopICEsIMEsStructure.setIMESuperFamilyFromBlastOfSPConjModule) != 0 and len(currOuterLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) == 0 and len(currInnerLoopICEsIMEsStructure.setIMESuperFamilyFromBlastOfSPConjModule) == 0 and len(currInnerLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) != 0:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : can not merge a distant anchors when 1 is related to IME and the other to ICE 1")
        return -1  # can not merge a distant anchors when 1 is related to IME and the other to ICE
    # currOuterLoopICEsIMEsStructure has no IME superfamily and ICE family
    # currInnerLoopICEsIMEsStructure has IME superfamily and no ICE family
    if len(currOuterLoopICEsIMEsStructure.setIMESuperFamilyFromBlastOfSPConjModule) == 0 and len(currOuterLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) != 0 and len(currInnerLoopICEsIMEsStructure.setIMESuperFamilyFromBlastOfSPConjModule) != 0 and len(currInnerLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) == 0:
        if DEBug_scoreMergeTwoICEsIMEsStructures :
            print("return -1 : can not merge a distant anchors when 1 is related to IME and the other to ICE 2")
        return -1  # can not merge a distant anchors when 1 is related to IME and the other to ICE

    # not necessary anymore, as more complex merging now
    # # check that numberOfICEsIMEsStructuresInBetween is different than 0.
    # if numberOfICEsIMEsStructuresInBetween == 0:
    #     if DEBug_scoreMergeTwoICEsIMEsStructures :
    #         print("return -1 : numberOfICEsIMEsStructuresInBetween == 0")
    #     return -1  # fusion with 0 other ICE/IME in between technically possible but do not take into account (it is an artefact if it is fused)

    #can not merge with a structure that is already full, not necessary
    # if len(currOuterLoopICEsIMEsStructure.listRelaxase) >= 1 and len(currOuterLoopICEsIMEsStructure.listCouplingProtein) >= 1 and len(currOuterLoopICEsIMEsStructure.listVirB4) >= 1:
    #     return -1
    # if len(currInnerLoopICEsIMEsStructure.listRelaxase) >= 1 and len(currInnerLoopICEsIMEsStructure.listCouplingProtein) >= 1 and len(currInnerLoopICEsIMEsStructure.listVirB4) >= 1:
    #     return -1


    # can not merge if can not rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure
    for currentSPInner in currInnerLoopICEsIMEsStructure.listOrderedSPs:
        if (currentSPInner.SPType == "Relaxase" or currentSPInner.SPType == "Coupling protein" or currentSPInner.SPType == "VirB4"):
            if len(currentSPInner.setICEsIMEsStructureInConflict) != 0:
                pass  # SP in conflict are regarded as being not there for the merging event
            else:
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currOuterLoopICEsIMEsStructure,
                        currentSPInner,
                        groupListSPintoICEsIMEsUsingFamilyInfo,
                        allowCheckingForMultipleDistantRelaxase,
                        allowCheckingForMultipleDistantCoupling,
                        allowCheckingForMultipleDistantVirB4)
                if not greenLightAddSPConjugaisonModule:
                    if DEBug_scoreMergeTwoICEsIMEsStructures :
                        print("return -1 : not greenLightAddSPConjugaisonModule for currOuterLoopICEsIMEsStructure {} ({}) and currentSPInner {}".format(currOuterLoopICEsIMEsStructure.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(currOuterLoopICEsIMEsStructure.listOrderedSPs), currentSPInner.locusTag))
                    return -1  # because the ICEsIMEsStructures are not compatible for merging

    if DEBug_scoreMergeTwoICEsIMEsStructures :
        print("return numberOfICEsIMEsStructuresInBetween = {}".format(str(numberOfICEsIMEsStructuresInBetween)))
    return numberOfICEsIMEsStructuresInBetween  # because the ICEsIMEsStructures are compatible for merging


# the more downsteam anchor will be merged into the more upstream anchor.
def getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge(
        keyICEIMEPrimarySent,
        setICEIMEStructSent,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        allowCheckingForMultipleDistantRelaxase,
        allowCheckingForMultipleDistantCoupling,
        allowCheckingForMultipleDistantVirB4):

    listICEIMEStruct = []
    listICEIMEStruct.extend(setICEIMEStructSent)
    listICEIMEStruct.append(keyICEIMEPrimarySent)
    CDSPositionInGenome2ICEIMEStruct = {}

    #print("getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge")

    for currICEIMEStructToMerge in listICEIMEStruct:

        #print("\t{}".format(currICEIMEStructToMerge.internalIdentifier))

        if currICEIMEStructToMerge.listOrderedSPs[0].CDSPositionInGenome in CDSPositionInGenome2ICEIMEStruct:
            raise RuntimeError("Error in getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge: currICEIMEStructToMerge.listOrderedSPs[0].CDSPositionInGenome ({}) already in CDSPositionInGenome2ICEIMEStruct".format(
                    str(currICEIMEStructToMerge.listOrderedSPs[0].CDSPositionInGenome)))
        else:
            CDSPositionInGenome2ICEIMEStruct[currICEIMEStructToMerge.listOrderedSPs[0].CDSPositionInGenome] = currICEIMEStructToMerge
    sortedListICEIMEStruct = []
    for currCDSPositionInGenome in sorted(CDSPositionInGenome2ICEIMEStruct.keys(), reverse=False):
        sortedListICEIMEStruct.append(CDSPositionInGenome2ICEIMEStruct[currCDSPositionInGenome])

    allICEIMEStructAreCompatibleToMerge = True
    for currIndexOUTER in range(len(sortedListICEIMEStruct)):  # LOOP_OUTER_COMBI_MERGE
        if not allICEIMEStructAreCompatibleToMerge:
            break
        currOuterLoopICEsIMEsStructure = sortedListICEIMEStruct[currIndexOUTER]
        currOuterLoopICEsIMEsStructure_idxInSeedListMostDownstreamICEsIMEsStructureMerged = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostDownstreamICEsIMEsStructureMerged()
        for currIndexINNER in range(currIndexOUTER + 1, len(sortedListICEIMEStruct)):  # LOOP_INNER_COMBI_MERGE
            currInnerLoopICEsIMEsStructure = sortedListICEIMEStruct[currIndexINNER]
            currInnerLoopICEsIMEsStructure_idxInSeedListMostUpstreamICEsIMEsStructureMerged = currInnerLoopICEsIMEsStructure.findIdxInSeedListMostUpstreamICEsIMEsStructureMerged()
            numberOfICEsIMEsStructuresInBetween = abs(currInnerLoopICEsIMEsStructure_idxInSeedListMostUpstreamICEsIMEsStructureMerged - currOuterLoopICEsIMEsStructure_idxInSeedListMostDownstreamICEsIMEsStructureMerged - 1)
            if numberOfICEsIMEsStructuresInBetween < 0:
                raise RuntimeError("Error in getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge: numberOfICEsIMEsStructuresInBetween <= 0 ({}) for ICEsIMEsStructures {} and {}".format(
                        str(numberOfICEsIMEsStructuresInBetween), currOuterLoopICEsIMEsStructure.internalIdentifier, currInnerLoopICEsIMEsStructure.internalIdentifier))
            scoreMerge = scoreMergeTwoICEsIMEsStructures(
                    currOuterLoopICEsIMEsStructure,
                    currInnerLoopICEsIMEsStructure,
                    numberOfICEsIMEsStructuresInBetween,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    allowCheckingForMultipleDistantRelaxase,
                    allowCheckingForMultipleDistantCoupling,
                    allowCheckingForMultipleDistantVirB4)
            if scoreMerge <= 0:
                #print("scoreMerge <= 0")
                allICEIMEStructAreCompatibleToMerge = False
                break

    if allICEIMEStructAreCompatibleToMerge:
        mostUpstreamICEIMEThatWillHoldTheMerge = None
        listICEIMECompatibleToMerge = []

        for currCDSPositionInGenome in sorted(
                CDSPositionInGenome2ICEIMEStruct.keys(), reverse=False):
            if mostUpstreamICEIMEThatWillHoldTheMerge is None:
                mostUpstreamICEIMEThatWillHoldTheMerge = CDSPositionInGenome2ICEIMEStruct[currCDSPositionInGenome]
            else:
                listICEIMECompatibleToMerge.append(CDSPositionInGenome2ICEIMEStruct[currCDSPositionInGenome])
        return (mostUpstreamICEIMEThatWillHoldTheMerge,
                listICEIMECompatibleToMerge)

    else:
        mostUpstreamICEIMEThatWillHoldTheMerge = None
        listICEIMECompatibleToMerge = []
        return (mostUpstreamICEIMEThatWillHoldTheMerge,
                listICEIMECompatibleToMerge)
