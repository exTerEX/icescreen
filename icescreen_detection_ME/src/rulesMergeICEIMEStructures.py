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

            for familyFromBlastIT in currSP.setSPIMEFamilyFromBlast:
                # print("setSPIMEFamilyFromBlast: {} for {}".format(
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
                        groupListSPintoICEsIMEsUsingFamilyInfo)
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
                        groupListSPintoICEsIMEsUsingFamilyInfo)
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
                        groupListSPintoICEsIMEsUsingFamilyInfo)
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
                # if len(currUpstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currUpstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currUpstreamSPIT.setSPIMEFamilyFromBlast) == 0:
                #    preconditionsCompatibleNeigbhorFamilyStatus = True
                preconditionsCompatibleNeigbhorFamilyStatus = True
                if len(currSameFamilyMergeStructure.listVirB4) > 0 and len(currUpstreamSPIT.setSPIMEFamilyFromBlast) > 0:
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
                                groupListSPintoICEsIMEsUsingFamilyInfo)
                        if greenLightAddSPConjugaisonModule:
                            setGreenLightedSPConjugaisonModuleWithoutFamily.add(currUpstreamSPIT)
                            # print("added currUpstreamSPIT "+currUpstreamSPIT.locusTag+"to setGreenLightedSPConjugaisonModuleWithoutFamily")

            # deal with downstream SP
            idxDownstreamSPIT = currSPameFamilyMergeStructure.idxInListSP + 1
            if idxDownstreamSPIT < len(currListSPs):
                currDownstreamSPIT = currListSPs[idxDownstreamSPIT]

                # preconditionsCompatibleNeigbhorFamilyStatus = False
                # if len(currDownstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPIMEFamilyFromBlast) == 0:
                #    preconditionsCompatibleNeigbhorFamilyStatus = True
                preconditionsCompatibleNeigbhorFamilyStatus = True
                if len(currSameFamilyMergeStructure.listVirB4) > 0 and len(currDownstreamSPIT.setSPIMEFamilyFromBlast) > 0:
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
                    # if len(currDownstreamSPIT.setSPICESuperFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPICEFamilyFromBlast) == 0 and len(currDownstreamSPIT.setSPIMEFamilyFromBlast) == 0:
                    if currDownstreamSPIT not in currSameFamilyMergeStructure.listOrderedSPs:
                        greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                                currSameFamilyMergeStructure,
                                currDownstreamSPIT,
                                groupListSPintoICEsIMEsUsingFamilyInfo)
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
                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)


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
    dictICEIME2otherICEIMETomerge2CDSInBetween = {}
    # for currIndexOUTER, currOuterLoopICEsIMEsStructure in enumerate(listICEsIMEsStructures):  # enumerate not working for nested indexes
    for currIndexOUTER in range(len(listICEsIMEsStructures)):  # LOOP_OUTER_COMBI_MERGE
        currOuterLoopICEsIMEsStructure = listICEsIMEsStructures[currIndexOUTER]

        # idxInSeedListMostUpstreamICEsIMEsStructureNested = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostUpstreamICEsIMEsStructureNested()
        idxInSeedListMostDownstreamICEsIMEsStructureNested = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostDownstreamICEsIMEsStructureNested()

        # for currIndexINNER, currInnerLoopICEsIMEsStructure in enumerate(listICEsIMEsStructures):  # LOOP_INNER_COMBI_MERGE
        for currIndexINNER in range(idxInSeedListMostDownstreamICEsIMEsStructureNested + 1, len(listICEsIMEsStructures)):  # LOOP_INNER_COMBI_MERGE
            currInnerLoopICEsIMEsStructure = listICEsIMEsStructures[currIndexINNER]

            numberOfICEsIMEsStructuresInBetween = currIndexINNER - idxInSeedListMostDownstreamICEsIMEsStructureNested - 1
            if numberOfICEsIMEsStructuresInBetween < 0:
                raise RuntimeError("Error in tryMergeNestedICEsIMEsStructures: numberOfICEsIMEsStructuresInBetween <= 0 ({}) for ICEsIMEsStructures {} and {}".format(
                        str(numberOfICEsIMEsStructuresInBetween), currOuterLoopICEsIMEsStructure.internalIdentifier, currInnerLoopICEsIMEsStructure.internalIdentifier))
            scoreMerge = scoreMergeTwoICEsIMEsStructures(
                    currOuterLoopICEsIMEsStructure,
                    currInnerLoopICEsIMEsStructure,
                    numberOfICEsIMEsStructuresInBetween,
                    groupListSPintoICEsIMEsUsingFamilyInfo)
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
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
                        # performedAMergeEvent = True
                    else:
                        performedAMergeEvent = keyICEIMESecondary.mergeWith(
                                keyICEIMEPrimary,
                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                locusTagMerge2Comment,
                                locusTagIntegrase2Comment,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
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
                setInteresctFamilies.update(keyICEIMEPrimary.setIMEFamilyFromBlastOfSPConjModule.intersection(keyICEIMESecondary.setIMEFamilyFromBlastOfSPConjModule))
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
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
                                    # performedAMergeEvent = True
                                else:
                                    performedAMergeEvent = currICEIMEStructToMerge.mergeWith(
                                            keyICEIMEPrimary,
                                            groupListSPintoICEsIMEsUsingFamilyInfo,
                                            locusTagMerge2Comment, locusTagIntegrase2Comment,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
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
                                        groupListSPintoICEsIMEsUsingFamilyInfo)
                                if mostUpstreamICEIMEThatWillHoldTheMerge and len(listICEIMECompatibleToMerge) > 0:
                                    for currICEIMEStructToMerge in listICEIMECompatibleToMerge:
                                        performedAMergeEvent_partOfMultipleMergeEvents = mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                                                currICEIMEStructToMerge,
                                                groupListSPintoICEsIMEsUsingFamilyInfo,
                                                locusTagMerge2Comment,
                                                locusTagIntegrase2Comment,
                                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                                useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
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
                                keyICEIMEPrimary, setICEIMEStruct, groupListSPintoICEsIMEsUsingFamilyInfo)
                        if mostUpstreamICEIMEThatWillHoldTheMerge and len(listICEIMECompatibleToMerge) > 0:
                            for currICEIMEStructToMerge in listICEIMECompatibleToMerge:
                                performedAMergeEvent_partOfMultipleMergeEvents = mostUpstreamICEIMEThatWillHoldTheMerge.mergeWith(
                                        currICEIMEStructToMerge,
                                        groupListSPintoICEsIMEsUsingFamilyInfo, locusTagMerge2Comment,
                                        locusTagIntegrase2Comment,
                                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
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
        groupListSPintoICEsIMEsUsingFamilyInfo):

    if currOuterLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome > currInnerLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:
        raise RuntimeError("Error in scoreMergeTwoICEsIMEsStructures: currOuterLoopICEsIMEsStructure {} is not upstream ({}) of currInnerLoopICEsIMEsStructure {} ({})".format(
                str(currOuterLoopICEsIMEsStructure.internalIdentifier),
                str(currOuterLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome),
                str(currInnerLoopICEsIMEsStructure.internalIdentifier),
                str(currInnerLoopICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome)))

    if currOuterLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
        return -1
    if currInnerLoopICEsIMEsStructure.delMerging_idxListUpstreamStructure >= 0:
        return -1

    if len(currOuterLoopICEsIMEsStructure.listRelaxase) >= 1 and len(currOuterLoopICEsIMEsStructure.listCouplingProtein) >= 1 and len(currOuterLoopICEsIMEsStructure.listVirB4) >= 1:
        return -1

    if len(currInnerLoopICEsIMEsStructure.listRelaxase) >= 1 and len(currInnerLoopICEsIMEsStructure.listCouplingProtein) >= 1 and len(currInnerLoopICEsIMEsStructure.listVirB4) >= 1:
        return -1

    if currOuterLoopICEsIMEsStructure.listIntegraseDownstream:
        return -1  # can not merge a upstream structure with downstream structure if registred downstream integrase
    if currInnerLoopICEsIMEsStructure.listIntegraseUpstream:
        return -1  # can not merge a downstream structure with upstream structure if registred upstream integrase

    # disallow merge of anchors of ICE family with distant anchor of IME family
    # currOuterLoopICEsIMEsStructure has IME family and no ICE family
    # currInnerLoopICEsIMEsStructure has no IME family and ICE family
    if len(currOuterLoopICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule) != 0 and len(currOuterLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) == 0 and len(currInnerLoopICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule) == 0 and len(currInnerLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) != 0:
        return -1  # can not merge a distant anchors when 1 is related to IME and the other to ICE
    # currOuterLoopICEsIMEsStructure has no IME family and ICE family
    # currInnerLoopICEsIMEsStructure has IME family and no ICE family
    if len(currOuterLoopICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule) == 0 and len(currOuterLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) != 0 and len(currInnerLoopICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule) != 0 and len(currInnerLoopICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) == 0:
        return -1  # can not merge a distant anchors when 1 is related to IME and the other to ICE

    for currentSPInner in currInnerLoopICEsIMEsStructure.listOrderedSPs:
        if (currentSPInner.SPType == "Relaxase" or currentSPInner.SPType == "Coupling protein" or currentSPInner.SPType == "VirB4"):
            if len(currentSPInner.setICEsIMEsStructureInConflict) != 0:
                pass  # SP in conflict are regarded as being not there for the merging event
            else:
                greenLightAddSPConjugaisonModule = rulesSeedSPExtension.tryAddingSPToConjugaisonModuleEMStructure(
                        currOuterLoopICEsIMEsStructure,
                        currentSPInner,
                        groupListSPintoICEsIMEsUsingFamilyInfo)
                if not greenLightAddSPConjugaisonModule:
                    return -1  # because the ICEsIMEsStructures are not compatible for merging
    # check that numberOfICEsIMEsStructuresInBetween is different than 0. 0 would mean that the ListSPs.seedICEsIMEsStructure method failed
    if numberOfICEsIMEsStructuresInBetween == 0:
        return -1  # fusion with 0 other ICE/IME in between technically possible but do not take into account (it is an artefact if it is fused)
        # raise RuntimeError("Error in scoreMergeTwoICEsIMEsStructures: structure could be merge but numberOfICEsIMEsStructuresInBetween == 0 for structures {} and {}".format(str(currOuterLoopICEsIMEsStructure.internalIdentifier), str(currInnerLoopICEsIMEsStructure.internalIdentifier)))
    return numberOfICEsIMEsStructuresInBetween  # because the ICEsIMEsStructures are compatible for merging


# the more downsteam anchor will be merged into the more upstream anchor.
def getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge(
        keyICEIMEPrimarySent,
        setICEIMEStructSent,
        groupListSPintoICEsIMEsUsingFamilyInfo):

    listICEIMEStruct = []
    listICEIMEStruct.extend(setICEIMEStructSent)
    listICEIMEStruct.append(keyICEIMEPrimarySent)
    CDSPositionInGenome2ICEIMEStruct = {}

    # print("getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge")

    for currICEIMEStructToMerge in listICEIMEStruct:

        # print("\t{}".format(currICEIMEStructToMerge.internalIdentifier))

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
        currOuterLoopICEsIMEsStructure_idxInSeedListMostDownstreamICEsIMEsStructureNested = currOuterLoopICEsIMEsStructure.findIdxInSeedListMostDownstreamICEsIMEsStructureNested()
        for currIndexINNER in range(currIndexOUTER + 1, len(sortedListICEIMEStruct)):  # LOOP_INNER_COMBI_MERGE
            currInnerLoopICEsIMEsStructure = sortedListICEIMEStruct[currIndexINNER]
            currInnerLoopICEsIMEsStructure_idxInSeedListMostUpstreamICEsIMEsStructureNested = currInnerLoopICEsIMEsStructure.findIdxInSeedListMostUpstreamICEsIMEsStructureNested()
            numberOfICEsIMEsStructuresInBetween = abs(currInnerLoopICEsIMEsStructure_idxInSeedListMostUpstreamICEsIMEsStructureNested - currOuterLoopICEsIMEsStructure_idxInSeedListMostDownstreamICEsIMEsStructureNested - 1)
            if numberOfICEsIMEsStructuresInBetween < 0:
                raise RuntimeError("Error in getmostUpstreamICEIMEThatWillHoldTheMergeAndlistSecondaryICEIMECompatibleToMerge: numberOfICEsIMEsStructuresInBetween <= 0 ({}) for ICEsIMEsStructures {} and {}".format(
                        str(numberOfICEsIMEsStructuresInBetween), currOuterLoopICEsIMEsStructure.internalIdentifier, currInnerLoopICEsIMEsStructure.internalIdentifier))
            scoreMerge = scoreMergeTwoICEsIMEsStructures(
                    currOuterLoopICEsIMEsStructure,
                    currInnerLoopICEsIMEsStructure,
                    numberOfICEsIMEsStructuresInBetween,
                    groupListSPintoICEsIMEsUsingFamilyInfo)
            if scoreMerge <= 0:
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
