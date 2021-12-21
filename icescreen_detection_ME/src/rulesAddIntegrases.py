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
import icescreen_OO
import EMStructure
import hit

########
# GLOBAL VARS  #
########

DEBUG = False
structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure = {}
dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd = {}


# noqa: E501 This method correct the wrongly attribution of a subsequent integrase that was initially thought to be an obvious choice for a conjugation module
def changeObviousIntegraseAttributionToUnsureBecauseOfBothUpstreamAndDownstreamEqualyPossible(
        currICEsIMEsStructure,
        otherICEsIMEsStructureToMerge,
        locusTagIntegrase2Comment):

    listUpstreamIntChanged = []
    listDownstreamIntChanged = []

    strCommentIT = "Following a merge, upstream integrase {} and downstream integrase {} appear to be both adjacent and have already been seen rattached to the conjugaison module family of structure {}. ".format(
            hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseUpstream),
            hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseDownstream),
            currICEsIMEsStructure.internalIdentifier)
    if strCommentIT not in currICEsIMEsStructure.comment:
        currICEsIMEsStructure.comment += strCommentIT
    for currSp in currICEsIMEsStructure.listIntegraseUpstream:
        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
        # take off key and clean up integrase and comment in already registred structure
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
        if strCommentToRemove in currICEsIMEsStructure.comment:
            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, otherICEsIMEsStructureToMerge.internalIdentifier)
        if strCommentToRemove in otherICEsIMEsStructureToMerge.comment:
            otherICEsIMEsStructureToMerge.comment = otherICEsIMEsStructureToMerge.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)

        # currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)
        listUpstreamIntChanged.append(currSp)
        if currSp in currICEsIMEsStructure.listOrderedSPs:
            currICEsIMEsStructure.listOrderedSPs.remove(currSp)
    currICEsIMEsStructure.listIntegraseUpstream.clear()

    for currSp in currICEsIMEsStructure.listIntegraseDownstream:
        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
        # take off key and clean up integrase and comment in already registred structure
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
        if strCommentToRemove in currICEsIMEsStructure.comment:
            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, otherICEsIMEsStructureToMerge.internalIdentifier)
        if strCommentToRemove in otherICEsIMEsStructureToMerge.comment:
            otherICEsIMEsStructureToMerge.comment = otherICEsIMEsStructureToMerge.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)

        # currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)
        listDownstreamIntChanged.append(currSp)
        if currSp in currICEsIMEsStructure.listOrderedSPs:
            currICEsIMEsStructure.listOrderedSPs.remove(currSp)
    currICEsIMEsStructure.listIntegraseDownstream.clear()

    return (listUpstreamIntChanged, listDownstreamIntChanged)


# valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase returned:
# 0: did not perform the test
# 1: Downstream integrase is significantly further away to ICE/IME structure
# 2: Downstream integrase is NOT significantly further away to ICE/IME structure
# 3: Upstream integrase is significantly further away to ICE/IME structure
# 4: Upstream integrase is NOT significantly further away to ICE/IME structure
def useIntegraseCDSDistanceToChooseBetweenUpstreamAndDownstreamStructures(
        currSp,
        EMStructurePreviouslyAdded,
        EMStructureToAdd,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
        booleanCommentAndCleanUpIfTestIsSignificant,
        locusTagIntegrase2Comment):

    if useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance < 0 or useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance < 0:
        return (0, -1, -1)
    else:
        nextUpstreamSpThatIsNotAnIntegrase = EMStructurePreviouslyAdded.findMostDownstrSpNotIntegr()
        nextDownstreamSpThatIsNotAnIntegrase = EMStructureToAdd.findMostUpstrSpNotIntegr()
        distCDSWithUpstreamStructure = abs(currSp.CDSPositionInGenome - nextUpstreamSpThatIsNotAnIntegrase.CDSPositionInGenome)
        distCDSWithDownstreamStructure = abs(nextDownstreamSpThatIsNotAnIntegrase.CDSPositionInGenome - currSp.CDSPositionInGenome)
        if distCDSWithUpstreamStructure < distCDSWithDownstreamStructure:
            # upstream structure is closer
            if distCDSWithUpstreamStructure <= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance and distCDSWithDownstreamStructure >= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance:
                # Downstream structure is significantly further away to CDS
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Integrase {} is significantly further away to downstream ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than upstream ICE/IME structure {} ({} CDSs which is <= to lower cutoff {} CDSs), the downstream integrase is therfore not attributed to ICE/IME structure {}. ".format(
                            currSp.locusTag,
                            str(EMStructureToAdd.internalIdentifier),
                            str(distCDSWithDownstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            str(EMStructurePreviouslyAdded.internalIdentifier),
                            str(distCDSWithUpstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance),
                            str(EMStructureToAdd.internalIdentifier))
                    if strCommentIT not in EMStructurePreviouslyAdded.comment:
                        EMStructurePreviouslyAdded.comment += strCommentIT
                    if strCommentIT not in EMStructureToAdd.comment:
                        EMStructureToAdd.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    EMStructurePreviouslyAdded.listIntegraseDownstream.remove(currSp)
                    EMStructurePreviouslyAdded.listOrderedSPs.remove(currSp)
                    EMStructureToAdd.listIntegraseUpstream.append(currSp)
                    EMStructureToAdd.listOrderedSPs.append(currSp)
                return (1, distCDSWithDownstreamStructure, distCDSWithUpstreamStructure)
            else:
                # Downstream structure is NOT significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Integrase {} is NOT significantly further away to downstream ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than upstream ICE/IME structure {} ({} CDSs which is <= to lower cutoff {} CDSs). ".format(
                            currSp.locusTag,
                            str(EMStructureToAdd.internalIdentifier),
                            str(distCDSWithDownstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            str(EMStructurePreviouslyAdded.internalIdentifier),
                            str(distCDSWithUpstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance))
                    if strCommentIT not in EMStructurePreviouslyAdded.comment:
                        EMStructurePreviouslyAdded.comment += strCommentIT
                    if strCommentIT not in EMStructureToAdd.comment:
                        EMStructureToAdd.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    EMStructurePreviouslyAdded.listIntegraseDownstream.remove(currSp)
                    EMStructurePreviouslyAdded.listOrderedSPs.remove(currSp)
                    EMStructurePreviouslyAdded.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)
                    EMStructureToAdd.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)
                return (2, distCDSWithDownstreamStructure, distCDSWithUpstreamStructure)
        else:
            # downstream structure is closer
            if distCDSWithDownstreamStructure <= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance and distCDSWithUpstreamStructure >= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance:
                # Upstream structure is significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Integrase {} is significantly further away to upstream ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than downstream ICE/IME structure {} ({} CDSs which is <= to lower cutoff {} CDSs), the upstream integrase is therfore not attributed to ICE/IME structure {}. ".format(
                            currSp.locusTag,
                            str(EMStructurePreviouslyAdded.internalIdentifier),
                            str(distCDSWithUpstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            str(EMStructureToAdd.internalIdentifier),
                            str(distCDSWithDownstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance),
                            str(EMStructurePreviouslyAdded.internalIdentifier))
                    if strCommentIT not in EMStructurePreviouslyAdded.comment:
                        EMStructurePreviouslyAdded.comment += strCommentIT
                    if strCommentIT not in EMStructureToAdd.comment:
                        EMStructureToAdd.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                return (3, distCDSWithDownstreamStructure, distCDSWithUpstreamStructure)
            else:
                # Upstream structure is NOT significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Integrase {} is NOT significantly further away to upstream ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than downstream ICE/IME structure {} ({} CDSs which is <= to lower cutoff {} CDSs). ".format(
                            currSp.locusTag,
                            str(EMStructurePreviouslyAdded.internalIdentifier),
                            str(distCDSWithUpstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            str(EMStructureToAdd.internalIdentifier),
                            str(distCDSWithDownstreamStructure),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance))
                    if strCommentIT not in EMStructurePreviouslyAdded.comment:
                        EMStructurePreviouslyAdded.comment += strCommentIT
                    if strCommentIT not in EMStructureToAdd.comment:
                        EMStructureToAdd.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                return (4, distCDSWithDownstreamStructure, distCDSWithUpstreamStructure)


# valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase returned:
# 0: did not perform the test
# 1: Downstream integrase is significantly further away to ICE/IME structure
# 2: Downstream integrase is NOT significantly further away to ICE/IME structure
# 3: Upstream integrase is significantly further away to ICE/IME structure
# 4: Upstream integrase is NOT significantly further away to ICE/IME structure
def useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase(
        listUpstreamIntToAdd,
        listDownstreamIntToAdd,
        currICEsIMEsStructure,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
        booleanCommentAndCleanUpIfTestIsSignificant,
        locusTagIntegrase2Comment):

    if useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance < 0 or useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance < 0:
        return (0, -1, -1)
    else:
        distCDSWithUpstreamIntegrase = abs(currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome - listUpstreamIntToAdd[-1].CDSPositionInGenome)
        distCDSWithDownstreamIntegrase = abs(listDownstreamIntToAdd[0].CDSPositionInGenome - currICEsIMEsStructure.listOrderedSPs[-1].CDSPositionInGenome)
        if distCDSWithUpstreamIntegrase < distCDSWithDownstreamIntegrase:
            # upstream integrase is closer
            if distCDSWithUpstreamIntegrase <= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance and distCDSWithDownstreamIntegrase >= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance:
                # Downstream integrase is significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    # disregard downstream integrase
                    strCommentIT = "Downstream integrase {} is significantly further away to ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than upstream integrase {} ({} CDSs which is <= to lower cutoff {} CDSs), the downstream integrase is therfore not attributed to ICE/IME structure {}. ".format(
                            hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd),
                            str(currICEsIMEsStructure.internalIdentifier),
                            str(distCDSWithDownstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd),
                            str(distCDSWithUpstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance),
                            str(currICEsIMEsStructure.internalIdentifier))
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    for currSp in listDownstreamIntToAdd:
                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    listDownstreamIntToAdd.clear()
                return (1, distCDSWithDownstreamIntegrase, distCDSWithUpstreamIntegrase)
            else:
                # Downstream integrase is NOT significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Downstream integrase {} is NOT significantly further away to ICE/IME structure {} ({} CDSs, higher cutoff is {} CDSs) than upstream integrase {} ({} CDSs, lower cutoff is {} CDSs). ".format(
                            hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd),
                            str(currICEsIMEsStructure.internalIdentifier),
                            str(distCDSWithDownstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd),
                            str(distCDSWithUpstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance))
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    for currSp in listDownstreamIntToAdd:
                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                return (2, distCDSWithDownstreamIntegrase, distCDSWithUpstreamIntegrase)
        else:
            # downstream integrase is closer
            if distCDSWithDownstreamIntegrase <= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance and distCDSWithUpstreamIntegrase >= useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance:
                # Upstream integrase is significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Upstream integrase {} is significantly further away to ICE/IME structure {} ({} CDSs which is >= to higher cutoff {} CDSs) than downstream integrase {} ({} CDSs which is <= to lower cutoff {} CDSs), the upstream integrase is therfore not attributed to ICE/IME structure {}. ".format(
                            hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd),
                            str(currICEsIMEsStructure.internalIdentifier),
                            str(distCDSWithUpstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd),
                            str(distCDSWithDownstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance),
                            str(currICEsIMEsStructure.internalIdentifier))
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    for currSp in listUpstreamIntToAdd:
                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    listUpstreamIntToAdd.clear()
                return (3, distCDSWithDownstreamIntegrase, distCDSWithUpstreamIntegrase)
            else:
                # Upstream integrase is NOT significantly further away to ICE/IME structure
                if booleanCommentAndCleanUpIfTestIsSignificant:
                    strCommentIT = "Upstream integrase {} is NOT significantly further away to ICE/IME structure {} ({} CDSs, higher cutoff is {} CDSs) than downstream integrase {} ({} CDSs, lower cutoff is {} CDSs). ".format(
                            hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd),
                            str(currICEsIMEsStructure.internalIdentifier), str(distCDSWithUpstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance),
                            hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd),
                            str(distCDSWithDownstreamIntegrase),
                            str(useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance))
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    for currSp in listUpstreamIntToAdd:
                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                return (4, distCDSWithDownstreamIntegrase, distCDSWithUpstreamIntegrase)


def checkIfStructureHasBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structureAndIsNotInSetLocusTagToNotConsiderAsManuallyCheck(
        ICEsIMEsStructureSent,
        setLocusTagToNotConsiderAseManuallyCheckSent,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
        locusTagIntegrase2Comment):
    booleanToReturn = False

    '''
    print("\n\n0 checkIfStructureHasBothUpstreamAn... $: ICEsIMEsStructureSent {}".format(ICEsIMEsStructureSent.internalIdentifier))
    for structureIT2 in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
        print("0.1 checkIfStructureHasBothUpstreamAn... $: structureIT2.internalIdentifier {}".format(structureIT2.internalIdentifier))
    '''

    if ICEsIMEsStructureSent in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:

        listOfListIntegraseUpstreamDownstreamIT = structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[ICEsIMEsStructureSent]
        setIntegraseUpstream = set()
        setIntegraseDownstream = set()
        for integraseUpstreamIT in listOfListIntegraseUpstreamDownstreamIT[0]:
            # print("checkIfStructureHasBothUpstreamAn...: integraseUpstreamIT.locusTag {} ; setLocusTagToNotConsiderAseManuallyCheckSent {} " \
            # .format(integraseUpstreamIT.locusTag, repr(setLocusTagToNotConsiderAseManuallyCheckSent)))

            if integraseUpstreamIT.locusTag not in setLocusTagToNotConsiderAseManuallyCheckSent and integraseUpstreamIT not in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                setIntegraseUpstream.add(integraseUpstreamIT)

        for integraseDownstreamIT in listOfListIntegraseUpstreamDownstreamIT[1]:
            # print("checkIfStructureHasBothUpstreamAn...: integraseUpstreamIT.locusTag {} ; setLocusTagToNotConsiderAseManuallyCheckSent {} " \
            # .format(integraseUpstreamIT.locusTag, repr(setLocusTagToNotConsiderAseManuallyCheckSent)))
            if integraseDownstreamIT.locusTag not in setLocusTagToNotConsiderAseManuallyCheckSent and integraseDownstreamIT not in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                setIntegraseDownstream.add(integraseDownstreamIT)
                # print("2 checkIfStructureHasBothUpstreamAn... $: len(setIntegraseUpstream) {} ; len(setIntegraseDownstream) {}" \
                # .format(len(setIntegraseUpstream), len(setIntegraseDownstream)))

        if len(setIntegraseUpstream) > 0 and len(setIntegraseDownstream) > 0:
            # check significant distance here too
            # valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase returned:
            # 0: did not perform the test
            # 1: Downstream integrase is significantly further away to ICE/IME structure
            # 2: Downstream integrase is NOT significantly further away to ICE/IME structure
            # 3: Upstream integrase is significantly further away to ICE/IME structure
            # 4: Upstream integrase is NOT significantly further away to ICE/IME structure
            (valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase,
             distCDSWithDownstreamIntegrase,
             distCDSWithUpstreamIntegrase) = useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase(
                     listOfListIntegraseUpstreamDownstreamIT[0],
                     listOfListIntegraseUpstreamDownstreamIT[1],
                     ICEsIMEsStructureSent,
                     useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                     useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                     False,
                     locusTagIntegrase2Comment)

            if valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase == 2:
                booleanToReturn = True
            elif valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase == 4:
                booleanToReturn = True
#    print("1 checkIfStructureHasBothUpstreamAn...: {} ; ICEsIMEsStructureSent {} ; setLocusTagToNotConsiderAseManuallyCheckSent {} " \
#          .format(booleanToReturn, ICEsIMEsStructureSent.internalIdentifier, repr(setLocusTagToNotConsiderAseManuallyCheckSent)))

    return booleanToReturn


# integrase can be associated to any conjugation module without restriction on families. This method just add a comment if the integrase has already be seen associated to a conjugation module but now that the family for integrase has been removed from the database, this method is not usefull anymore.
def markIntegraseNeverPreviouslySeenAssociatedToConjugationModuleFamilyAsToManuallyCheck(
        listUpstreamIntToAdd,
        listDownstreamIntToAdd,
        currICEsIMEsStructure,
        locusTagIntegrase2Comment,
        setLocusTagToNotConsiderAseManuallyCheckSent,
        ICEsIMEsStructure2setIntegraseManuallyCheck2comment):

    keepIntegraseAssignAsSureButAddCommentAboutFamily = True
    # if len(listUpstreamIntToAdd) > 0 and len(listDownstreamIntToAdd) == 0 and currICEsIMEsStructure not in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
    #    keepIntegraseAssignAsSureButAddCommentAboutFamily = True
    # elif len(listUpstreamIntToAdd) == 0 and len(listDownstreamIntToAdd) > 0 and currICEsIMEsStructure not in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
    #    keepIntegraseAssignAsSureButAddCommentAboutFamily = True

    listLocusTagInCommentNoFamilyMatch = []
    setIntegraseToManuallyCheckNoFamilyMatch = set()
    if listUpstreamIntToAdd:
        for currSp in reversed(listUpstreamIntToAdd):
            # slight change, for the best ?
            if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 0:
                if keepIntegraseAssignAsSureButAddCommentAboutFamily:
                    # if change text, change below and also in EMStructure.py around l 761
                    commentITToAdd = "The integrase {} has been previously seen associated with the conjugation module superfamily {} but not yet with the superfamily {} of structure {}. ".format(
                            currSp.locusTag,
                            ", ".join(str(i) for i in sorted(currSp.setSPICESuperFamilyFromBlast)),
                            ", ".join(str(i) for i in sorted(currICEsIMEsStructure.setICESuperFamilyFromBlastOfSPConjModule)),
                            currICEsIMEsStructure.internalIdentifier)
                    if commentITToAdd not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += commentITToAdd
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                else:
                    listUpstreamIntToAdd.remove(currSp)
                    if currSp.locusTag not in setLocusTagToNotConsiderAseManuallyCheckSent:
                        setIntegraseToManuallyCheckNoFamilyMatch.add(currSp)
                        listLocusTagInCommentNoFamilyMatch.append(currSp.locusTag)
    if listDownstreamIntToAdd:
        for currSp in reversed(listDownstreamIntToAdd):
            # slight change, for the best ?
            if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 0:
                if keepIntegraseAssignAsSureButAddCommentAboutFamily:
                    # if change text, change up and also in EMStructure.py around l 761
                    # commentITToAdd = "The integrase {} has never been previously seen associated to this conjugation module family. ".format(currSp.locusTag)
                    commentITToAdd = "The integrase {} has been previously seen associated with the conjugation module superfamily {} but not yet with the superfamily {} of structure {}. ".format(
                            currSp.locusTag,
                            ", ".join(str(i) for i in sorted(currSp.setSPICESuperFamilyFromBlast)),
                            ", ".join(str(i) for i in sorted(currICEsIMEsStructure.setICESuperFamilyFromBlastOfSPConjModule)),
                            currICEsIMEsStructure.internalIdentifier)
                    if commentITToAdd not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += commentITToAdd
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                else:
                    listDownstreamIntToAdd.remove(currSp)
                    if currSp.locusTag not in setLocusTagToNotConsiderAseManuallyCheckSent:
                        setIntegraseToManuallyCheckNoFamilyMatch.add(currSp)
                        listLocusTagInCommentNoFamilyMatch.append(currSp.locusTag)
    if len(listLocusTagInCommentNoFamilyMatch) != 0:
        reprSetIntegraseToManuallyCheckNoFamilyMatch = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseToManuallyCheckNoFamilyMatch)
        if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
            setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
            if reprSetIntegraseToManuallyCheckNoFamilyMatch in setIntegraseManuallyCheck2comment:
                raise RuntimeError("Error in markIntegraseNeverPreviouslySeenAssociatedToConjugationModuleFamilyAsToManuallyCheck: reprSetIntegraseToManuallyCheckNoFamilyMatch ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                        hit.ListSPs.GetListProtIdsFromSetSP(reprSetIntegraseToManuallyCheckNoFamilyMatch), currICEsIMEsStructure.internalIdentifier))
            else:
                commentITToAdd = "The integrase(s) {} could be added to the structure {} but it has never been previously seen associated to this conjugation module super family ({}), please manually check. ".format(
                        " and/or ".join(str(i) for i in sorted(listLocusTagInCommentNoFamilyMatch)), currICEsIMEsStructure.internalIdentifier, ", ".join(str(i) for i in sorted(currICEsIMEsStructure.setICESuperFamilyFromBlastOfSPConjModule)))
                setIntegraseManuallyCheck2comment[reprSetIntegraseToManuallyCheckNoFamilyMatch] = commentITToAdd
                # for locusTagToAddCommentIT in listLocusTagInCommentNoFamilyMatch:
                #    icescreen_OO.addCommentToLocusTag2Comment(locusTagToAddCommentIT, commentITToAdd, locusTagIntegrase2Comment)
        else:
            setIntegraseManuallyCheck2comment = {}
            commentITToAdd = "The integrase(s) {} could be added to the structure {} but it has never been previously seen associated to this conjugation module super family ({}), please manually check. ".format(
                    " and/or ".join(str(i) for i in sorted(listLocusTagInCommentNoFamilyMatch)), currICEsIMEsStructure.internalIdentifier, ", ".join(str(i) for i in sorted(currICEsIMEsStructure.setICESuperFamilyFromBlastOfSPConjModule)))
            setIntegraseManuallyCheck2comment[reprSetIntegraseToManuallyCheckNoFamilyMatch] = commentITToAdd
            ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment
            # for locusTagToAddCommentIT in listLocusTagInCommentNoFamilyMatch:
            #        icescreen_OO.addCommentToLocusTag2Comment(locusTagToAddCommentIT, commentITToAdd, locusTagIntegrase2Comment)

        # currICEsIMEsStructure.comment += "The integrase(s) {} could be added to this ICE IME structure but it is unsure if the families would match, please manually check. "\
        #        .format(" and/or ".join(str(i) for i in sorted(listLocusTagInCommentNoFamilyMatch)))


# return True if this downstream integrase is a better fit for next struct
def checkThisDownstreamIntegraseIsABeABetterFitForNextStruct(currSp, currICEsIMEsStructure, nextICEsIMEsStructure):

    currDownstreamIntegraseIsABeABetterFitForNextStruct = False

    if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, nextICEsIMEsStructure, True) == 2:

        # famillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure)
        # famillyIntegraseSameAsConjModuleNextICEsIMEsStructure = isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, nextICEsIMEsStructure)
        # if famillyIntegraseSameAsConjModuleNextICEsIMEsStructure == 2 and famillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
        #    currDownstreamIntegraseIsABeABetterFitForNextStruct = True

        ICEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure)
        ICEFamillyIntegraseSameAsConjModuleNextICEsIMEsStructure = isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, nextICEsIMEsStructure)
        if ICEFamillyIntegraseSameAsConjModuleNextICEsIMEsStructure == 2 and ICEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
            currDownstreamIntegraseIsABeABetterFitForNextStruct = True

        IMEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure)
        IMEFamillyIntegraseSameAsConjModuleNextICEsIMEsStructure = isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, nextICEsIMEsStructure)
        if IMEFamillyIntegraseSameAsConjModuleNextICEsIMEsStructure == 2 and IMEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
            currDownstreamIntegraseIsABeABetterFitForNextStruct = True

    if DEBUG:
        if currDownstreamIntegraseIsABeABetterFitForNextStruct:
            print("TRUE checkThisDownstreamIntegraseIsABeABetterFitForNextStruct for currSp = {}, currICEsIMEsStructure = {}, nextICEsIMEsStructure = {}".format(
                    currSp.locusTag, currICEsIMEsStructure.internalIdentifier, nextICEsIMEsStructure.internalIdentifier))

    return currDownstreamIntegraseIsABeABetterFitForNextStruct


# return True if this upstream integrase is a better fit for previous struct
def checkThisUpstreamIntegraseIsABeABetterFitForPreviousStruct(currSp, currICEsIMEsStructure, previousICEsIMEsStructure):

    currUpstreamIntegraseIsABeABetterFitForPreviousStruct = False
    if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, previousICEsIMEsStructure, False) == 2:

        # famillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure)
        # famillyIntegraseSameAsConjModulePreviousICEsIMEsStructure = isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, previousICEsIMEsStructure)
        # if famillyIntegraseSameAsConjModulePreviousICEsIMEsStructure == 2 and famillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
        #    currUpstreamIntegraseIsABeABetterFitForPreviousStruct = True

        ICEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure)
        ICEFamillyIntegraseSameAsConjModulePreviousICEsIMEsStructure = isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, previousICEsIMEsStructure)
        if ICEFamillyIntegraseSameAsConjModulePreviousICEsIMEsStructure == 2 and ICEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
            currUpstreamIntegraseIsABeABetterFitForPreviousStruct = True

        IMEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure = isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure)
        IMEFamillyIntegraseSameAsConjModulePreviousICEsIMEsStructure = isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, previousICEsIMEsStructure)
        if IMEFamillyIntegraseSameAsConjModulePreviousICEsIMEsStructure == 2 and IMEFamillyIntegraseSameAsConjModuleCurrICEsIMEsStructure < 2:
            currUpstreamIntegraseIsABeABetterFitForPreviousStruct = True

    if DEBUG:
        if currUpstreamIntegraseIsABeABetterFitForPreviousStruct:
            print("TRUE checkThisUpstreamIntegraseIsABeABetterFitForPreviousStruct for currSp = {}, currICEsIMEsStructure = {}, previousICEsIMEsStructure = {}".format(
                    currSp.locusTag, currICEsIMEsStructure.internalIdentifier, previousICEsIMEsStructure.internalIdentifier))
    return currUpstreamIntegraseIsABeABetterFitForPreviousStruct


# integrase can be associated to any conjugation module without restriction on families. This method build a dictionary for integrases that have already be seen associated to a conjugation module but now that the family for integrase has been removed from the database, this method is not usefull anymore.
def buildamFily2SetPutativeIntegrases(listSPs):
    family2SetPutativeIntegrases = {}
    for currSP in listSPs:
        if (currSP.SPType in icescreen_OO.setIntegraseNames):
            if currSP.SPDetectedByBlast == 1 and len(currSP.setSPICEFamilyFromBlast) != 0:
                for currFamilyFromBlast in currSP.setSPICEFamilyFromBlast:
                    if currFamilyFromBlast in family2SetPutativeIntegrases:  # key already there
                        currSetIntegrases = family2SetPutativeIntegrases[currFamilyFromBlast]
                        currSetIntegrases.add(currSP)
                    else:  # key not there
                        currSetIntegrases = set()
                        currSetIntegrases.add(currSP)
                        family2SetPutativeIntegrases[currFamilyFromBlast] = currSetIntegrases
            if currSP.SPDetectedByBlast == 1 and len(currSP.setSPIMEFamilyFromBlast) != 0:
                for currFamilyFromBlast in currSP.setSPIMEFamilyFromBlast:
                    if currFamilyFromBlast in family2SetPutativeIntegrases:  # key already there
                        currSetIntegrases = family2SetPutativeIntegrases[currFamilyFromBlast]
                        currSetIntegrases.add(currSP)
                    else:  # key not there
                        currSetIntegrases = set()
                        currSetIntegrases.add(currSP)
                        family2SetPutativeIntegrases[currFamilyFromBlast] = currSetIntegrases
        elif (currSP.SPType == "Relaxase" or currSP.SPType == "Coupling protein" or currSP.SPType == "VirB4"):
            pass
        else:
            raise RuntimeError("Error in buildamFily2SetPutativeIntegrases: unrecognized currSP.SPType = {}".format(currSP.SPType))
    return family2SetPutativeIntegrases


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
        for currFamilyFromBlast in currICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule:
            if currFamilyFromBlast in family2SetICEsIMEsStructures:  # key already there
                currSetICEsIMEsStructures = family2SetICEsIMEsStructures[currFamilyFromBlast]
                currSetICEsIMEsStructures.add(currICEsIMEsStructure)
            else:  # key not there
                currSetICEsIMEsStructures = set()
                currSetICEsIMEsStructures.add(currICEsIMEsStructure)
                family2SetICEsIMEsStructures[currFamilyFromBlast] = currSetICEsIMEsStructures

    return family2SetICEsIMEsStructures


# return:
# - 0 if integrase not correctly oriented and this fact is independant of SP in conflict
# - 1 the set of SP in conflict that trigger the fact that the integrase is not correctly oriented
# - 2 if integrase correctly oriented and this fact is independant of SP in conflict
def isIntegraseCorrectlyOrientedForICEsIMEsStructure(
        integraseSent,
        ICEsIMEsStructureSent,
        integraseIsUpstreamOfICEsIMEsStructure):
    # Upstream ICE integrase needs to be on the - strand
    # Downstream ICE integrase needs to be on the + strand
    valueToReturn = 2

    setVirB4InConflictInICEsIMEsStructureSent = set()
    setVirB4NotInConflictInICEsIMEsStructureSent = set()
    for currVirB4 in ICEsIMEsStructureSent.listVirB4:
        if currVirB4 in ICEsIMEsStructureSent.setSPInConflict:
            setVirB4InConflictInICEsIMEsStructureSent.add(currVirB4)
        else:
            setVirB4NotInConflictInICEsIMEsStructureSent.add(currVirB4)

    # at least 1 virB' not in conflict
    if len(setVirB4NotInConflictInICEsIMEsStructureSent) >= 1:
        if integraseIsUpstreamOfICEsIMEsStructure and integraseSent.strand == "+":
            valueToReturn = 0
        if not integraseIsUpstreamOfICEsIMEsStructure and integraseSent.strand == "-":
            valueToReturn = 0
    elif len(setVirB4InConflictInICEsIMEsStructureSent) >= 1:
        if integraseIsUpstreamOfICEsIMEsStructure and integraseSent.strand == "+":
            valueToReturn = setVirB4InConflictInICEsIMEsStructureSent
        if not integraseIsUpstreamOfICEsIMEsStructure and integraseSent.strand == "-":
            valueToReturn = setVirB4InConflictInICEsIMEsStructureSent

    if DEBUG:
        print("{} isIntegraseCorrectlyOrientedForICEsIMEsStructure: integraseSent = {} (len(setVirB4InConflictInICEsIMEsStructureSent) = {}) ; ICEsIMEsStructureSent = {} ; integraseIsUpstreamOfICEsIMEsStructure = {} ".format(
                repr(valueToReturn), integraseSent.locusTag, len(setVirB4InConflictInICEsIMEsStructureSent), ICEsIMEsStructureSent.internalIdentifier, integraseIsUpstreamOfICEsIMEsStructure))

    # if valueToReturn == -1:
    #    raise RuntimeError("Error in isIntegraseCorrectlyOrientedForICEsIMEsStructure: valueToReturn == -1 ; integraseSent = {} (len(setVirB4InConflictInICEsIMEsStructureSent) = {}) ; ICEsIMEsStructureSent = {} ; integraseIsUpstreamOfICEsIMEsStructure = {} " \
    #          .format(integraseSent.locusTag, len(setVirB4InConflictInICEsIMEsStructureSent), ICEsIMEsStructureSent.internalIdentifier, integraseIsUpstreamOfICEsIMEsStructure))

    return valueToReturn


# return 0 if different, 1 if one or either do not have family info, and 2 if both are same family
def isIMEFamillyIntegraseSameAsIMEFamilyConjModule(integraseSent, ICEsIMEsStructureSent):
    IMEFamillyIntegraseSameAsIMEFamilyConjModule = -1
    if integraseSent.SPDetectedByBlast == 1 and len(integraseSent.setSPIMEFamilyFromBlast) != 0:
        if len(ICEsIMEsStructureSent.setIMEFamilyFromBlastOfSPConjModule) != 0:
            setInteresctSPFamilies = integraseSent.setSPIMEFamilyFromBlast.intersection(ICEsIMEsStructureSent.setIMEFamilyFromBlastOfSPConjModule)
            if len(setInteresctSPFamilies) != 0:
                IMEFamillyIntegraseSameAsIMEFamilyConjModule = 2
            else:
                IMEFamillyIntegraseSameAsIMEFamilyConjModule = 0
        else:
            IMEFamillyIntegraseSameAsIMEFamilyConjModule = 1
    else:
        IMEFamillyIntegraseSameAsIMEFamilyConjModule = 1
    if IMEFamillyIntegraseSameAsIMEFamilyConjModule == -1:
        raise RuntimeError("Error in isIMEFamillyIntegraseSameAsIMEFamilyConjModule: IMEFamillyIntegraseSameAsIMEFamilyConjModule == -1 for integraseSent.SPType = {} and ICEsIMEsStructureSent = {}".format(
                integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    # if DEBUG:
    #    print("isIMEFamillyIntegraseSameAsIMEFamilyConjModule = {} for integraseSent = {}, ICEsIMEsStructureSent = {}" \
    #              .format(str(IMEFamillyIntegraseSameAsIMEFamilyConjModule), integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    return IMEFamillyIntegraseSameAsIMEFamilyConjModule


# return 0 if different, 1 if one or either do not have family info, and 2 if both are same family
def isICEFamillyIntegraseSameAsICEFamilyConjModule(integraseSent, ICEsIMEsStructureSent):
    ICEFamillyIntegraseSameAsICEFamilyConjModule = -1
    if integraseSent.SPDetectedByBlast == 1 and len(integraseSent.setSPICEFamilyFromBlast) != 0:
        if len(ICEsIMEsStructureSent.setICEFamilyFromBlastOfSPConjModule) != 0:
            setInteresctSPFamilies = integraseSent.setSPICEFamilyFromBlast.intersection(ICEsIMEsStructureSent.setICEFamilyFromBlastOfSPConjModule)
            if len(setInteresctSPFamilies) != 0:
                ICEFamillyIntegraseSameAsICEFamilyConjModule = 2
            else:
                ICEFamillyIntegraseSameAsICEFamilyConjModule = 0
        else:
            ICEFamillyIntegraseSameAsICEFamilyConjModule = 1
    else:
        ICEFamillyIntegraseSameAsICEFamilyConjModule = 1
    if ICEFamillyIntegraseSameAsICEFamilyConjModule == -1:
        raise RuntimeError("Error in isICEFamillyIntegraseSameAsICEFamilyConjModule: ICEFamillyIntegraseSameAsICEFamilyConjModule == -1 for integraseSent.SPType = {} and ICEsIMEsStructureSent = {}".format(
                integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    # if DEBUG:
    #    print("isICEFamillyIntegraseSameAsICEFamilyConjModule = {} for integraseSent = {}, ICEsIMEsStructureSent = {}" \
    #              .format(str(ICEFamillyIntegraseSameAsICEFamilyConjModule), integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    return ICEFamillyIntegraseSameAsICEFamilyConjModule


# return 0 if different, 1 if one or either do not have superfamily info, and 2 if both are same superfamily
def isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(integraseSent, ICEsIMEsStructureSent):
    ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = -1
    if integraseSent.SPDetectedByBlast == 1 and len(integraseSent.setSPICESuperFamilyFromBlast) != 0:
        if len(ICEsIMEsStructureSent.setICESuperFamilyFromBlastOfSPConjModule) != 0:
            setInteresctSPSuperFamilies = integraseSent.setSPICESuperFamilyFromBlast.intersection(ICEsIMEsStructureSent.setICESuperFamilyFromBlastOfSPConjModule)
            if len(setInteresctSPSuperFamilies) != 0:
                ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 2
            else:
                ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 0
            # for currFamilyFromBlast in integraseSent.setSPICESuperFamilyFromBlast:
            #    if currFamilyFromBlast in ICEsIMEsStructureSent.setICESuperFamilyFromBlastOfSPConjModule:
            #        ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 2
            # if ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule == -1:
            #    ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 0
        else:
            ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 1
    else:
        ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule = 1
    if ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule == -1:
        raise RuntimeError("Error in isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule: ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule == -1 for integraseSent.SPType = {} and ICEsIMEsStructureSent = {}".format(
                integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    # if DEBUG:
    #    print("isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule = {} for integraseSent = {}, ICEsIMEsStructureSent = {}" \
    #              .format(str(ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule), integraseSent.locusTag, ICEsIMEsStructureSent.internalIdentifier))
    return ICESuperFamillyIntegraseSameAsICESuperFamilyConjModule


# look for upstream or downstream subsequent integrase that have been preivously seen associated with the conj module
def checkForObviousIntegraseUpstreamAndDownstreamToAdd(
        EMStructureToAdd,
        listSPs,
        idxCurrentSP,
        allowAdjacentIntegraseOnlyForSer,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):

    setIntegraseTypeToCheck = set()
    if allowAdjacentIntegraseOnlyForSer == "YES":
        setIntegraseTypeToCheck = icescreen_OO.setIntegraseSerNames
    elif allowAdjacentIntegraseOnlyForSer == "NO":
        setIntegraseTypeToCheck = icescreen_OO.setIntegraseNames
    else:
        raise RuntimeError("Error in checkForObviousIntegraseUpstreamAndDownstreamToAdd: unrecognized allowAdjacentIntegraseOnlyForSer = {}".format(allowAdjacentIntegraseOnlyForSer))

    listUpstreamIntToAdd = []
    idxSPsUpstreamToCheck = idxCurrentSP - len(EMStructureToAdd.listOrderedSPs)
    if idxSPsUpstreamToCheck >= 0:
        listSPsUpstream = listSPs[:(idxSPsUpstreamToCheck + 1)]
        for currSp in reversed(listSPsUpstream):
            # print("HERE up {}".format(currSp.locusTag))
            if listUpstreamIntToAdd:
                if (abs(listUpstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                    if currSp.SPType in setIntegraseTypeToCheck:
                        if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, EMStructureToAdd, True) != 2:
                            # if EMStructureToAdd.listVirB4 and currSp.strand == "+":
                            break  # int for ICE need to be on the right strand
                        else:
                            if currSp.SPType == listUpstreamIntToAdd[-1].SPType and currSp.strand == listUpstreamIntToAdd[-1].strand:
                                listUpstreamIntToAdd.append(currSp)
                                continue
                            else:
                                # EMStructureToAdd.comment += "The upstream integrases {} and {} are adjacent but not of the same type, {} has been discarded.".format(listUpstreamIntToAdd[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                break
                    else:
                        break
                else:  # SP are not adjacent in genome
                    break
            else:  # listIntegraseUpstream is empty
                if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                    if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, EMStructureToAdd, True) != 2:
                        # if EMStructureToAdd.listVirB4 and currSp.strand == "+":
                        break  # int for ICE need to be on the right strand
                    else:
                        # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #    listUpstreamIntToAdd.append(currSp)
                        #    continue
                        if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                            listUpstreamIntToAdd.append(currSp)
                            continue
                        elif isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                            listUpstreamIntToAdd.append(currSp)
                            continue
                        else:
                            break
                else:
                    break

    listDownstreamIntToAdd = []
    idxSPsDownstreamToCheck = idxCurrentSP + 1
    if idxSPsDownstreamToCheck < len(listSPs):
        listSPsDownstream = listSPs[idxSPsDownstreamToCheck:]
        for currSp in listSPsDownstream:
            # print("HERE down {}".format(currSp.locusTag))
            if listDownstreamIntToAdd:  # listIntegraseDownstream is not empty
                if (abs(listDownstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                    if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                        # if EMStructureToAdd.listVirB4 and currSp.strand == "-":
                        if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, EMStructureToAdd, False) != 2:
                            # EMStructureToAdd.comment += "The downstream integrase {} is on the - strand, it has been discarded. ".format(currSp.locusTag)
                            break
                        else:
                            if currSp.SPType == listDownstreamIntToAdd[-1].SPType and currSp.strand == listDownstreamIntToAdd[-1].strand:
                                listDownstreamIntToAdd.append(currSp)
                                continue
                            else:
                                # EMStructureToAdd.comment += "The downstream integrases {} and {} are adjacent but not of the same type, {} has been discarded. ".format(listDownstreamIntToAdd[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                break
                    else:
                        break
                else:  # SP are not adjacent in genome
                    break
            else:  # listIntegraseDownstream is empty
                if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                    # if EMStructureToAdd.listVirB4 and currSp.strand == "-":
                    if isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, EMStructureToAdd, False) != 2:
                        # EMStructureToAdd.comment += "The downstream integrase {} is on the - strand, it has been discarded. ".format(currSp.locusTag)
                        break
                    else:
                        # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #    listDownstreamIntToAdd.append(currSp)
                        #    continue
                        if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                            listDownstreamIntToAdd.append(currSp)
                            continue
                        elif isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                            listDownstreamIntToAdd.append(currSp)
                            continue
                        else:
                            break
                else:
                    break

    # if both int up and downstream, use info on family to disciminate if possible
    if listUpstreamIntToAdd and listDownstreamIntToAdd:
        # can not decide, integrase situation will be dealt with latter
        if EMStructureToAdd in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
            raise RuntimeError("Error in checkForObviousIntegraseUpstreamAndDownstreamToAdd: EMStructureToAdd {} in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure".format(EMStructureToAdd.internalIdentifier))
        else:
            listOfListUpstreamAndDownstreamIntToAdd = []
            listOfListUpstreamAndDownstreamIntToAdd.append(listUpstreamIntToAdd)
            listOfListUpstreamAndDownstreamIntToAdd.append(listDownstreamIntToAdd)
            structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[EMStructureToAdd] = listOfListUpstreamAndDownstreamIntToAdd

    elif listUpstreamIntToAdd:
        for currSp in listUpstreamIntToAdd:
            # if next structure shares also this obvious integrase, remove it
            if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                # raise RuntimeError("Error in checkForObviousIntegraseUpstreamAndDownstreamToAdd: integrase {} can not be rattached to two structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure structures {} and {}." \
                #           .format(currSp.locusTag, structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[currSp].internalIdentifier, EMStructureToAdd.internalIdentifier))
                EMStructurePreviouslyAdded = dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd[currSp]

                strCommentIT = "Integrase {} could be associated to both structures {} and {} because they are adjacent and this integrase has already be seen rattached to their conjugaison module family. ".format(
                        currSp.locusTag, EMStructureToAdd.internalIdentifier, EMStructurePreviouslyAdded.internalIdentifier)
                if strCommentIT not in EMStructureToAdd.comment:
                    EMStructureToAdd.comment += strCommentIT
                if strCommentIT not in EMStructurePreviouslyAdded.comment:
                    EMStructurePreviouslyAdded.comment += strCommentIT
                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                # take off key and clean up integrase and comment in already registred structure
                strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                        currSp.locusTag, EMStructurePreviouslyAdded.internalIdentifier)

                # print("{}".format(strCommentIT))
                # print("{}".format(strCommentToRemove))
                # print("EMStructureToAdd.listIntegraseUpstream: {}".format(hit.ListSPs.GetListProtIdsFromListSP(EMStructureToAdd.listIntegraseUpstream)))
                # print("EMStructureToAdd.listIntegraseDownstream: {}".format(hit.ListSPs.GetListProtIdsFromListSP(EMStructureToAdd.listIntegraseDownstream)))
                # print("EMStructurePreviouslyAdded.listIntegraseUpstream: {}".format(hit.ListSPs.GetListProtIdsFromListSP(EMStructurePreviouslyAdded.listIntegraseUpstream)))
                # print("EMStructurePreviouslyAdded.listIntegraseDownstream: {}".format(hit.ListSPs.GetListProtIdsFromListSP(EMStructurePreviouslyAdded.listIntegraseDownstream)))

                if strCommentToRemove in EMStructurePreviouslyAdded.comment:
                    EMStructurePreviouslyAdded.comment = EMStructurePreviouslyAdded.comment.replace(strCommentToRemove, "")
                icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)

                useIntegraseCDSDistanceToChooseBetweenUpstreamAndDownstreamStructures(
                        currSp,
                        EMStructurePreviouslyAdded,
                        EMStructureToAdd,
                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                        True,
                        locusTagIntegrase2Comment)
#                EMStructurePreviouslyAdded.listIntegraseDownstream.remove(currSp)
#                EMStructurePreviouslyAdded.listOrderedSPs.remove(currSp)
#                EMStructurePreviouslyAdded.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)
#                EMStructureToAdd.setIntegraseLocusTagsToManuallyCheck.add(currSp.locusTag)

            else:
                # structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[currSp] = EMStructureToAdd
                dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd[currSp] = EMStructureToAdd
                strCommentIT = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                        currSp.locusTag, EMStructureToAdd.internalIdentifier)
                if strCommentIT not in EMStructureToAdd.comment:
                    EMStructureToAdd.comment += strCommentIT
                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                EMStructureToAdd.listIntegraseUpstream.append(currSp)
                EMStructureToAdd.listOrderedSPs.append(currSp)

    elif listDownstreamIntToAdd:
        for currSp in listDownstreamIntToAdd:
            # if next structure shares also this obvious integrase, remove it
            if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                raise RuntimeError("Error in checkForObviousIntegraseUpstreamAndDownstreamToAdd: integrase {} can not be rattached to two structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure structures {} and {}.".format(
                        currSp.locusTag, structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[currSp].internalIdentifier, EMStructureToAdd.internalIdentifier))

            else:
                # structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure[currSp] = EMStructureToAdd
                dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd[currSp] = EMStructureToAdd
                strCommentIT = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. ".format(
                        currSp.locusTag, EMStructureToAdd.internalIdentifier)
                if strCommentIT not in EMStructureToAdd.comment:
                    EMStructureToAdd.comment += strCommentIT
                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                EMStructureToAdd.listIntegraseDownstream.append(currSp)
                EMStructureToAdd.listOrderedSPs.append(currSp)

    if listUpstreamIntToAdd or listDownstreamIntToAdd:
        EMStructureToAdd.refreshListIdxOrderedSPs()


# add both upstream and downstream valid integrase to a conj module
def addSPIntegraseUpstreamAndDownstream(
        listSPs,
        listICEsIMEsStructures,
        maxNumberCDSForFilterIMESize,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        setLocusTagToNotConsiderAseManuallyCheckSent,
        countIterAddSPIntegraseUpstreamAndDownstream,
        allowAdjacentIntegraseOnlyForSer,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):

    doNotTakeIntoAccountIfManuallyCheckButAnIntegraseMoreSureIsRegistred = True

    countIterAddSPIntegraseUpstreamAndDownstream += 1

    family2SetICEsIMEsStructures = buildFamily2SetICEsIMEsStructures(listICEsIMEsStructures)
    family2SetPutativeIntegrases = buildamFily2SetPutativeIntegrases(listSPs)

    integraseAttributedSureUpstream2setICEsIMEsStructure = {}
    integraseAttributedSureDownstream2setICEsIMEsStructure = {}
    ICEsIMEsStructure2setIntegraseMoreSure = {}
    ICEsIMEsStructure2setIntegraseManuallyCheck2comment = {}

    for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):

        if DEBUG:
            print("\n\n ** Starting currICEsIMEsStructure = {} ;  idxCurrICEsIMEsStructure = {} ; setLocusTagToNotConsiderAseManuallyCheckSent = {} ; countIterAddSPIntegraseUpstreamAndDownstream = {} ".format(
                    currICEsIMEsStructure.internalIdentifier, str(idxCurrICEsIMEsStructure), repr(setLocusTagToNotConsiderAseManuallyCheckSent), str(countIterAddSPIntegraseUpstreamAndDownstream)))

        # integrase have already been added with checkForObviousIntegraseUpstreamAndDownstreamToAdd
        if len(currICEsIMEsStructure.listIntegraseUpstream) >= 1 or len(currICEsIMEsStructure.listIntegraseDownstream) >= 1:
            for currIntegraseUpstream in currICEsIMEsStructure.listIntegraseUpstream:
                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseMoreSure:
                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure]
                    setIntegraseMoreSure.add(currIntegraseUpstream)
                else:
                    setIntegraseMoreSure = set()
                    setIntegraseMoreSure.add(currIntegraseUpstream)
                    ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure] = setIntegraseMoreSure
            for currIntegraseDownstream in currICEsIMEsStructure.listIntegraseDownstream:
                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseMoreSure:
                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure]
                    setIntegraseMoreSure.add(currIntegraseDownstream)
                else:
                    setIntegraseMoreSure = set()
                    setIntegraseMoreSure.add(currIntegraseDownstream)
                    ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure] = setIntegraseMoreSure
            continue  # integrase have already been added with checkForObviousIntegraseUpstreamAndDownstreamToAdd

        # flagDoNotCheckIntegraseStrand = False
        # if len(currICEsIMEsStructure.listVirB4) >= 1:
        #    flagDoNotCheckIntegraseStrand = True
        setIntegraseTypeToCheck = set()
        if allowAdjacentIntegraseOnlyForSer == "YES":
            setIntegraseTypeToCheck = icescreen_OO.setIntegraseSerNames
        elif allowAdjacentIntegraseOnlyForSer == "NO":
            setIntegraseTypeToCheck = icescreen_OO.setIntegraseNames
        else:
            raise RuntimeError("Error in addSPIntegraseUpstreamAndDownstream: unrecognized allowAdjacentIntegraseOnlyForSer = {}".format(allowAdjacentIntegraseOnlyForSer))

        setIntegraseDiscarded = set()
        if currICEsIMEsStructure.delMerging_idxListUpstreamStructure == -1:  # skip if deleted after merge

            listUpstreamIntToAdd = []
            foundUpstreamIntWithMatchFamily = False
            idxSPsUpstreamToCheck = currICEsIMEsStructure.listOrderedSPs[0].idxInListSP - 1
            if idxSPsUpstreamToCheck >= 0:
                listSPsUpstream = listSPs[:(idxSPsUpstreamToCheck + 1)]
                for currSp in reversed(listSPsUpstream):

                    if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                        # strCommentIT = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. " \
                        #    .format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                        # if strCommentIT not in currICEsIMEsStructure.comment:
                        #    currICEsIMEsStructure.comment += strCommentIT
                        # icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                        continue

                    # print("HERE up {}".format(currSp.locusTag))
                    if listUpstreamIntToAdd:
                        if (abs(listUpstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                                valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, currICEsIMEsStructure, True)
                                if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                    setIntegraseDiscarded.add(currSp.locusTag)
                                    strCommentIT = "Integrase {} can not be associated with the downstream structure {}, it needs to be on the - strand. ".format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    break
                                elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                    strCommentIT = "If {} is part of the structure {}, the upstream integrase {} would not be associated with this structure as it needs to be on the - strand. ".format(
                                            hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                        icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                if currSp.SPType == listUpstreamIntToAdd[-1].SPType and currSp.strand == listUpstreamIntToAdd[-1].strand:

                                    if idxCurrICEsIMEsStructure != 0:  # not the first element
                                        if checkThisUpstreamIntegraseIsABeABetterFitForPreviousStruct(currSp, currICEsIMEsStructure, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1]):
                                            if checkIfStructureHasBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structureAndIsNotInSetLocusTagToNotConsiderAsManuallyCheck(
                                                    listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1],
                                                    setLocusTagToNotConsiderAseManuallyCheckSent,
                                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                                    locusTagIntegrase2Comment):
                                                # if listICEsIMEsStructures[idxCurrICEsIMEsStructure-1] in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
                                                strCommentIT = "Integrase {} has a family that was already seen associated with the upstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                        currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            else:
                                                strCommentIT = "Integrase {} is not associated with the downstream structure {} because its family has already be seen associated with the upstream structure {}. ".format(
                                                        currSp.locusTag, currICEsIMEsStructure.internalIdentifier, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                                strCommentIT = "Integrase {} has a family that was already seen associated with the upstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                        currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                                if strCommentIT in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentIT, "")
                                                icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                                break

                                    listUpstreamIntToAdd.append(currSp)
                                    # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                    if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                        foundUpstreamIntWithMatchFamily = True
                                    if isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                        foundUpstreamIntWithMatchFamily = True
                                    continue
                                else:
                                    setIntegraseDiscarded.add(currSp.locusTag)
                                    strCommentIT = "The upstream integrases {} and {} are adjacent but not of the same type or strand, {} is not considered as part of a tandem.".format(listUpstreamIntToAdd[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    break
                            else:
                                break
                        else:  # SP are not adjacent in genome
                            break
                    else:  # listUpstreamIntToAdd is empty
                        if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"
                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, currICEsIMEsStructure, True)
                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "Integrase {} can not be associated with the downstream structure {}, it needs to be on the - strand. ".format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                                if strCommentIT not in currICEsIMEsStructure.comment:
                                    currICEsIMEsStructure.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                strCommentIT = "If {} is part of the structure {}, the upstream integrase {} would not be associated with this structure as it needs to be on the - strand. ".format(
                                        hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                if strCommentIT not in currICEsIMEsStructure.comment:
                                    currICEsIMEsStructure.comment += strCommentIT
                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                            if idxCurrICEsIMEsStructure != 0:  # not the first element
                                if checkThisUpstreamIntegraseIsABeABetterFitForPreviousStruct(currSp, currICEsIMEsStructure, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1]):
                                    if checkIfStructureHasBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structureAndIsNotInSetLocusTagToNotConsiderAsManuallyCheck(
                                            listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1],
                                            setLocusTagToNotConsiderAseManuallyCheckSent,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                            locusTagIntegrase2Comment):
                                        # if listICEsIMEsStructures[idxCurrICEsIMEsStructure-1] in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
                                        strCommentIT = "Integrase {} has a family that was already seen associated with the upstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                        if strCommentIT not in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment += strCommentIT
                                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    else:
                                        strCommentIT = "Integrase {} is not associated with the downstream structure {} because its family has already be seen associated with the upstream structure {}. ".format(
                                                currSp.locusTag, currICEsIMEsStructure.internalIdentifier, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                        if strCommentIT not in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment += strCommentIT
                                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                        strCommentIT = "Integrase {} has a family that was already seen associated with the upstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure - 1].internalIdentifier)
                                        if strCommentIT in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentIT, "")
                                        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                        break
                            listUpstreamIntToAdd.append(currSp)
                            # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                            if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                foundUpstreamIntWithMatchFamily = True
                            if isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                foundUpstreamIntWithMatchFamily = True
                            continue
                        else:
                            break

            if DEBUG:
                print(" - 1. listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))

            listDownstreamIntToAdd = []
            foundDownstreamIntWithMatchFamily = False
            idxSPsDownstreamToCheck = currICEsIMEsStructure.listOrderedSPs[-1].idxInListSP + 1
            if idxSPsDownstreamToCheck < len(listSPs):
                listSPsDownstream = listSPs[idxSPsDownstreamToCheck:]
                for currSp in listSPsDownstream:

                    if currSp in dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd:
                        # strCommentIT = "Integrase {} has been associated to the structure {} because they are adjacent, this integrase has already be seen rattached to this conjugaison module family, and there is no upstream/downstream ambiguity. " \
                        #    .format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                        # if strCommentIT not in currICEsIMEsStructure.comment:
                        #    currICEsIMEsStructure.comment += strCommentIT
                        # icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                        continue

                    #  print("HERE down {}".format(currSp.locusTag))
                    if listDownstreamIntToAdd:  # listIntegraseDownstream is not empty
                        if (abs(listDownstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2):  # Integrases are adjacent in genome or separated by a CDS
                            if currSp.SPType in setIntegraseTypeToCheck:  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"

                                # if not isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, currICEsIMEsStructure, False):
                                #    setIntegraseDiscarded.add(currSp.locusTag)
                                #    strCommentIT = "Integrase {} can not be associated with the upstream ICE {}, it needs to be on the + strand. " \
                                #            .format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                                #    if strCommentIT not in currICEsIMEsStructure.comment:
                                #        currICEsIMEsStructure.comment += strCommentIT
                                #    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                #    break

                                valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, currICEsIMEsStructure, False)
                                if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                    setIntegraseDiscarded.add(currSp.locusTag)
                                    strCommentIT = "Integrase {} can not be associated with the upstream structure {}, it needs to be on the + strand. ".format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    break
                                elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                    strCommentIT = "If {} is part of the structure {}, the downstream integrase {} would not be associated with this structure as it needs to be on the + strand. ".format(
                                            hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                        icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                if currSp.SPType == listDownstreamIntToAdd[-1].SPType and currSp.strand == listDownstreamIntToAdd[-1].strand:

                                    if idxCurrICEsIMEsStructure != len(listICEsIMEsStructures) - 1:  # not the last element

                                        if checkThisDownstreamIntegraseIsABeABetterFitForNextStruct(currSp, currICEsIMEsStructure, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1]):
                                            if checkIfStructureHasBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structureAndIsNotInSetLocusTagToNotConsiderAsManuallyCheck(
                                                    listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1],
                                                    setLocusTagToNotConsiderAseManuallyCheckSent,
                                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                                    locusTagIntegrase2Comment):
                                                # if listICEsIMEsStructures[idxCurrICEsIMEsStructure+1] in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
                                                strCommentIT = "Integrase {} has a family that was already seen associated with the downstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                        currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            else:
                                                strCommentIT = "Integrase {} is not associated with the upstream structure {} because its family has already be seen associated with the downstream structure {}. ".format(
                                                        currSp.locusTag, currICEsIMEsStructure.internalIdentifier, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                                strCommentIT = "Integrase {} has a family that was already seen associated with the downstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                        currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                                if strCommentIT in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentIT, "")
                                                icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                                break

                                    listDownstreamIntToAdd.append(currSp)
                                    # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                    if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                        foundDownstreamIntWithMatchFamily = True
                                    if isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                        foundDownstreamIntWithMatchFamily = True
                                    continue
                                else:
                                    setIntegraseDiscarded.add(currSp.locusTag)
                                    strCommentIT = "The downstream integrases {} and {} are adjacent but not of the same type or strand, {} is not considered as part of a tandem. ".format(listDownstreamIntToAdd[-1].locusTag, currSp.locusTag, currSp.locusTag)
                                    if strCommentIT not in currICEsIMEsStructure.comment:
                                        currICEsIMEsStructure.comment += strCommentIT
                                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    break
                            else:
                                break
                        else:  # SP are not adjacent in genome
                            break
                    else:  # listDownstreamIntToAdd is empty
                        if (currSp.SPType in icescreen_OO.setIntegraseNames):  # currSp.SPType == "IntTyr" or currSp.SPType == "IntSer" or currSp.SPType == "DDE"

                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currSp, currICEsIMEsStructure, False)
                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                setIntegraseDiscarded.add(currSp.locusTag)
                                strCommentIT = "Integrase {} can not be associated with the upstream structure {}, it needs to be on the + strand. ".format(currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                                if strCommentIT not in currICEsIMEsStructure.comment:
                                    currICEsIMEsStructure.comment += strCommentIT
                                icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                break
                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                strCommentIT = "If {} is part of the structure {}, the downstream integrase {} would not be associated with this structure as it needs to be on the + strand. ".format(
                                        hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                if strCommentIT not in currICEsIMEsStructure.comment:
                                    currICEsIMEsStructure.comment += strCommentIT
                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                            if idxCurrICEsIMEsStructure != len(listICEsIMEsStructures) - 1:  # not the last element
                                if checkThisDownstreamIntegraseIsABeABetterFitForNextStruct(currSp, currICEsIMEsStructure, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1]):
                                    if checkIfStructureHasBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structureAndIsNotInSetLocusTagToNotConsiderAsManuallyCheck(
                                            listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1],
                                            setLocusTagToNotConsiderAseManuallyCheckSent,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                                            useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                                            locusTagIntegrase2Comment):
                                        # if listICEsIMEsStructures[idxCurrICEsIMEsStructure+1] in structureWithBothUpstreamAndDownstreamIntegraseSameFamilyAsConjModule2structure:
                                        strCommentIT = "Integrase {} has a family that was already seen associated with the downstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                        if strCommentIT not in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment += strCommentIT
                                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    else:
                                        strCommentIT = "Integrase {} is not associated with the upstream structure {} because its family has already be seen associated with the downstream structure {}. ".format(
                                                currSp.locusTag, currICEsIMEsStructure.internalIdentifier, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                        if strCommentIT not in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment += strCommentIT
                                        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                        strCommentIT = "Integrase {} has a family that was already seen associated with the downstream structure {} but this latter has both an upstream and downstream integrase already seen associated to this family. ".format(
                                                currSp.locusTag, listICEsIMEsStructures[idxCurrICEsIMEsStructure + 1].internalIdentifier)
                                        if strCommentIT in currICEsIMEsStructure.comment:
                                            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentIT, "")
                                        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                        break

                            listDownstreamIntToAdd.append(currSp)
                            # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                            if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                foundDownstreamIntWithMatchFamily = True
                            if isIMEFamillyIntegraseSameAsIMEFamilyConjModule(currSp, currICEsIMEsStructure) == 2:
                                foundDownstreamIntWithMatchFamily = True
                            continue
                        else:
                            break

            if DEBUG:
                print(" - 1. listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            # if both int up and downstream, use info on family to disciminate if possible
            if listUpstreamIntToAdd and listDownstreamIntToAdd and (len(currICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule) != 0 or len(currICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule) != 0):
                if foundUpstreamIntWithMatchFamily and foundDownstreamIntWithMatchFamily:
                    pass
                elif foundUpstreamIntWithMatchFamily:
                    listDownstreamIntToAdd.clear()
                elif foundDownstreamIntWithMatchFamily:
                    listUpstreamIntToAdd.clear()
                else:
                    pass

            if DEBUG:
                print(" - 1.5. listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
            if DEBUG:
                print(" - 1.5. listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            if foundUpstreamIntWithMatchFamily:
                pass  # can not find any better match
            elif foundDownstreamIntWithMatchFamily:
                pass  # can not find any better match
            else:
                # use family2SetICEsIMEsStructures and family2SetPutativeIntegrases to complement info on integrase
                tmpSetBothICEIMEFamilyFromBlastOfSPConjModule = set()
                tmpSetBothICEIMEFamilyFromBlastOfSPConjModule.update(currICEsIMEsStructure.setICEFamilyFromBlastOfSPConjModule)
                tmpSetBothICEIMEFamilyFromBlastOfSPConjModule.update(currICEsIMEsStructure.setIMEFamilyFromBlastOfSPConjModule)
                for currFamilyFromBlast in tmpSetBothICEIMEFamilyFromBlastOfSPConjModule:  # currICEsIMEsStructure.setFamilyFromBlastOfSPConjModule

                    onlyThisICEIMEStructIsRegistredForThisFamily = False
                    if currFamilyFromBlast in family2SetICEsIMEsStructures:
                        if currICEsIMEsStructure not in family2SetICEsIMEsStructures[currFamilyFromBlast]:
                            raise RuntimeError("Error in addSPIntegraseUpstreamAndDownstream: currICEsIMEsStructure {} not in family2SetICEsIMEsStructures[currFamilyFromBlast] {}".format(
                                    currICEsIMEsStructure.internalIdentifier, EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(family2SetICEsIMEsStructures[currFamilyFromBlast])))
                        if len(family2SetICEsIMEsStructures[currFamilyFromBlast]) == 1:
                            onlyThisICEIMEStructIsRegistredForThisFamily = True

                    setIntegraseCloseAndWithSimilarFamilyFromBlast = set()
                    if currFamilyFromBlast in family2SetPutativeIntegrases:  # this family has at least a registred Integrase
                        currSetPutativeIntegrasesIT = family2SetPutativeIntegrases[currFamilyFromBlast]
                        currSetICEsIMEsStructuresIT = family2SetICEsIMEsStructures[currFamilyFromBlast]
                        for currPutativeIntegrasesIT in currSetPutativeIntegrasesIT:
                            # is this integrase closer to our ICEsIMEsStructure ?
                            distWithOurICEsIMEsStructure = -1
                            if currPutativeIntegrasesIT.CDSPositionInGenome < currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:
                                distWithOurICEsIMEsStructure = abs(currPutativeIntegrasesIT.CDSPositionInGenome - currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome)
                            else:
                                distWithOurICEsIMEsStructure = abs(currPutativeIntegrasesIT.CDSPositionInGenome - currICEsIMEsStructure.listOrderedSPs[-1].CDSPositionInGenome)
                            ourICEsIMEsStructureIsClosest = True
                            for currICEsIMEsStructuresIT in currSetICEsIMEsStructuresIT:
                                if currICEsIMEsStructuresIT == currICEsIMEsStructure:
                                    continue
                                distWithOurOtherICEsIMEsStructure = -1
                                if currPutativeIntegrasesIT.CDSPositionInGenome < currICEsIMEsStructuresIT.listOrderedSPs[0].CDSPositionInGenome:
                                    distWithOurOtherICEsIMEsStructure = abs(currPutativeIntegrasesIT.CDSPositionInGenome - currICEsIMEsStructuresIT.listOrderedSPs[0].CDSPositionInGenome)
                                else:
                                    distWithOurOtherICEsIMEsStructure = abs(currPutativeIntegrasesIT.CDSPositionInGenome - currICEsIMEsStructuresIT.listOrderedSPs[-1].CDSPositionInGenome)
                                if distWithOurOtherICEsIMEsStructure < distWithOurICEsIMEsStructure:
                                    ourICEsIMEsStructureIsClosest = False
                                    break
                            if ourICEsIMEsStructureIsClosest:
                                if currPutativeIntegrasesIT.locusTag not in setIntegraseDiscarded:  # do not add currPutativeIntegrasesIT.locusTag not in setLocusTagToNotConsiderAseManuallyCheckSent and because of integrase that will be merged
                                    if len(currICEsIMEsStructure.listVirB4) > 0:
                                        if currPutativeIntegrasesIT.CDSPositionInGenome < currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:  # integrase is upstream

                                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currPutativeIntegrasesIT, currICEsIMEsStructure, True)
                                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                                strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the downstream structure {}, it needs to be on the - strand. ".format(
                                                        currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                                strCommentIT = "If {} is part of the structure {}, the upstream integrase {} would not be associated with this structure as it needs to be on the - strand. ".format(
                                                        hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                            # if not isIntegraseCorrectlyOrientedForICEsIMEsStructure(currPutativeIntegrasesIT, currICEsIMEsStructure, True):
                                            #    strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the downstream ICE {}, it needs to be on the - strand. " \
                                            #        .format(currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                            #    if strCommentIT not in currICEsIMEsStructure.comment:
                                            #        currICEsIMEsStructure.comment += strCommentIT
                                            #    icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            else:
                                                setIntegraseCloseAndWithSimilarFamilyFromBlast.add(currPutativeIntegrasesIT)
                                        elif currPutativeIntegrasesIT.CDSPositionInGenome > currICEsIMEsStructure.listOrderedSPs[-1].CDSPositionInGenome:  # integrase is downstream
                                            valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure = isIntegraseCorrectlyOrientedForICEsIMEsStructure(currPutativeIntegrasesIT, currICEsIMEsStructure, False)
                                            if valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure == 0:
                                                strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the upstream structure {}, it needs to be on the + strand. ".format(
                                                        currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            elif type(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) is set and len(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure) > 0:
                                                strCommentIT = "If {} is part of the structure {}, the downstream integrase {} would not be associated with this structure as it needs to be on the + strand. ".format(
                                                        hit.ListSPs.GetListProtIdsFromSetSP(valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure), currICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                                                if strCommentIT not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += strCommentIT
                                                for conflictVirB4IT in valueIsIntegraseCorrectlyOrientedForICEsIMEsStructure:
                                                    icescreen_OO.addCommentToLocusTag2Comment(conflictVirB4IT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                                            # if not isIntegraseCorrectlyOrientedForICEsIMEsStructure(currPutativeIntegrasesIT, currICEsIMEsStructure, False):
                                            #    strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the upstream ICE {}, it needs to be on the + strand. " \
                                            #        .format(currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                            #    if strCommentIT not in currICEsIMEsStructure.comment:
                                            #        currICEsIMEsStructure.comment += strCommentIT
                                            #    icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                            else:
                                                setIntegraseCloseAndWithSimilarFamilyFromBlast.add(currPutativeIntegrasesIT)
                                        else:
                                            strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the element {}, it needs to be upstream or downstream. ".format(
                                                    currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                            if strCommentIT not in currICEsIMEsStructure.comment:
                                                currICEsIMEsStructure.comment += strCommentIT
                                            icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)
                                    else:
                                        if currPutativeIntegrasesIT.CDSPositionInGenome < currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:  # integrase is upstream
                                            setIntegraseCloseAndWithSimilarFamilyFromBlast.add(currPutativeIntegrasesIT)
                                        elif currPutativeIntegrasesIT.CDSPositionInGenome > currICEsIMEsStructure.listOrderedSPs[-1].CDSPositionInGenome:  # integrase is downstream
                                            setIntegraseCloseAndWithSimilarFamilyFromBlast.add(currPutativeIntegrasesIT)
                                        else:
                                            strCommentIT = "Integrase {} previously seen associated to conjugation module family {} can not be associated with the element {}, it needs to be upstream or downstream. ".format(
                                                    currPutativeIntegrasesIT.locusTag, currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier)
                                            if strCommentIT not in currICEsIMEsStructure.comment:
                                                currICEsIMEsStructure.comment += strCommentIT
                                            icescreen_OO.addCommentToLocusTag2Comment(currPutativeIntegrasesIT.locusTag, strCommentIT, locusTagIntegrase2Comment)

                    if len(setIntegraseCloseAndWithSimilarFamilyFromBlast) != 0:  # this family has at least a registred Integrase that is closest

                        # setPutativeIntegraseWithSimilarFamilyFromBlast = family2SetPutativeIntegrases[currFamilyFromBlast]
                        if not listUpstreamIntToAdd and not listDownstreamIntToAdd:  # no Upstream or Downstream integrase found
                            # print("HERE 0 {}: setIntegraseCloseAndWithSimilarFamilyFromBlast = {} ; onlyThisICEIMEStructIsRegistredForThisFamily {} "\
                            #      .format(currICEsIMEsStructure.internalIdentifier, repr(setIntegraseCloseAndWithSimilarFamilyFromBlast), onlyThisICEIMEStructIsRegistredForThisFamily))

                            # merge ?
                            boolMergeIntWasPerformed = False
                            if len(setIntegraseCloseAndWithSimilarFamilyFromBlast) == 1 and onlyThisICEIMEStructIsRegistredForThisFamily:

                                for currIntegraseCloseAndWithSimilarFamilyFromBlast in setIntegraseCloseAndWithSimilarFamilyFromBlast:

                                    numberMaxIcesIMEsStructInBetween = 2

                                    if currIntegraseCloseAndWithSimilarFamilyFromBlast.CDSPositionInGenome < currICEsIMEsStructure.listOrderedSPs[0].CDSPositionInGenome:
                                        # merge upstream ?
                                        if idxCurrICEsIMEsStructure > 0:
                                            listICEsIMEsStructuresUpstream = listICEsIMEsStructures[:idxCurrICEsIMEsStructure]
                                            listICEsIMEsStructuresGuests = []
                                            for currentICEsIMEsStructuresUpstream in reversed(listICEsIMEsStructuresUpstream):
                                                if currentICEsIMEsStructuresUpstream.listOrderedSPs[0].CDSPositionInGenome < currIntegraseCloseAndWithSimilarFamilyFromBlast.CDSPositionInGenome:
                                                    break
                                                listICEsIMEsStructuresGuests.append(currentICEsIMEsStructuresUpstream)
                                            if len(listICEsIMEsStructuresGuests) == 0:
                                                srtCommentToAdd = "The upstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) but is not being merged into the structure {} because the number of would-be inbetween guest structures is null. ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier)
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                            elif len(listICEsIMEsStructuresGuests) > numberMaxIcesIMEsStructInBetween:
                                                srtCommentToAdd = "The upstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) but is not being merged into the structure {} because the number of would-be inbetween guest structures is too high ({} greater than max parameter {}). ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier,
                                                        str(len(listICEsIMEsStructuresGuests)),
                                                        str(numberMaxIcesIMEsStructInBetween))
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                            else:
                                                boolMergeIntWasPerformed = True
                                                # if countIterAddSPIntegraseUpstreamAndDownstream == 1:
                                                srtCommentToAdd = "The upstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) and is being merged into the structure {}. ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier)
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                                currICEsIMEsStructure.idxListDownstrStructMerged.append("I<" + str(listICEsIMEsStructuresGuests[-1].idxInSeedList))
                                                listUpstreamIntToAdd.append(currIntegraseCloseAndWithSimilarFamilyFromBlast)
                                                # host and guest will be done in fillUpColocalizedOtherICEsIMEsStructures afterward
                                                # pass
                                    else:
                                        # merge downstream ?
                                        if idxCurrICEsIMEsStructure < len(listICEsIMEsStructures) - 1:
                                            listICEsIMEsStructuresDownstream = listICEsIMEsStructures[(idxCurrICEsIMEsStructure + 1):]
                                            listICEsIMEsStructuresGuests = []
                                            for currentICEsIMEsStructuresDownstream in listICEsIMEsStructuresDownstream:
                                                if currentICEsIMEsStructuresDownstream.listOrderedSPs[0].CDSPositionInGenome > currIntegraseCloseAndWithSimilarFamilyFromBlast.CDSPositionInGenome:
                                                    break
                                                listICEsIMEsStructuresGuests.append(currentICEsIMEsStructuresDownstream)

                                            if len(listICEsIMEsStructuresGuests) == 0:
                                                srtCommentToAdd = "The downstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) but is not being merged into the structure {} because the number of would-be inbetween guest structures is null. ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier)
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                            elif len(listICEsIMEsStructuresGuests) > numberMaxIcesIMEsStructInBetween:
                                                srtCommentToAdd = "The downstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) but is not being merged into the structure {} because the number of would-be inbetween guest structures is too high ({} greater than max parameter {}). ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier,
                                                        str(len(listICEsIMEsStructuresGuests)),
                                                        str(numberMaxIcesIMEsStructInBetween))
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                            else:
                                                # if len(listICEsIMEsStructuresGuests) > 0  and len(listICEsIMEsStructuresGuests) < numberMaxIcesIMEsStructInBetween + 1:
                                                boolMergeIntWasPerformed = True
                                                # if countIterAddSPIntegraseUpstreamAndDownstream == 1:
                                                srtCommentToAdd = "The downstream integrase {} has been seen associated with the conjugaison module family (ICE {} ; IME {}) and is being merged into the structure {}. ".format(
                                                        currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag,
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPICEFamilyFromBlast)),
                                                        ", ".join(str(i) for i in sorted(currIntegraseCloseAndWithSimilarFamilyFromBlast.setSPIMEFamilyFromBlast)),
                                                        currICEsIMEsStructure.internalIdentifier)
                                                if srtCommentToAdd not in currICEsIMEsStructure.comment:
                                                    currICEsIMEsStructure.comment += srtCommentToAdd
                                                icescreen_OO.addCommentToLocusTag2Comment(currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag, srtCommentToAdd, locusTagIntegrase2Comment)
                                                currICEsIMEsStructure.idxListDownstrStructMerged.append("I>" + str(listICEsIMEsStructuresGuests[-1].idxInSeedList + 1))  # plus 1 because we else the last strucutre will not be taken into account
                                                listDownstreamIntToAdd.append(currIntegraseCloseAndWithSimilarFamilyFromBlast)
                                                # host and guest will be done in fillUpColocalizedOtherICEsIMEsStructures afterward
                                                # pass
                                    break

                            if not boolMergeIntWasPerformed:

                                setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel = set()
                                for currIntegraseCloseAndWithSimilarFamilyFromBlast in setIntegraseCloseAndWithSimilarFamilyFromBlast:
                                    if currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag in setLocusTagToNotConsiderAseManuallyCheckSent:
                                        if DEBUG:
                                            print("setIntegraseCloseAndWithSimilarFamilyFromBlast no Upstream or Downstream integrase found not boolMergeIntWasPerformed: {}".format(
                                                    currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag))
                                        setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel.add(currIntegraseCloseAndWithSimilarFamilyFromBlast)
                                    if currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag in currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck:
                                        # already taken car in method checkForObviousIntegrase propably
                                        setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel.add(currIntegraseCloseAndWithSimilarFamilyFromBlast)

                                for currIntegraseCloseAndWithSimilarFamilyFromBlastToDel in setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel:
                                    setIntegraseCloseAndWithSimilarFamilyFromBlast.remove(currIntegraseCloseAndWithSimilarFamilyFromBlastToDel)

                                if len(setIntegraseCloseAndWithSimilarFamilyFromBlast) > 0:
                                    currStrComment = ""
                                    # if len(setPutativeIntegraseWithSimilarFamilyFromBlast) != 0:
                                    if currStrComment:
                                        currStrComment += ", "
                                    currStrComment += " " + hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseCloseAndWithSimilarFamilyFromBlast)  # setPutativeIntegraseWithSimilarFamilyFromBlast
                                    # else:
                                    #    raise RuntimeError("Error in addSPIntegraseUpstreamAndDownstream use family2SetICEsIMEsStructuresIdentifiers and family2SetPutativeIntLocusTag to complement info on integrase 1:"\
                                    #                       + "len(setPutativeIntegraseWithSimilarFamilyFromBlast) == 0 for currFamilyFromBlast {} and currICEsIMEsStructure.internalIdentifier {}".format(currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier))

                                    currStrCommentBis = "Integrase(s) previously seen associated to a conjugation module family but not directly upstream or downstream has/have been found: " + currStrComment + ". "
                                    currStrCommentBis += "Please manually check if this/those integrase(s) should be merged to the structure {}. ".format(currICEsIMEsStructure.internalIdentifier)
                                    if not onlyThisICEIMEStructIsRegistredForThisFamily:
                                        currStrCommentBis += "Please note that there are other ICE /IME structures with similar family: " + EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(family2SetICEsIMEsStructures[currFamilyFromBlast]) + ". "
                                    # for integraseToAddCommentIT in setIntegraseCloseAndWithSimilarFamilyFromBlast:
                                    #    icescreen_OO.addCommentToLocusTag2Comment(integraseToAddCommentIT.locusTag, currStrCommentBis, locusTagIntegrase2Comment)

                                    reprSetPutativeIntegraseWithSimilarFamilyFromBlast = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseCloseAndWithSimilarFamilyFromBlast)  # setPutativeIntegraseWithSimilarFamilyFromBlast
                                    if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
                                        setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
                                        if reprSetPutativeIntegraseWithSimilarFamilyFromBlast in setIntegraseManuallyCheck2comment:
                                            raise RuntimeError(
                                                    "Error in addSPIntegraseUpstreamAndDownstream reprSetPutativeIntegraseWithSimilarFamilyFromBlast: reprSetPutativeIntegraseWithSimilarFamilyFromBlast ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                                                            hit.ListSPs.GetListProtIdsFromSetSP(reprSetPutativeIntegraseWithSimilarFamilyFromBlast), currICEsIMEsStructure.internalIdentifier))
                                        else:
                                            setIntegraseManuallyCheck2comment[reprSetPutativeIntegraseWithSimilarFamilyFromBlast] = currStrCommentBis
                                    else:
                                        setIntegraseManuallyCheck2comment = {}
                                        setIntegraseManuallyCheck2comment[reprSetPutativeIntegraseWithSimilarFamilyFromBlast] = currStrCommentBis
                                        ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment

                        else:  # some Upstream and/or Downstream integrase found
                            setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel = set()
                            for currIntegraseCloseAndWithSimilarFamilyFromBlast in setIntegraseCloseAndWithSimilarFamilyFromBlast:
                                if currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag in setLocusTagToNotConsiderAseManuallyCheckSent:
                                    if DEBUG:
                                        print("setIntegraseCloseAndWithSimilarFamilyFromBlast some Upstream and/or Downstream integrase found: {}".format(
                                                currIntegraseCloseAndWithSimilarFamilyFromBlast.locusTag))
                                    setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel.add(currIntegraseCloseAndWithSimilarFamilyFromBlast)
                            for currIntegraseCloseAndWithSimilarFamilyFromBlastToDel in setCurrIntegraseCloseAndWithSimilarFamilyFromBlastToDel:
                                setIntegraseCloseAndWithSimilarFamilyFromBlast.remove(currIntegraseCloseAndWithSimilarFamilyFromBlastToDel)

                            if len(setIntegraseCloseAndWithSimilarFamilyFromBlast) > 0:

                                currStrComment = ""
                                if listUpstreamIntToAdd:
                                    currStrComment += "upstream integrase(s) " + hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)
                                if listDownstreamIntToAdd:
                                    if currStrComment:
                                        currStrComment += ", "
                                    currStrComment += "downstream integrase(s) " + hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)
                                # if len(setPutativeIntegraseWithSimilarFamilyFromBlast) != 0:
                                if currStrComment:
                                    currStrComment += ", "
                                currStrComment += "integrase previously seen associated to a conjugation module family but not directly upstream or downstream "\
                                    + hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseCloseAndWithSimilarFamilyFromBlast)  # setPutativeIntegraseWithSimilarFamilyFromBlast
                                # else:
                                #    raise RuntimeError("Error in addSPIntegraseUpstreamAndDownstream use family2SetICEsIMEsStructuresIdentifiers and family2SetPutativeIntLocusTag to complement info on integrase 2:"\
                                #                       + "len(setPutativeIntegraseWithSimilarFamilyFromBlast) == 0 for currFamilyFromBlast {} and currICEsIMEsStructure.internalIdentifier {}".format(currFamilyFromBlast, currICEsIMEsStructure.internalIdentifier))

                                setIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound = set()
                                for currInt in listUpstreamIntToAdd:
                                    # currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currInt.locusTag)
                                    setIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound.add(currInt)
                                for currInt in listDownstreamIntToAdd:
                                    # currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currInt.locusTag)
                                    setIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound.add(currInt)
                                for currInt in setIntegraseCloseAndWithSimilarFamilyFromBlast:  # setPutativeIntegraseWithSimilarFamilyFromBlast
                                    # currICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currInt.locusTag)
                                    setIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound.add(currInt)

                                currStrCommentBis = "Could not choose among the following candidate integrases: " + currStrComment + ". "
                                currStrCommentBis += "Please manually check if this/those integrase(s) should be added to the structure {}. ".format(currICEsIMEsStructure.internalIdentifier)
                                if not onlyThisICEIMEStructIsRegistredForThisFamily:
                                    currStrCommentBis += "Please note that there are other ICE /IME structures with similar family: " + EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(family2SetICEsIMEsStructures[currFamilyFromBlast]) + ". "

                                # for integraseToAddCommentIT in listUpstreamIntToAdd:
                                #    icescreen_OO.addCommentToLocusTag2Comment(integraseToAddCommentIT.locusTag, currStrCommentBis, locusTagIntegrase2Comment)
                                # for integraseToAddCommentIT in listDownstreamIntToAdd:
                                #    icescreen_OO.addCommentToLocusTag2Comment(integraseToAddCommentIT.locusTag, currStrCommentBis, locusTagIntegrase2Comment)
                                # for integraseToAddCommentIT in setIntegraseCloseAndWithSimilarFamilyFromBlast:
                                #    icescreen_OO.addCommentToLocusTag2Comment(integraseToAddCommentIT.locusTag, currStrCommentBis, locusTagIntegrase2Comment)

                                reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound)
                                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
                                    setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
                                    if reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound in setIntegraseManuallyCheck2comment:
                                        raise RuntimeError(
                                                "Error in addSPIntegraseUpstreamAndDownstream reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound: reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                                                        hit.ListSPs.GetListProtIdsFromSetSP(reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound), currICEsIMEsStructure.internalIdentifier))
                                    else:
                                        setIntegraseManuallyCheck2comment[reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound] = currStrCommentBis
                                else:
                                    setIntegraseManuallyCheck2comment = {}
                                    setIntegraseManuallyCheck2comment[reprSetIntegraseToManuallyCheckUpAndOrDownAndOrSimilarFamilyFound] = currStrCommentBis
                                    ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment

                                listUpstreamIntToAdd.clear()
                                listDownstreamIntToAdd.clear()

            if DEBUG:
                print(" - 2. listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
            if DEBUG:
                print(" - 2. listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            markIntegraseNeverPreviouslySeenAssociatedToConjugationModuleFamilyAsToManuallyCheck(
                    listUpstreamIntToAdd,
                    listDownstreamIntToAdd,
                    currICEsIMEsStructure,
                    locusTagIntegrase2Comment,
                    setLocusTagToNotConsiderAseManuallyCheckSent,
                    ICEsIMEsStructure2setIntegraseManuallyCheck2comment)

            if DEBUG:
                print(" - 3. listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
            if DEBUG:
                print(" - 3. listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            # print("HERE 0 listUpstreamIntToAdd and listDownstreamIntToAdd !!!!: {} and {} ".format(str(len(listUpstreamIntToAdd)), str(len(listDownstreamIntToAdd))))
            if listUpstreamIntToAdd and listDownstreamIntToAdd:  # remove integrase that should not be commented
                for currSp in reversed(listUpstreamIntToAdd):
                    if currSp.locusTag in setLocusTagToNotConsiderAseManuallyCheckSent:
                        listUpstreamIntToAdd.remove(currSp)
                for currSp in reversed(listDownstreamIntToAdd):
                    if currSp.locusTag in setLocusTagToNotConsiderAseManuallyCheckSent:
                        listDownstreamIntToAdd.remove(currSp)
            # print("HERE 1 listUpstreamIntToAdd and listDownstreamIntToAdd !!!!: {} and {} ".format(str(len(listUpstreamIntToAdd)), str(len(listDownstreamIntToAdd))))

            if DEBUG:
                print(" - 4. listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
            if DEBUG:
                print(" - 4. listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            # check if 1 integrase is way closer than the other
            if listUpstreamIntToAdd and listDownstreamIntToAdd:
                # and useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance > 0 and useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance > 0:

                # valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase returned:
                # 0: did not perform the test
                # 1: Downstream integrase is significantly further away to ICE/IME structure
                # 2: Downstream integrase is NOT significantly further away to ICE/IME structure
                # 3: Upstream integrase is significantly further away to ICE/IME structure
                # 4: Upstream integrase is NOT significantly further away to ICE/IME structure
                (valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase,
                 distCDSWithDownstreamIntegrase,
                 distCDSWithUpstreamIntegrase) = useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase(
                         listUpstreamIntToAdd,
                         listDownstreamIntToAdd,
                         currICEsIMEsStructure,
                         useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                         useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
                         True,
                         locusTagIntegrase2Comment)

            if listUpstreamIntToAdd and listDownstreamIntToAdd:
                listLocusTagInComment = []

                setIntegraseBothUpstreamAndDownstream = set()
                for currSp in reversed(listUpstreamIntToAdd):
                    setIntegraseBothUpstreamAndDownstream.add(currSp)
                    listUpstreamIntToAdd.remove(currSp)
                    listLocusTagInComment.append(currSp.locusTag)
                for currSp in reversed(listDownstreamIntToAdd):
                    setIntegraseBothUpstreamAndDownstream.add(currSp)
                    listDownstreamIntToAdd.remove(currSp)
                    listLocusTagInComment.append(currSp.locusTag)

                reprSetIntegraseBothUpstreamAndDownstream = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseBothUpstreamAndDownstream)
                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
                    setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
                    if reprSetIntegraseBothUpstreamAndDownstream in setIntegraseManuallyCheck2comment:
                        raise RuntimeError(
                                "Error in addSPIntegraseUpstreamAndDownstream reprSetIntegraseBothUpstreamAndDownstream: reprSetIntegraseBothUpstreamAndDownstream ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                                        hit.ListSPs.GetListProtIdsFromSetSP(reprSetIntegraseBothUpstreamAndDownstream), currICEsIMEsStructure.internalIdentifier))
                    else:
                        commentITToAdd = "Could not choose between upstream and downstream putative integrase: {}. Please manually check. ".format(
                                " and/or ".join(str(i) for i in sorted(listLocusTagInComment)))
                        setIntegraseManuallyCheck2comment[reprSetIntegraseBothUpstreamAndDownstream] = commentITToAdd
                        # for locusTagToAddCommentIT in listLocusTagInComment:
                        #    icescreen_OO.addCommentToLocusTag2Comment(locusTagToAddCommentIT, commentITToAdd, locusTagIntegrase2Comment)
                else:
                    setIntegraseManuallyCheck2comment = {}
                    commentITToAdd = "Could not choose between upstream and downstream putative integrase: {}. Please manually check. ".format(
                            " and/or ".join(str(i) for i in sorted(listLocusTagInComment)))
                    setIntegraseManuallyCheck2comment[reprSetIntegraseBothUpstreamAndDownstream] = commentITToAdd
                    # for locusTagToAddCommentIT in listLocusTagInComment:
                    #    icescreen_OO.addCommentToLocusTag2Comment(locusTagToAddCommentIT, commentITToAdd, locusTagIntegrase2Comment)

                    ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment

            elif listUpstreamIntToAdd:

                for currSp in listUpstreamIntToAdd:

                    if currSp in integraseAttributedSureUpstream2setICEsIMEsStructure:
                        currSetICEsIMEsStructure = integraseAttributedSureUpstream2setICEsIMEsStructure[currSp]
                        currSetICEsIMEsStructure.add(currICEsIMEsStructure)
                    else:
                        currSetICEsIMEsStructure = set()
                        currSetICEsIMEsStructure.add(currICEsIMEsStructure)
                        integraseAttributedSureUpstream2setICEsIMEsStructure[currSp] = currSetICEsIMEsStructure

                    if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseMoreSure:
                        setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure]
                        setIntegraseMoreSure.add(currSp)
                    else:
                        setIntegraseMoreSure = set()
                        setIntegraseMoreSure.add(currSp)
                        ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure] = setIntegraseMoreSure

            elif listDownstreamIntToAdd:

                for currSp in listDownstreamIntToAdd:

                    if currSp in integraseAttributedSureDownstream2setICEsIMEsStructure:
                        currSetICEsIMEsStructure = integraseAttributedSureDownstream2setICEsIMEsStructure[currSp]
                        currSetICEsIMEsStructure.add(currICEsIMEsStructure)
                    else:
                        currSetICEsIMEsStructure = set()
                        currSetICEsIMEsStructure.add(currICEsIMEsStructure)
                        integraseAttributedSureDownstream2setICEsIMEsStructure[currSp] = currSetICEsIMEsStructure

                    if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseMoreSure:
                        setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure]
                        setIntegraseMoreSure.add(currSp)
                    else:
                        setIntegraseMoreSure = set()
                        setIntegraseMoreSure.add(currSp)
                        ICEsIMEsStructure2setIntegraseMoreSure[currICEsIMEsStructure] = setIntegraseMoreSure

    if doNotTakeIntoAccountIfManuallyCheckButAnIntegraseMoreSureIsRegistred:
        setLocusTagToNotConsiderAseManuallyCheck = set()
        for keyICEsIMEsStructureManuallyCheck, valueReprSetIntegraseManuallyCheck2comment in ICEsIMEsStructure2setIntegraseManuallyCheck2comment.items():
            for keyReprSetIntegraseManuallyCheck, valueComment in valueReprSetIntegraseManuallyCheck2comment.items():
                # keyICEsIMEsStructureManuallyCheck.comment += valueComment
                listIntegraseLocusTag = keyReprSetIntegraseManuallyCheck.split(", ")
                for currIntegraseLocusTagManuallyCheck in listIntegraseLocusTag:

                    for keyICEsIMEsStructureMoreSure, valueSetIntegraseMoreSure in ICEsIMEsStructure2setIntegraseMoreSure.items():
                        for currIntegraseMoreSure in valueSetIntegraseMoreSure:
                            if currIntegraseMoreSure.locusTag == currIntegraseLocusTagManuallyCheck:
                                # print("keyReprSetIntegraseManuallyCheck {}: collusion locus tag {} and keyICEsIMEsStructureManuallyCheck {} - keyICEsIMEsStructureMoreSure {}"\
                                #      .format(repr(keyReprSetIntegraseManuallyCheck), currIntegraseMoreSure.locusTag, keyICEsIMEsStructureManuallyCheck.internalIdentifier, keyICEsIMEsStructureMoreSure.internalIdentifier))
                                setLocusTagToNotConsiderAseManuallyCheck.add(currIntegraseLocusTagManuallyCheck)
                                break

        setUnionLocusTagToNotConsiderAseManuallyCheck = setLocusTagToNotConsiderAseManuallyCheckSent.union(setLocusTagToNotConsiderAseManuallyCheck)
        if len(setLocusTagToNotConsiderAseManuallyCheck) >= 0 and countIterAddSPIntegraseUpstreamAndDownstream < 5 and len(setUnionLocusTagToNotConsiderAseManuallyCheck) > len(setLocusTagToNotConsiderAseManuallyCheckSent):
            if DEBUG:
                print("\n !!! rerun addSPIntegraseUpstreamAndDownstream with new setUnionLocusTagToNotConsiderAseManuallyCheck = {}\n".format(
                        setUnionLocusTagToNotConsiderAseManuallyCheck))

            addSPIntegraseUpstreamAndDownstream(
                    listSPs,
                    listICEsIMEsStructures,
                    maxNumberCDSForFilterIMESize,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    setUnionLocusTagToNotConsiderAseManuallyCheck,
                    countIterAddSPIntegraseUpstreamAndDownstream,
                    allowAdjacentIntegraseOnlyForSer,
                    locusTagIntegrase2Comment,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
            return

    for keyIntegraseAttributedSureUpstream, valueSetICEsIMEsStructure in integraseAttributedSureUpstream2setICEsIMEsStructure.items():
        currSetICEsIMEsStructure = set()
        currSetICEsIMEsStructure.update(valueSetICEsIMEsStructure)
        if keyIntegraseAttributedSureUpstream in integraseAttributedSureDownstream2setICEsIMEsStructure:
            currSetICEsIMEsStructure.update(integraseAttributedSureDownstream2setICEsIMEsStructure[keyIntegraseAttributedSureUpstream])
            del integraseAttributedSureDownstream2setICEsIMEsStructure[keyIntegraseAttributedSureUpstream]
        if len(currSetICEsIMEsStructure) >= 2:
            # move to To check
            setIntegraseMoveToCheck = set()
            setIntegraseMoveToCheck.add(keyIntegraseAttributedSureUpstream)

            reprSetIntegraseMoveToCheck = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseMoveToCheck)
            for currICEsIMEsStructure in currSetICEsIMEsStructure:
                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
                    setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
                    if reprSetIntegraseMoveToCheck == setIntegraseManuallyCheck2comment:
                        raise RuntimeError(
                                "Error in addSPIntegraseUpstreamAndDownstream reprSetIntegraseMoveToCheck: reprSetIntegraseMoveToCheck ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                                        hit.ListSPs.GetListProtIdsFromSetSP(reprSetIntegraseMoveToCheck), currICEsIMEsStructure.internalIdentifier))
                    else:
                        commentITToAdd = "Integrase {} could be attributed to multiple ICE IME structures ({}), please manually check. ".format(
                                keyIntegraseAttributedSureUpstream.locusTag, EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                        setIntegraseManuallyCheck2comment[reprSetIntegraseMoveToCheck] = commentITToAdd
                        # icescreen_OO.addCommentToLocusTag2Comment(keyIntegraseAttributedSureUpstream.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                else:
                    setIntegraseManuallyCheck2comment = {}

                    commentITToAdd = "Integrase {} could be attributed to multiple ICE IME structures ({}), please manually check. ".format(
                            keyIntegraseAttributedSureUpstream.locusTag, EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                    setIntegraseManuallyCheck2comment[reprSetIntegraseMoveToCheck] = commentITToAdd
                    # icescreen_OO.addCommentToLocusTag2Comment(keyIntegraseAttributedSureUpstream.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                    ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment
        else:
            # register integrase to ICEsIMEsStructure
            for currValueICEsIMEsStructure in currSetICEsIMEsStructure:
                currValueICEsIMEsStructure.listIntegraseUpstream.append(keyIntegraseAttributedSureUpstream)
                currValueICEsIMEsStructure.listOrderedSPs.append(keyIntegraseAttributedSureUpstream)
                currValueICEsIMEsStructure.refreshListIdxOrderedSPs()
                break

    for keyIntegraseAttributedSureDownstream, valueSetICEsIMEsStructure in integraseAttributedSureDownstream2setICEsIMEsStructure.items():
        currSetICEsIMEsStructure = set()
        currSetICEsIMEsStructure.update(valueSetICEsIMEsStructure)
        if keyIntegraseAttributedSureDownstream in integraseAttributedSureUpstream2setICEsIMEsStructure:
            currSetICEsIMEsStructure.update(integraseAttributedSureUpstream2setICEsIMEsStructure[keyIntegraseAttributedSureDownstream])
            del integraseAttributedSureUpstream2setICEsIMEsStructure[keyIntegraseAttributedSureDownstream]
        if len(currSetICEsIMEsStructure) >= 2:
            # move to To check
            setIntegraseMoveToCheck = set()
            setIntegraseMoveToCheck.add(keyIntegraseAttributedSureDownstream)
            reprSetIntegraseMoveToCheck = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseMoveToCheck)
            for currICEsIMEsStructure in currSetICEsIMEsStructure:
                if currICEsIMEsStructure in ICEsIMEsStructure2setIntegraseManuallyCheck2comment:
                    setIntegraseManuallyCheck2comment = ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure]
                    if reprSetIntegraseMoveToCheck == setIntegraseManuallyCheck2comment:
                        raise RuntimeError(
                                "Error in addSPIntegraseUpstreamAndDownstream reprSetIntegraseMoveToCheck: reprSetIntegraseMoveToCheck ({}) in setIntegraseManuallyCheck2comment for currICEsIMEsStructure.internalIdentifier {}".format(
                                        hit.ListSPs.GetListProtIdsFromSetSP(reprSetIntegraseMoveToCheck), currICEsIMEsStructure.internalIdentifier))
                    else:

                        commentITToAdd = "Integrase {} could be attributed to multiple ICE IME structures ({}), please manually check. ".format(
                                keyIntegraseAttributedSureDownstream.locusTag, EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                        setIntegraseManuallyCheck2comment[reprSetIntegraseMoveToCheck] = commentITToAdd
                        # icescreen_OO.addCommentToLocusTag2Comment(keyIntegraseAttributedSureDownstream.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                else:
                    setIntegraseManuallyCheck2comment = {}
                    commentITToAdd = "Integrase {} could be attributed to multiple ICE IME structures ({}), please manually check. ".format(
                            keyIntegraseAttributedSureDownstream.locusTag, EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                    setIntegraseManuallyCheck2comment[reprSetIntegraseMoveToCheck] = commentITToAdd
                    # icescreen_OO.addCommentToLocusTag2Comment(keyIntegraseAttributedSureDownstream.locusTag, commentITToAdd, locusTagIntegrase2Comment)
                    ICEsIMEsStructure2setIntegraseManuallyCheck2comment[currICEsIMEsStructure] = setIntegraseManuallyCheck2comment
        else:
            # register integrase to ICEsIMEsStructure
            for currValueICEsIMEsStructure in currSetICEsIMEsStructure:
                currValueICEsIMEsStructure.listIntegraseDownstream.append(keyIntegraseAttributedSureDownstream)
                currValueICEsIMEsStructure.listOrderedSPs.append(keyIntegraseAttributedSureDownstream)
                currValueICEsIMEsStructure.refreshListIdxOrderedSPs()
                break

    # ICEsIMEsStructure2setIntegraseManuallyCheck2comment
    for keyICEsIMEsStructure, valueReprSetIntegraseManuallyCheck2comment in ICEsIMEsStructure2setIntegraseManuallyCheck2comment.items():
        for keyReprSetIntegraseManuallyCheck, valueComment in valueReprSetIntegraseManuallyCheck2comment.items():
            keyICEsIMEsStructure.comment += valueComment
            listIntegraseLocusTag = keyReprSetIntegraseManuallyCheck.split(", ")
            for currIntegraseLocusTagManuallyCheck in listIntegraseLocusTag:
                keyICEsIMEsStructure.setIntegraseLocusTagsToManuallyCheck.add(currIntegraseLocusTagManuallyCheck)
                icescreen_OO.addCommentToLocusTag2Comment(
                        currIntegraseLocusTagManuallyCheck,
                        valueComment,
                        locusTagIntegrase2Comment)

    return locusTagIntegrase2Comment
