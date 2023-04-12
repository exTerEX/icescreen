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

#structureWithBothUpstreamAndDownstreamIntegraseCanNotChoose2structure = {}
dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd = {}


# noqa: E501 This method correct the wrongly attribution of a subsequent integrase that was initially thought to be an obvious choice for a conjugation module
def changeObviousIntegraseAttributionToUnsureBecauseOfBothUpstreamAndDownstreamEqualyPossible(
        currICEsIMEsStructure,
        otherICEsIMEsStructureToMerge,
        locusTagIntegrase2Comment):

    listUpstreamIntChanged = []
    listDownstreamIntChanged = []

    strCommentIT = "Following a merge, both upstream integrase {} and downstream integrase {} could possibility be associated to the conjugaison module family of structure {}. ".format(
            hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseUpstream),
            hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseDownstream),
            currICEsIMEsStructure.internalIdentifier)
    if strCommentIT not in currICEsIMEsStructure.comment:
        currICEsIMEsStructure.comment += strCommentIT
    for currSp in currICEsIMEsStructure.listIntegraseUpstream:
        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
        # take off key and clean up integrase and comment in already registred structure
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
        if strCommentToRemove in currICEsIMEsStructure.comment:
            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, otherICEsIMEsStructureToMerge.internalIdentifier)
        if strCommentToRemove in otherICEsIMEsStructureToMerge.comment:
            otherICEsIMEsStructureToMerge.comment = otherICEsIMEsStructureToMerge.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)

        # currICEsIMEsStructure.setIntegraseToManuallyCheck.add(currSp.locusTag)
        listUpstreamIntChanged.append(currSp)
        if currSp in currICEsIMEsStructure.listOrderedSPs:
            currICEsIMEsStructure.listOrderedSPs.remove(currSp)
    currICEsIMEsStructure.listIntegraseUpstream.clear()

    for currSp in currICEsIMEsStructure.listIntegraseDownstream:
        icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
        # take off key and clean up integrase and comment in already registred structure
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
        if strCommentToRemove in currICEsIMEsStructure.comment:
            currICEsIMEsStructure.comment = currICEsIMEsStructure.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)
        strCommentToRemove = "Integrase {} has been associated to the structure {} because they are adjacent and there is no upstream/downstream ambiguity. ".format(
                currSp.locusTag, otherICEsIMEsStructureToMerge.internalIdentifier)
        if strCommentToRemove in otherICEsIMEsStructureToMerge.comment:
            otherICEsIMEsStructureToMerge.comment = otherICEsIMEsStructureToMerge.comment.replace(strCommentToRemove, "")
        icescreen_OO.removeCommentToLocusTag2Comment(currSp.locusTag, strCommentToRemove, locusTagIntegrase2Comment)

        # currICEsIMEsStructure.setIntegraseToManuallyCheck.add(currSp.locusTag)
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

    # print("{} isIntegraseCorrectlyOrientedForICEsIMEsStructure: integraseSent = {} (len(setVirB4InConflictInICEsIMEsStructureSent) = {}) ; ICEsIMEsStructureSent = {} ; integraseIsUpstreamOfICEsIMEsStructure = {} ".format(
    #             repr(valueToReturn), integraseSent.locusTag, len(setVirB4InConflictInICEsIMEsStructureSent), ICEsIMEsStructureSent.internalIdentifier, integraseIsUpstreamOfICEsIMEsStructure))

    # if valueToReturn == -1:
    #    raise RuntimeError("Error in isIntegraseCorrectlyOrientedForICEsIMEsStructure: valueToReturn == -1 ; integraseSent = {} (len(setVirB4InConflictInICEsIMEsStructureSent) = {}) ; ICEsIMEsStructureSent = {} ; integraseIsUpstreamOfICEsIMEsStructure = {} " \
    #          .format(integraseSent.locusTag, len(setVirB4InConflictInICEsIMEsStructureSent), ICEsIMEsStructureSent.internalIdentifier, integraseIsUpstreamOfICEsIMEsStructure))

    return valueToReturn


def addObviousIntegraseUpstreamAndDownstream_priorMerging(
        listICEsIMEsStructures,
        currListSPs,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance,
        allowAdjacentIntegraseOnlyForSer
    ):

    global dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd
    for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):

        #print("NEW addObviousIntegraseUpstreamAndDownstream_priorMerging for currICEsIMEsStructure = {} ({})".format(currICEsIMEsStructure.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listOrderedSPs)))

        upstreamICEsIMEsStructure = None
        downstreamICEsIMEsStructure = None
        if idxCurrICEsIMEsStructure > 0:
            upstreamICEsIMEsStructure = listICEsIMEsStructures[idxCurrICEsIMEsStructure-1]
        if idxCurrICEsIMEsStructure < (len(listICEsIMEsStructures) - 1):
            downstreamICEsIMEsStructure  = listICEsIMEsStructures[idxCurrICEsIMEsStructure+1]

        idxlastSP_currICEsIMEsStructure = hit.ListSPs.GetIdxSPInList(currICEsIMEsStructure.listOrderedSPs[len(currICEsIMEsStructure.listOrderedSPs)-1], currListSPs)
        (listUpstreamObviousIntegrase_currICEsIMEsStructure, listDownstreamObviousIntegrase_currICEsIMEsStructure) = checkForObviousIntegraseUpstreamAndDownstreamToAdd(
                currICEsIMEsStructure,
                currListSPs,
                idxlastSP_currICEsIMEsStructure,
                allowAdjacentIntegraseOnlyForSer,
                locusTagIntegrase2Comment
                )

        listUpstreamObviousIntegrase_upstreamICEsIMEsStructure = None
        listDownstreamObviousIntegrase_upstreamICEsIMEsStructure = None
        if upstreamICEsIMEsStructure is not None:
            idxlastSP_upstreamICEsIMEsStructure = hit.ListSPs.GetIdxSPInList(upstreamICEsIMEsStructure.listOrderedSPs[len(upstreamICEsIMEsStructure.listOrderedSPs)-1], currListSPs)
            (listUpstreamObviousIntegrase_upstreamICEsIMEsStructure, listDownstreamObviousIntegrase_upstreamICEsIMEsStructure) = checkForObviousIntegraseUpstreamAndDownstreamToAdd(
                upstreamICEsIMEsStructure,
                currListSPs,
                idxlastSP_upstreamICEsIMEsStructure,
                allowAdjacentIntegraseOnlyForSer,
                locusTagIntegrase2Comment
                )

        listUpstreamObviousIntegrase_downstreamICEsIMEsStructure = None
        listDownstreamObviousIntegrase_downstreamICEsIMEsStructure = None
        if downstreamICEsIMEsStructure is not None:
            idxlastSP_downstreamICEsIMEsStructure = hit.ListSPs.GetIdxSPInList(downstreamICEsIMEsStructure.listOrderedSPs[len(downstreamICEsIMEsStructure.listOrderedSPs)-1], currListSPs)
            (listUpstreamObviousIntegrase_downstreamICEsIMEsStructure, listDownstreamObviousIntegrase_downstreamICEsIMEsStructure) = checkForObviousIntegraseUpstreamAndDownstreamToAdd(
                downstreamICEsIMEsStructure,
                currListSPs,
                idxlastSP_downstreamICEsIMEsStructure,
                allowAdjacentIntegraseOnlyForSer,
                locusTagIntegrase2Comment
                )

        if len(listUpstreamObviousIntegrase_currICEsIMEsStructure) > 0 and len(listDownstreamObviousIntegrase_currICEsIMEsStructure) > 0 :
            # no obvious integrase as both upstream and downstream
            pass
        elif len(listUpstreamObviousIntegrase_currICEsIMEsStructure) > 0:
            setSPToRemove = set()
            for currSp in listUpstreamObviousIntegrase_currICEsIMEsStructure:
                # if upstream structure shares also this obvious integrase, do not take it into account
                if listDownstreamObviousIntegrase_upstreamICEsIMEsStructure is not None and len(listDownstreamObviousIntegrase_upstreamICEsIMEsStructure) > 0 :
                    if currSp in listDownstreamObviousIntegrase_upstreamICEsIMEsStructure :
                        # upstream structure shares also this obvious integrase, do not take it into account
                        setSPToRemove.add(currSp)

            # check if downstream structure has no integrase at all around it, if so this structure could be a nested structure between the module conj and integrase of the structure downstream, do not take it into account
            if downstreamICEsIMEsStructure :
                if not listUpstreamObviousIntegrase_downstreamICEsIMEsStructure and not listDownstreamObviousIntegrase_downstreamICEsIMEsStructure :
                    for currSp in listUpstreamObviousIntegrase_currICEsIMEsStructure:
                        setSPToRemove.add(currSp)
                        # strCommentIT = "Upstream integrase {} was not associated to the structure {} pre-merging because there is no integrase around the downstream structure {}, therefore there is a possibility that the structure {} is nested between the conj module of structure {} and the integrase {}. ".format(
                        #         currSp.locusTag, currICEsIMEsStructure.internalIdentifier, downstreamICEsIMEsStructure.internalIdentifier, currICEsIMEsStructure.internalIdentifier, downstreamICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                        # if strCommentIT not in currICEsIMEsStructure.comment:
                        #     currICEsIMEsStructure.comment += strCommentIT
                        # icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)

            for SPToRemoveIT in setSPToRemove :
                listUpstreamObviousIntegrase_currICEsIMEsStructure.remove(SPToRemoveIT)

            # register the integrase as obvious if it survives the tests before
            if len(listUpstreamObviousIntegrase_currICEsIMEsStructure) > 0:
                for currSp in listUpstreamObviousIntegrase_currICEsIMEsStructure :
                    strCommentIT = "Upstream integrase {} has been associated to the structure {} pre-merging because it is adjacent to the conj module and there is no upstream/downstream ambiguity. ".format(
                            currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd[currSp] = currICEsIMEsStructure
                    currICEsIMEsStructure.listIntegraseUpstream.append(currSp)
                    currICEsIMEsStructure.listOrderedSPs.append(currSp)
                    currICEsIMEsStructure.refreshListIdxOrderedSPs()


        elif len(listDownstreamObviousIntegrase_currICEsIMEsStructure) > 0:
            setSPToRemove = set()
            for currSp in listDownstreamObviousIntegrase_currICEsIMEsStructure:
                # if downstream structure shares also this obvious integrase, do not take it into account
                if listUpstreamObviousIntegrase_downstreamICEsIMEsStructure is not None and len(listUpstreamObviousIntegrase_downstreamICEsIMEsStructure) > 0 :
                    if currSp in listUpstreamObviousIntegrase_downstreamICEsIMEsStructure :
                        # upstream structure shares also this obvious integrase, do not take it into account
                        setSPToRemove.add(currSp)

            #check if upstream structure has no integrase at all around it, if so this structure could be a nested structure between the module conj and integrase of the structure upstream, do not take it into account
            if upstreamICEsIMEsStructure:
                if not listUpstreamObviousIntegrase_upstreamICEsIMEsStructure and not listDownstreamObviousIntegrase_upstreamICEsIMEsStructure :
                    for currSp in listDownstreamObviousIntegrase_currICEsIMEsStructure:
                        setSPToRemove.add(currSp)
                        # strCommentIT = "Downstream integrase {} was not associated to the structure {} pre-merging because there is no integrase around the upstream structure {}, therefore there is a possibility that the structure {} is nested between the conj module of structure {} and the integrase {}. ".format(
                        #         currSp.locusTag, currICEsIMEsStructure.internalIdentifier, upstreamICEsIMEsStructure.internalIdentifier, currICEsIMEsStructure.internalIdentifier, upstreamICEsIMEsStructure.internalIdentifier, currSp.locusTag)
                        # if strCommentIT not in currICEsIMEsStructure.comment:
                        #     currICEsIMEsStructure.comment += strCommentIT
                        # icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)

            for SPToRemoveIT in setSPToRemove :
                listDownstreamObviousIntegrase_currICEsIMEsStructure.remove(SPToRemoveIT)
            #register the integrase if it survives the tests before
            if len(listDownstreamObviousIntegrase_currICEsIMEsStructure) > 0:
                for currSp in listDownstreamObviousIntegrase_currICEsIMEsStructure :
                    dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd[currSp] = currICEsIMEsStructure
                    strCommentIT = "Downstream integrase {} has been associated to the structure {} pre-merging because it is adjacent to the conj module and there is no upstream/downstream ambiguity. ".format(
                            currSp.locusTag, currICEsIMEsStructure.internalIdentifier)
                    if strCommentIT not in currICEsIMEsStructure.comment:
                        currICEsIMEsStructure.comment += strCommentIT
                    icescreen_OO.addCommentToLocusTag2Comment(currSp.locusTag, strCommentIT, locusTagIntegrase2Comment)
                    currICEsIMEsStructure.listIntegraseDownstream.append(currSp)
                    currICEsIMEsStructure.listOrderedSPs.append(currSp)
                    currICEsIMEsStructure.refreshListIdxOrderedSPs()




# look for upstream or downstream subsequent integrase that have been preivously seen associated with the conj module
def checkForObviousIntegraseUpstreamAndDownstreamToAdd(
        EMStructureToAdd,
        listSPs,
        idxCurrentSP,
        allowAdjacentIntegraseOnlyForSer,
        locusTagIntegrase2Comment
        # useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        # useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance
        ):

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
        #print("listSPsUpstream for EMStructureToAdd {} : {}".format(EMStructureToAdd.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listSPsUpstream)))
        for currSp in reversed(listSPsUpstream):
            # print("HERE up {}".format(currSp.locusTag))
            if listUpstreamIntToAdd:
                #should be consistent with hit.getListIntegraseGroupJustUpstreamOfThisSP
                #TODO: TEST also consider fragment as integrase that follow up
                if (abs(listUpstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2) or (listUpstreamIntToAdd[-1] in currSp.listSiblingFragmentedSP and currSp in listUpstreamIntToAdd[-1].listSiblingFragmentedSP):  # Integrases are adjacent in genome or separated by a CDS
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
                        listUpstreamIntToAdd.append(currSp)
                        continue
                        # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #    listUpstreamIntToAdd.append(currSp)
                        #    continue
                        # if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #     listUpstreamIntToAdd.append(currSp)
                        #     continue
                        # elif isIMEFamillyIntegraseSameAsIMESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #     listUpstreamIntToAdd.append(currSp)
                        #     continue
                        # else:
                        #     break
                else:
                    break

    listDownstreamIntToAdd = []
    idxSPsDownstreamToCheck = idxCurrentSP + 1
    if idxSPsDownstreamToCheck < len(listSPs):
        listSPsDownstream = listSPs[idxSPsDownstreamToCheck:]
        #print("listSPsDownstream for EMStructureToAdd {} : {}".format(EMStructureToAdd.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listSPsDownstream)))
        for currSp in listSPsDownstream:
            #print("HERE down {}".format(currSp.locusTag))
            if listDownstreamIntToAdd:  # listIntegraseDownstream is not empty
                #TODO: TEST also consider fragment as integrase that follow up
                if (abs(listDownstreamIntToAdd[-1].CDSPositionInGenome - currSp.CDSPositionInGenome) <= 2) or (listDownstreamIntToAdd[-1] in currSp.listSiblingFragmentedSP and currSp in listDownstreamIntToAdd[-1].listSiblingFragmentedSP):  # Integrases are adjacent in genome or separated by a CDS
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
                        #print("\tHERE listDownstreamIntToAdd.append(currSp)")

                        listDownstreamIntToAdd.append(currSp)
                        continue
                        # if isICESuperFamilyIntegraseSameAsICESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #    listDownstreamIntToAdd.append(currSp)
                        #    continue
                        # if isICEFamillyIntegraseSameAsICEFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #     listDownstreamIntToAdd.append(currSp)
                        #     continue
                        # elif isIMEFamillyIntegraseSameAsIMESuperFamilyConjModule(currSp, EMStructureToAdd) == 2:
                        #     listDownstreamIntToAdd.append(currSp)
                        #     continue
                        # else:
                        #     break
                else:
                    break
    return (listUpstreamIntToAdd, listDownstreamIntToAdd)


def addEntryToICEsIMEsStructure2IntegraseManuallyCheck2comment(ICEsIMEsStructureIT, setIntegraseIT, commentITToAdd, ICEsIMEsStructure2IntegraseManuallyCheck2comment) :

    hashIntegraseManuallyCheck2comment = {}
    if ICEsIMEsStructureIT in ICEsIMEsStructure2IntegraseManuallyCheck2comment:
        hashIntegraseManuallyCheck2comment = ICEsIMEsStructure2IntegraseManuallyCheck2comment[ICEsIMEsStructureIT]
    for integraseIT in setIntegraseIT :
        if integraseIT in hashIntegraseManuallyCheck2comment :
            commentAlreadyInhashIntegraseManuallyCheck2comment = hashIntegraseManuallyCheck2comment[integraseIT]
            if commentITToAdd not in commentAlreadyInhashIntegraseManuallyCheck2comment :
                newCommentToPutInhashIntegraseManuallyCheck2comment = commentAlreadyInhashIntegraseManuallyCheck2comment + commentITToAdd
                hashIntegraseManuallyCheck2comment[integraseIT] = newCommentToPutInhashIntegraseManuallyCheck2comment
        else:
            hashIntegraseManuallyCheck2comment[integraseIT] = commentITToAdd
    if ICEsIMEsStructureIT not in ICEsIMEsStructure2IntegraseManuallyCheck2comment:
        ICEsIMEsStructure2IntegraseManuallyCheck2comment[ICEsIMEsStructureIT] = hashIntegraseManuallyCheck2comment

        


def integraseCouldBeAttributedToMultipleICEIMEStructures(integraseIT, setMultipleICEsIMEsStructure, possibleNested, locusTagIntegrase2Comment): #, ICEsIMEsStructure2IntegraseManuallyCheck2comment

    strCalrifyComment = "adjacent ICE IME structures"
    if possibleNested:
        strCalrifyComment = "possibly nested ICE IME structures"

    commentITToAdd = "Integrase {} could be attributed to multiple ".format(integraseIT.locusTag) + strCalrifyComment + " ({}), please manually check. ".format(
                        EMStructure.BasicEMStructure.GetListInternIdFromSetEMStructure(setMultipleICEsIMEsStructure))
    icescreen_OO.addCommentToLocusTag2Comment(
            integraseIT.locusTag,
            commentITToAdd,
            locusTagIntegrase2Comment)
    for currICEsIMEsStructure in setMultipleICEsIMEsStructure:
        if commentITToAdd not in currICEsIMEsStructure.comment:
            currICEsIMEsStructure.comment += commentITToAdd
        # currICEsIMEsStructure.setIntegraseToManuallyCheck.add(integraseIT.locusTag)
        currICEsIMEsStructure.setIntegraseToManuallyCheck.add(integraseIT)


def addStructuresChainedByIntegrase(setIntersectionCommonKeysSent, outer_integraseAttributed2setICEsIMEsStructure, inner_integraseAttributed2setICEsIMEsStructure, listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturnIT):

    for commonKeysBetweenSetsSureUpstreamAndSureDownstreamIT in setIntersectionCommonKeysSent :
        valueSetICEsIMEsStructureOuter = outer_integraseAttributed2setICEsIMEsStructure[commonKeysBetweenSetsSureUpstreamAndSureDownstreamIT]
        valueSetICEsIMEsStructureInner = inner_integraseAttributed2setICEsIMEsStructure[commonKeysBetweenSetsSureUpstreamAndSureDownstreamIT]
        symmetricDifferenceIT = valueSetICEsIMEsStructureOuter.symmetric_difference(valueSetICEsIMEsStructureInner)
        if len(symmetricDifferenceIT) > 0:
            foundASetChainsOfStructuresLinkedBySharedIntegraseITWithAlreadyOneOfOurStructureInIt = False
            for setChainsOfStructuresLinkedBySharedIntegraseIT in listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturnIT :
                # intersect setChainsOfStructuresLinkedBySharedIntegraseIT with valueSetICEsIMEsStructureOuter and valueSetICEsIMEsStructureInner
                intersectWithOuter = setChainsOfStructuresLinkedBySharedIntegraseIT.intersection(valueSetICEsIMEsStructureOuter) 
                intersectWithInner = setChainsOfStructuresLinkedBySharedIntegraseIT.intersection(valueSetICEsIMEsStructureInner)
                if len(intersectWithOuter) > 0 or len(intersectWithInner) > 0 :
                    foundASetChainsOfStructuresLinkedBySharedIntegraseITWithAlreadyOneOfOurStructureInIt = True
                    #setUnionIT = setChainsOfStructuresLinkedBySharedIntegraseIT.union(valueSetICEsIMEsStructureOuter, valueSetICEsIMEsStructureInner)
                    setChainsOfStructuresLinkedBySharedIntegraseIT.update(valueSetICEsIMEsStructureOuter)
                    setChainsOfStructuresLinkedBySharedIntegraseIT.update(valueSetICEsIMEsStructureInner)
            if not foundASetChainsOfStructuresLinkedBySharedIntegraseITWithAlreadyOneOfOurStructureInIt :
                # store this new chain
                setUnionIT = valueSetICEsIMEsStructureOuter.union(valueSetICEsIMEsStructureInner)
                listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturnIT.append(setUnionIT)


def createListOfListChainsOfStructuresLinkedBySharedIntegrase(integraseAttributedSureUpstream2setICEsIMEsStructure, integraseAttributedSureDownstream2setICEsIMEsStructure, IntegraseManuallyCheckUpstream2setICEsIMEsStructure, IntegraseManuallyCheckDownstream2setICEsIMEsStructure, DEBUG):
    listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn = []

    if DEBUG:
        print("\n ** Starting createListOfListChainsOfStructuresLinkedBySharedIntegrase ")

    # get all the integrase keys and find the one that are in more than 1 dictionary
    setKeysSureUpstream = set(integraseAttributedSureUpstream2setICEsIMEsStructure)
    setKeysSureDownstream = set(integraseAttributedSureDownstream2setICEsIMEsStructure)
    setKeysManuallyUpstream = set(IntegraseManuallyCheckUpstream2setICEsIMEsStructure)
    setKeysManuallyDownstream = set(IntegraseManuallyCheckDownstream2setICEsIMEsStructure)
    
    setIntersectionCommonKeysBetweenSetsSureUpstreamAndSureDownstream = setKeysSureUpstream.intersection(setKeysSureDownstream)
    setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyUpstream = setKeysSureUpstream.intersection(setKeysManuallyUpstream)
    setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyDownstream = setKeysSureUpstream.intersection(setKeysManuallyDownstream)
    setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyUpstream = setKeysSureDownstream.intersection(setKeysManuallyUpstream)
    setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyDownstream = setKeysSureDownstream.intersection(setKeysManuallyDownstream)
    setIntersectionCommonKeysBetweenSetKeysManuallyUpstreamAndSetKeysManuallyDownstream = setKeysManuallyUpstream.intersection(setKeysManuallyDownstream)

    if len(setIntersectionCommonKeysBetweenSetsSureUpstreamAndSureDownstream) > 0 :

        addStructuresChainedByIntegrase(setIntersectionCommonKeysBetweenSetsSureUpstreamAndSureDownstream, integraseAttributedSureUpstream2setICEsIMEsStructure, integraseAttributedSureDownstream2setICEsIMEsStructure, listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn)

    if len(setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyUpstream) > 0 :
        raise RuntimeError("Error in createListOfListChainsOfStructuresLinkedBySharedIntegrase: len(setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyUpstream) > 0 for {}".format(hit.ListSPs.GetListProtIdsFromSetSP(setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyUpstream)))

    if len(setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyDownstream) > 0 :
        addStructuresChainedByIntegrase(setIntersectionCommonKeysBetweenSetKeysSureUpstreamAndSetKeysManuallyDownstream, integraseAttributedSureUpstream2setICEsIMEsStructure, IntegraseManuallyCheckDownstream2setICEsIMEsStructure, listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn)

    if len(setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyUpstream) > 0 :
        addStructuresChainedByIntegrase(setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyUpstream, integraseAttributedSureDownstream2setICEsIMEsStructure, IntegraseManuallyCheckUpstream2setICEsIMEsStructure, listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn)

    if len(setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyDownstream) > 0 :
        raise RuntimeError("Error in createListOfListChainsOfStructuresLinkedBySharedIntegrase: (setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyDownstream) > 0 : for {}".format(hit.ListSPs.GetListProtIdsFromSetSP(setIntersectionCommonKeysBetweenSetKeysSureDownstreamAndSetKeysManuallyDownstream)))

    if len(setIntersectionCommonKeysBetweenSetKeysManuallyUpstreamAndSetKeysManuallyDownstream) > 0 :
        addStructuresChainedByIntegrase(setIntersectionCommonKeysBetweenSetKeysManuallyUpstreamAndSetKeysManuallyDownstream, IntegraseManuallyCheckUpstream2setICEsIMEsStructure, IntegraseManuallyCheckDownstream2setICEsIMEsStructure, listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn)


    return listOfSetChainsOfStructuresLinkedBySharedIntegraseToReturn


# add both upstream and downstream valid integrase to a conj module
def addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure(
        listSPs,
        listICEsIMEsStructures,
        maxNumberCDSForFilterIMESize,
        groupListSPintoICEsIMEsUsingFamilyInfo,
        setLocusTagToNotConsiderAseManuallyCheckSent,
        countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure,
        allowAdjacentIntegraseOnlyForSer,
        locusTagIntegrase2Comment,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
        useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance):

    # is False initially
    DEBUG = False


    doNotTakeIntoAccountIfManuallyCheckButAnIntegraseMoreSureIsRegistred = True # allow iterative run of this method in case we found some sure integrase in latter structure that affect the way integrase are handled in upstream structures

    if DEBUG:
        print("\n****** DEBUG addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure :"
            + "\nlistSPs = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listSPs))
            + "\nlistICEsIMEsStructures = {}".format(EMStructure.ICEsIMEsStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures))
            + "\ncountIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure = {}".format(str(countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure))
            )


    countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure += 1
    if countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure > 1000 :
        raise RuntimeError(
            "Error in addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure  countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure > 1000 for listICEsIMEsStructures = {}".format(
                    EMStructure.BasicEMStructure.GetListInternIdFromListEMStructure(listICEsIMEsStructures)))

    #family2SetICEsIMEsStructures = EMStructure.BasicEMStructure.buildFamily2SetICEsIMEsStructures(listICEsIMEsStructures)
    #family2SetPutativeIntegrases = buildamFily2SetPutativeIntegrases(listSPs)

    integraseAttributedSureUpstream2setICEsIMEsStructure = {}
    integraseAttributedSureDownstream2setICEsIMEsStructure = {}
    ICEsIMEsStructure2setIntegraseMoreSure = {}
    IntegraseManuallyCheckUpstream2setICEsIMEsStructure = {}
    IntegraseManuallyCheckDownstream2setICEsIMEsStructure = {}
    ICEsIMEsStructure2IntegraseManuallyCheck2comment = {}

    if DEBUG:
        print("\n ** start for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):")
    for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):
        if DEBUG:
            print(" - Starting currICEsIMEsStructure = {} ;  idxCurrICEsIMEsStructure = {} ; setLocusTagToNotConsiderAseManuallyCheckSent = {} ; countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure = {} ".format(
                    currICEsIMEsStructure.internalIdentifier, str(idxCurrICEsIMEsStructure), repr(setLocusTagToNotConsiderAseManuallyCheckSent), str(countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure)))

        # integrase have already been added with addObviousIntegraseUpstreamAndDownstream_priorMerging
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
            
            if DEBUG:
                print("Continue, integrase have already been added with addObviousIntegraseUpstreamAndDownstream_priorMerging : listIntegraseUpstream = {} and listIntegraseDownstream = {}".format(hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseUpstream), hit.ListSPs.GetListProtIdsFromListSP(currICEsIMEsStructure.listIntegraseDownstream)))
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
            raise RuntimeError("Error in addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure: unrecognized allowAdjacentIntegraseOnlyForSer = {}".format(allowAdjacentIntegraseOnlyForSer))


        #setIntegraseDiscarded = set()
        if currICEsIMEsStructure.delMerging_idxListUpstreamStructure == -1:  # skip if deleted after merge

            listUpstreamIntToAdd = hit.ListSPs.getListIntegraseGroupJustUpstreamOfThisSP(currICEsIMEsStructure.listOrderedSPs[0], listSPs, dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd, setIntegraseTypeToCheck, currICEsIMEsStructure, locusTagIntegrase2Comment)

            listDownstreamIntToAdd = hit.ListSPs.getListIntegraseGroupJustDownstreamOfThisSP(currICEsIMEsStructure.listOrderedSPs[-1], listSPs, dictIntegraseAttributedByCheckForObviousIntegraseUpstreamAndDownstreamToAdd, setIntegraseTypeToCheck, currICEsIMEsStructure, locusTagIntegrase2Comment)

            if DEBUG:
                print(" - 1. getListIntegraseGroupJustUpstreamOfThisSP : listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
                print(" - 1. getListIntegraseGroupJustDownstreamOfThisSP: listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))


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
                print(" - 2. remove integrase that should not be commented (in setLocusTagToNotConsiderAseManuallyCheckSent) : listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
                print(" - 2. remove integrase that should not be commented (in setLocusTagToNotConsiderAseManuallyCheckSent) : listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))

            # if listUpstreamIntToAdd and listDownstreamIntToAdd, check if 1 integrase is way closer than the other
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
                # actually do something with valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase
                if valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase == 1 :
                    listDownstreamIntToAdd.clear()
                elif valueTestUseCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase == 3 :
                    listUpstreamIntToAdd.clear()

                if DEBUG:
                    print(" - 3. both listUpstreamIntToAdd and listDownstreamIntToAdd, check if 1 integrase is way closer than the other (useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase) : listUpstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listUpstreamIntToAdd)))
                    print(" - 3. both listUpstreamIntToAdd and listDownstreamIntToAdd, check if 1 integrase is way closer than the other (useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase) : listDownstreamIntToAdd = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listDownstreamIntToAdd)))


            if listUpstreamIntToAdd and listDownstreamIntToAdd:
                #listLocusTagInComment = []

                setIntegraseBothUpstreamAndDownstream = set()
                for currSp in reversed(listUpstreamIntToAdd):
                    setIntegraseBothUpstreamAndDownstream.add(currSp)
                    # fill up IntegraseManuallyCheckUpstream2setICEsIMEsStructure
                    if currSp in IntegraseManuallyCheckUpstream2setICEsIMEsStructure :
                        setICEsIMEsStructureIT = IntegraseManuallyCheckUpstream2setICEsIMEsStructure[currSp]
                        setICEsIMEsStructureIT.add(currICEsIMEsStructure)
                    else :
                        setICEsIMEsStructureIT = set()
                        setICEsIMEsStructureIT.add(currICEsIMEsStructure)
                        IntegraseManuallyCheckUpstream2setICEsIMEsStructure[currSp] = setICEsIMEsStructureIT
                    listUpstreamIntToAdd.remove(currSp)
                    #listLocusTagInComment.append(currSp.locusTag)
                for currSp in reversed(listDownstreamIntToAdd):
                    setIntegraseBothUpstreamAndDownstream.add(currSp)
                    # fill up IntegraseManuallyCheckDownstream2setICEsIMEsStructure
                    if currSp in IntegraseManuallyCheckDownstream2setICEsIMEsStructure :
                        setICEsIMEsStructureIT = IntegraseManuallyCheckDownstream2setICEsIMEsStructure[currSp]
                        setICEsIMEsStructureIT.add(currICEsIMEsStructure)
                    else :
                        setICEsIMEsStructureIT = set()
                        setICEsIMEsStructureIT.add(currICEsIMEsStructure)
                        IntegraseManuallyCheckDownstream2setICEsIMEsStructure[currSp] = setICEsIMEsStructureIT
                    listDownstreamIntToAdd.remove(currSp)
                    #listLocusTagInComment.append(currSp.locusTag)

                #reprSetIntegraseBothUpstreamAndDownstream = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseBothUpstreamAndDownstream)
                commentITToAdd = "Could not choose assigning between upstream or downstream integrases {} to structure {}, please manually check. ".format(hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseBothUpstreamAndDownstream), currICEsIMEsStructure.internalIdentifier)
                addEntryToICEsIMEsStructure2IntegraseManuallyCheck2comment(currICEsIMEsStructure, setIntegraseBothUpstreamAndDownstream, commentITToAdd, ICEsIMEsStructure2IntegraseManuallyCheck2comment)
                # commented below because can not easily remove comment if the integrase is sure for another structure
                # if commentITToAdd not in currICEsIMEsStructure.comment :
                #     currICEsIMEsStructure.comment += commentITToAdd
                # for integraseManuallyCheckIT in setIntegraseBothUpstreamAndDownstream:
                #     currICEsIMEsStructure.setIntegraseToManuallyCheck.add(integraseManuallyCheckIT.locusTag)
                #     icescreen_OO.addCommentToLocusTag2Comment(
                #         integraseManuallyCheckIT.locusTag,
                #         commentITToAdd,
                #         locusTagIntegrase2Comment)



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

    # end of 1st loop -> for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):

    #OLD WAS HERE BEFORE, got moved after
    # if doNotTakeIntoAccountIfManuallyCheckButAnIntegraseMoreSureIsRegistred:
    #     ...


    if DEBUG:
        strToPrintintegraseAttributedSureUpstream2setICEsIMEsStructure = "{"
        for keyIntegraseAttributedSureUpstream, valueSetICEsIMEsStructure in integraseAttributedSureUpstream2setICEsIMEsStructure.items():
            strToPrintintegraseAttributedSureUpstream2setICEsIMEsStructure += keyIntegraseAttributedSureUpstream.locusTag + " : " + EMStructure.ICEsIMEsStructure.GetListInternIdFromSetEMStructure(valueSetICEsIMEsStructure) + " / "
        strToPrintintegraseAttributedSureUpstream2setICEsIMEsStructure = "}"
        strToPrintintegraseAttributedSureDownstream2setICEsIMEsStructure = "{"
        for keyIntegraseAttributedSureDownstream, valueSetICEsIMEsStructure in integraseAttributedSureDownstream2setICEsIMEsStructure.items():
            strToPrintintegraseAttributedSureDownstream2setICEsIMEsStructure += keyIntegraseAttributedSureDownstream.locusTag + " : " + EMStructure.ICEsIMEsStructure.GetListInternIdFromSetEMStructure(valueSetICEsIMEsStructure) + " / "
        strToPrintintegraseAttributedSureDownstream2setICEsIMEsStructure = "}"
        print("\n ** Done for idxCurrICEsIMEsStructure, currICEsIMEsStructure in enumerate(listICEsIMEsStructures):\n\tintegraseAttributedSureDownstream2setICEsIMEsStructure = {}\n\tintegraseAttributedSureUpstream2setICEsIMEsStructure = {}".format(strToPrintintegraseAttributedSureUpstream2setICEsIMEsStructure, strToPrintintegraseAttributedSureDownstream2setICEsIMEsStructure))




    # integrase more sure after phase 1 are in integraseAttributedSureUpstream2setICEsIMEsStructure and ICEsIMEsStructure2setIntegraseMoreSure
    
    afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister = {}
    for keyIntegraseAttributedSureUpstream, valueSetICEsIMEsStructure in integraseAttributedSureUpstream2setICEsIMEsStructure.items():
        currSetICEsIMEsStructure = set()
        currSetICEsIMEsStructure.update(valueSetICEsIMEsStructure)
        if keyIntegraseAttributedSureUpstream in integraseAttributedSureDownstream2setICEsIMEsStructure:
            currSetICEsIMEsStructure.update(integraseAttributedSureDownstream2setICEsIMEsStructure[keyIntegraseAttributedSureUpstream])
            
        if len(currSetICEsIMEsStructure) >= 2:

            if DEBUG:
                print("UpstreamToRegister len(currSetICEsIMEsStructure) >= 2 for integrase {} which is shared between adjacent structures {}".format(keyIntegraseAttributedSureUpstream.locusTag, EMStructure.ICEsIMEsStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                )

            # move to To check
            del integraseAttributedSureDownstream2setICEsIMEsStructure[keyIntegraseAttributedSureUpstream]
            integraseCouldBeAttributedToMultipleICEIMEStructures(keyIntegraseAttributedSureUpstream, currSetICEsIMEsStructure, False, locusTagIntegrase2Comment) #, ICEsIMEsStructure2IntegraseManuallyCheck2comment

        else:
            # register integrase to be added to the ICEsIMEsStructure
            for currValueICEsIMEsStructure in currSetICEsIMEsStructure:
                if currValueICEsIMEsStructure in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister:
                    # raise RuntimeError(
                    #     "Error in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister: want to register value keyIntegraseAttributedSureUpstream = {} as key currValueICEsIMEsStructure = {} but this key is already in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister with value {}".format(
                    #         keyIntegraseAttributedSureUpstream.locusTag,
                    #         currValueICEsIMEsStructure.internalIdentifier,
                    #         afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[currValueICEsIMEsStructure].locusTag
                    #         ))
                    setIntegraseAttributedSureUpstream = afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[currValueICEsIMEsStructure]
                    setIntegraseAttributedSureUpstream.add(keyIntegraseAttributedSureUpstream)
                else :
                    setIntegraseAttributedSureUpstream = set()
                    setIntegraseAttributedSureUpstream.add(keyIntegraseAttributedSureUpstream)
                    afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[currValueICEsIMEsStructure] = setIntegraseAttributedSureUpstream
                # currValueICEsIMEsStructure.listIntegraseUpstream.append(keyIntegraseAttributedSureUpstream)
                # currValueICEsIMEsStructure.listOrderedSPs.append(keyIntegraseAttributedSureUpstream)
                # currValueICEsIMEsStructure.refreshListIdxOrderedSPs()
                break

    afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister = {}
    for keyIntegraseAttributedSureDownstream, valueSetICEsIMEsStructure in integraseAttributedSureDownstream2setICEsIMEsStructure.items():
        currSetICEsIMEsStructure = set()
        currSetICEsIMEsStructure.update(valueSetICEsIMEsStructure)
        if keyIntegraseAttributedSureDownstream in integraseAttributedSureUpstream2setICEsIMEsStructure:
            currSetICEsIMEsStructure.update(integraseAttributedSureUpstream2setICEsIMEsStructure[keyIntegraseAttributedSureDownstream])
            
        if len(currSetICEsIMEsStructure) >= 2:
            if DEBUG:
                print("DownstreamToRegister len(currSetICEsIMEsStructure) >= 2 for integrase {} which is shared between adjacent structures {}".format(keyIntegraseAttributedSureDownstream.locusTag, EMStructure.ICEsIMEsStructure.GetListInternIdFromSetEMStructure(currSetICEsIMEsStructure))
                )

            # move to To check
            del integraseAttributedSureUpstream2setICEsIMEsStructure[keyIntegraseAttributedSureDownstream]
            integraseCouldBeAttributedToMultipleICEIMEStructures(keyIntegraseAttributedSureDownstream, currSetICEsIMEsStructure, False, locusTagIntegrase2Comment) #, ICEsIMEsStructure2IntegraseManuallyCheck2comment

        else:
            # register integrase to be added to the ICEsIMEsStructure
            for currValueICEsIMEsStructure in currSetICEsIMEsStructure:
                if currValueICEsIMEsStructure in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister:
                    # raise RuntimeError(
                    #     # "Error in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister: currValueICEsIMEsStructure = {} already in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister".format(
                    #     #     currValueICEsIMEsStructure.internalIdentifier))
                    #     "Error in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister: want to register value keyIntegraseAttributedSureDownstream = {} as key currValueICEsIMEsStructure = {} but this key is already in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister with value {}".format(
                    #         keyIntegraseAttributedSureDownstream.locusTag,
                    #         currValueICEsIMEsStructure.internalIdentifier,
                    #         afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[currValueICEsIMEsStructure].locusTag
                    #         ))
                    setIntegraseAttributedSureDownstream = afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[currValueICEsIMEsStructure]
                    setIntegraseAttributedSureDownstream.add(keyIntegraseAttributedSureDownstream)
                else :
                    setIntegraseAttributedSureDownstream = set()
                    setIntegraseAttributedSureDownstream.add(keyIntegraseAttributedSureDownstream)
                    afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[currValueICEsIMEsStructure] = setIntegraseAttributedSureDownstream
                # currValueICEsIMEsStructure.listIntegraseDownstream.append(keyIntegraseAttributedSureDownstream)
                # currValueICEsIMEsStructure.listOrderedSPs.append(keyIntegraseAttributedSureDownstream)
                # currValueICEsIMEsStructure.refreshListIdxOrderedSPs()
                break



    if DEBUG:
        print("\n** starting listOfSetChainsOfStructuresLinkedBySharedIntegrase")

    # create a chain of structures that are linked by integrase in between them. Basis for latter test to check if upstream nd downstream chain reaction of changes regarding integrase assignment
    listOfSetChainsOfStructuresLinkedBySharedIntegrase = createListOfListChainsOfStructuresLinkedBySharedIntegrase(integraseAttributedSureUpstream2setICEsIMEsStructure, integraseAttributedSureDownstream2setICEsIMEsStructure, IntegraseManuallyCheckUpstream2setICEsIMEsStructure, IntegraseManuallyCheckDownstream2setICEsIMEsStructure, DEBUG)


    # test a chain of structures that are linked by possibly shared integrase in between them : check if upstream and downstream chain reaction of changes regarding integrase assignment
    for setChainOfStructuresLinkedBySharedIntegrase in listOfSetChainsOfStructuresLinkedBySharedIntegrase :
        if DEBUG:
                print(" - looping through listOfSetChainsOfStructuresLinkedBySharedIntegrase : setChainOfStructuresLinkedBySharedIntegrase = {}".format(EMStructure.ICEsIMEsStructure.GetListInternIdFromSetEMStructure(setChainOfStructuresLinkedBySharedIntegrase)))
        upstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase = False
        downstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase = False
        for structureToCheckInChainIT in setChainOfStructuresLinkedBySharedIntegrase :
            if structureToCheckInChainIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister :
                upstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase = True
            if structureToCheckInChainIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister :
                downstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase = True
        if upstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase and downstreamSureIntegraseWillAffectThisChainOfStructuresLinkedBySharedIntegrase :
            # This chain of structure will be affected by both upstream and downstream ure integrase, it is not safe to assume sure integrase anymore for it
            if DEBUG:
                print("This chain of structure will be affected by both upstream and downstream ure integrase, it is not safe to assume sure integrase anymore for it")
            setIntegraseAffectedByContradictoryChainReaction = set()
            for structureToCheckInChainIT in setChainOfStructuresLinkedBySharedIntegrase :
                if structureToCheckInChainIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister :
                    setIntegraseAffectedByContradictoryChainReaction.update(afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[structureToCheckInChainIT])
                    del afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[structureToCheckInChainIT]
                if structureToCheckInChainIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister :
                    setIntegraseAffectedByContradictoryChainReaction.update(afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[structureToCheckInChainIT])
                    del afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[structureToCheckInChainIT]
                if structureToCheckInChainIT in ICEsIMEsStructure2setIntegraseMoreSure :
                    setIntegraseAffectedByContradictoryChainReaction.update(ICEsIMEsStructure2setIntegraseMoreSure[structureToCheckInChainIT])
                    #del after
            # add comment to structureToCheckInChainIT
            commentITToAdd_contradictoryChainReaction = "Integrases {} are not automatically assigned to their surrounding structures of conjugation module as some contradictory upstream and downstream integrase assignments were detected. ".format(hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseAffectedByContradictoryChainReaction))
            for structureToCheckInChainIT in setChainOfStructuresLinkedBySharedIntegrase :
                if structureToCheckInChainIT in ICEsIMEsStructure2setIntegraseMoreSure :
                    # add the integrase as not sure if they were registred as sure
                    addEntryToICEsIMEsStructure2IntegraseManuallyCheck2comment(structureToCheckInChainIT, ICEsIMEsStructure2setIntegraseMoreSure[structureToCheckInChainIT], commentITToAdd_contradictoryChainReaction, ICEsIMEsStructure2IntegraseManuallyCheck2comment)
                    # if commentITToAdd_contradictoryChainReaction not in structureToCheckInChainIT.comment :
                    #     structureToCheckInChainIT.comment += commentITToAdd_contradictoryChainReaction
                    # for integraseManuallyCheckIT in ICEsIMEsStructure2setIntegraseMoreSure[structureToCheckInChainIT]:
                    #     structureToCheckInChainIT.setIntegraseToManuallyCheck.add(integraseManuallyCheckIT.locusTag)
                    #     icescreen_OO.addCommentToLocusTag2Comment(
                    #         integraseManuallyCheckIT.locusTag,
                    #         commentITToAdd_contradictoryChainReaction,
                    #         locusTagIntegrase2Comment)
                    del ICEsIMEsStructure2setIntegraseMoreSure[structureToCheckInChainIT]



    # check for possible nested element between the conj module and the inegrase. if 2 module conj seems to be in accretion and 1 do not have an integrase, we can not assign the integrase for sure to the closest structure as it could be nested between the inegrase and the conj module of the structure next to it 
    
    #print("\n****** check for possible nested element between the conj module and the inegrase : upstream integrase ; size = {}".format(str(len(afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.keys()))))
    setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister = set()
    for keyICEsIMEsStructureIT, valueSetIntegraseUpstreamToRegisterIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.items() :
        #print("aICEsIMEsStructure2IntegraseUpstreamToRegister: keyICEsIMEsStructureIT={} ({}) -> {}".format(keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(keyICEsIMEsStructureIT.listOrderedSPs), hit.ListSPs.GetListProtIdsFromSetSP(valueSetIntegraseUpstreamToRegisterIT)))
        #get the ICEsIMEsStructureJustDownstreamIT to our currICEsIMEsStructure
        #ICEsIMEsStructureJustDownstreamIT will be None if deleted after merge
        ICEsIMEsStructureJustDownstreamIT = EMStructure.ICEsIMEsStructure.getICEsIMEsStructureDownstreamOfICEsIMEsStructure(keyICEsIMEsStructureIT, listICEsIMEsStructures, True)
        #print("aICEsIMEsStructure2IntegraseUpstreamToRegister: ICEsIMEsStructureJustDownstreamIT={} ({})".format(ICEsIMEsStructureJustDownstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureJustDownstreamIT.listOrderedSPs)))
        if ICEsIMEsStructureJustDownstreamIT is not None :
            #print("ICEsIMEsStructureJustDownstreamIT is not None")
            #change to can be MasterMergeStructure (False), check if encompass our keyICEsIMEsStructureIT (has SP upstream and downstream of it), if so then fail test just as if ICEsIMEsStructureJustDownstreamIT was None
            mostUpstrSpNotIntegr_keyICEsIMEsStructureIT = keyICEsIMEsStructureIT.findMostUpstrSpNotIntegr()
            mostDownstrSpNotIntegr_keyICEsIMEsStructureIT = keyICEsIMEsStructureIT.findMostDownstrSpNotIntegr()
            mostUpstrSpNotIntegr_ICEsIMEsStructureJustDownstreamIT = ICEsIMEsStructureJustDownstreamIT.findMostUpstrSpNotIntegr()
            mostDownstrSpNotIntegr_ICEsIMEsStructureJustDownstreamIT = ICEsIMEsStructureJustDownstreamIT.findMostDownstrSpNotIntegr()
            keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustDownstreamIT = False
            if mostUpstrSpNotIntegr_keyICEsIMEsStructureIT.CDSPositionInGenome > mostUpstrSpNotIntegr_ICEsIMEsStructureJustDownstreamIT.CDSPositionInGenome and mostDownstrSpNotIntegr_keyICEsIMEsStructureIT.CDSPositionInGenome < mostDownstrSpNotIntegr_ICEsIMEsStructureJustDownstreamIT.CDSPositionInGenome :
                keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustDownstreamIT = True
            if not keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustDownstreamIT :
                #ICEsIMEsStructureJustDownstreamIT could have some integrase assign pre-merging
                if len(ICEsIMEsStructureJustDownstreamIT.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureJustDownstreamIT.listIntegraseDownstream) == 0:
                    if ICEsIMEsStructureJustDownstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister :
                        #print("ICEsIMEsStructureJustDownstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister")
                        if ICEsIMEsStructureJustDownstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister:
                            #print("ICEsIMEsStructureJustDownstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister")

                            if not EMStructure.ICEsIMEsStructure.IsThereAnIntegraseBetweenThoseTwoConjModule(keyICEsIMEsStructureIT, ICEsIMEsStructureJustDownstreamIT, listSPs, setIntegraseTypeToCheck) :

                                convertListIntegraseUpstreamToRegisterIT = list(valueSetIntegraseUpstreamToRegisterIT)
                                convertListIntegraseUpstreamToRegisterIT.sort(key=lambda x: (x.genomeAccessionRank, x.start), reverse=False) #key=lambda x: x.start

                                if isIntegraseCorrectlyOrientedForICEsIMEsStructure(convertListIntegraseUpstreamToRegisterIT[0], ICEsIMEsStructureJustDownstreamIT, True) == 2 :

                                    # test the possibility of another integrase that follow upstream of this one and that are upstream of convertListIntegraseUpstreamToRegisterIT and that is not associated with any other structure

                                    booleanAnotherSingleIntegraseFoundUpstream = False

                                    listIntegraseGroupUpstreamOfSP = hit.ListSPs.getListIntegraseGroupJustUpstreamOfThisSP(convertListIntegraseUpstreamToRegisterIT[0], listSPs, None, setIntegraseTypeToCheck, ICEsIMEsStructureJustDownstreamIT, locusTagIntegrase2Comment)
                                    # #print("TEST listIntegraseGroupUpstreamOfSP = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP)))

                                    if len(listIntegraseGroupUpstreamOfSP) > 0 :
                                        idxStructureUpstreamToCheck = keyICEsIMEsStructureIT.idxInSeedList - 1
                                        booleanHasStructureWithoutRegistredIntergaseUpstreamIT = False
                                        if idxStructureUpstreamToCheck >= 0:
                                            listStucturesUpstreamOfKeyICEsIMEsStructureIT = listICEsIMEsStructures[:(idxStructureUpstreamToCheck + 1)]
                                            if not listStucturesUpstreamOfKeyICEsIMEsStructureIT[-1].structureHasThoseIntegraseRegistred(listIntegraseGroupUpstreamOfSP) :
                                                booleanHasStructureWithoutRegistredIntergaseUpstreamIT = EMStructure.BasicEMStructure.listStucturesHasAtLeastOneStructureWithoutRegistredIntergase(listStucturesUpstreamOfKeyICEsIMEsStructureIT)

                                        if not booleanHasStructureWithoutRegistredIntergaseUpstreamIT :
                                            
                                            ## add the integrase as sure. If want to add the integrase as to check instead, see the commented code below
                                            setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.add(keyICEsIMEsStructureIT)
                                            keyICEsIMEsStructureIT.addListIntegraseUpstream(convertListIntegraseUpstreamToRegisterIT)
                                            ICEsIMEsStructureJustDownstreamIT.addListIntegraseUpstream(listIntegraseGroupUpstreamOfSP)
                                            # add to ICEsIMEsStructure2setIntegraseMoreSure ?
                                            for currIntegraseUpstreamIT in convertListIntegraseUpstreamToRegisterIT:
                                                if keyICEsIMEsStructureIT in ICEsIMEsStructure2setIntegraseMoreSure:
                                                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[keyICEsIMEsStructureIT]
                                                    setIntegraseMoreSure.add(currIntegraseUpstreamIT)
                                                else:
                                                    setIntegraseMoreSure = set()
                                                    setIntegraseMoreSure.add(currIntegraseUpstreamIT)
                                                    ICEsIMEsStructure2setIntegraseMoreSure[keyICEsIMEsStructureIT] = setIntegraseMoreSure
                                            for currIntegraseUpstreamIT in listIntegraseGroupUpstreamOfSP:
                                                if ICEsIMEsStructureJustDownstreamIT in ICEsIMEsStructure2setIntegraseMoreSure:
                                                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[ICEsIMEsStructureJustDownstreamIT]
                                                    setIntegraseMoreSure.add(currIntegraseUpstreamIT)
                                                else:
                                                    setIntegraseMoreSure = set()
                                                    setIntegraseMoreSure.add(currIntegraseUpstreamIT)
                                                    ICEsIMEsStructure2setIntegraseMoreSure[ICEsIMEsStructureJustDownstreamIT] = setIntegraseMoreSure
                                            # add comments
                                            commentITToAdd_keyICEsIMEsStructureIT = "Upstream integrase {} can be associated to structure {} as nested within structure {} associated to integrase {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseUpstreamToRegisterIT), keyICEsIMEsStructureIT.internalIdentifier, ICEsIMEsStructureJustDownstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP))
                                            if commentITToAdd_keyICEsIMEsStructureIT not in keyICEsIMEsStructureIT.comment :
                                                keyICEsIMEsStructureIT.comment += commentITToAdd_keyICEsIMEsStructureIT
                                            for currIntegraseUpstreamIT in convertListIntegraseUpstreamToRegisterIT:
                                                icescreen_OO.addCommentToLocusTag2Comment(
                                                    currIntegraseUpstreamIT.locusTag,
                                                    commentITToAdd_keyICEsIMEsStructureIT,
                                                    locusTagIntegrase2Comment)
                                            commentITToAdd_ICEsIMEsStructureJustDownstreamIT = "Upstream integrase {} can be associated to structure {} and form a host structure of structure {} associated to integrase {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP), ICEsIMEsStructureJustDownstreamIT.internalIdentifier, keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseUpstreamToRegisterIT))
                                            if commentITToAdd_ICEsIMEsStructureJustDownstreamIT not in ICEsIMEsStructureJustDownstreamIT.comment :
                                                ICEsIMEsStructureJustDownstreamIT.comment += commentITToAdd_ICEsIMEsStructureJustDownstreamIT
                                            for currIntegraseUpstreamIT in listIntegraseGroupUpstreamOfSP:
                                                icescreen_OO.addCommentToLocusTag2Comment(
                                                    currIntegraseUpstreamIT.locusTag,
                                                    commentITToAdd_ICEsIMEsStructureJustDownstreamIT,
                                                    locusTagIntegrase2Comment)
                                            ## If I want to add the integrase as "to manually check" instead, uncomment bolow and comment up
                                            # setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.add(keyICEsIMEsStructureIT)
                                            # add comments
                                            # setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT = {}
                                            # commentITToAdd_keyICEsIMEsStructureIT = "Upstream integrase {} could be associated to structure {} as nested within structure {} associated to integrase {}, please check. ".format(hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseUpstreamToRegisterIT), keyICEsIMEsStructureIT.internalIdentifier, ICEsIMEsStructureJustDownstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP))
                                            # #reprSetIntegraseBothDownstreamAndUpstream = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseBothDownstreamAndUpstream)
                                            # setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT[hit.ListSPs.GetListProtIdsFromSetSP(valueSetIntegraseUpstreamToRegisterIT)] = commentITToAdd_keyICEsIMEsStructureIT
                                            # ICEsIMEsStructure2IntegraseManuallyCheck2comment[keyICEsIMEsStructureIT] = setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT
                                            # setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustDownstreamIT = {}
                                            # commentITToAdd_ICEsIMEsStructureJustDownstreamIT = "Upstream integrase {} could be associated to structure {} and form a host structure of structure {} associated to integrase {}, please check. ".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP), ICEsIMEsStructureJustDownstreamIT.internalIdentifier, keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseUpstreamToRegisterIT))
                                            # setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustDownstreamIT[hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupUpstreamOfSP)] = commentITToAdd_ICEsIMEsStructureJustDownstreamIT
                                            # ICEsIMEsStructure2IntegraseManuallyCheck2comment[ICEsIMEsStructureJustDownstreamIT] = setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustDownstreamIT

                                            booleanAnotherSingleIntegraseFoundUpstream = True

                                    if not booleanAnotherSingleIntegraseFoundUpstream :
                                        
                                        for valueIntegraseUpstreamToRegisterIT in convertListIntegraseUpstreamToRegisterIT:

                                            #print("valueIntegraseUpstreamToRegisterIT={}".format(valueIntegraseUpstreamToRegisterIT.locusTag))

                                            #print("isIntegraseCorrectlyOrientedForICEsIMEsStructure(valueIntegraseUpstreamToRegisterIT, ICEsIMEsStructureJustDownstreamIT, True) == 2")
                                            # this integrase could also be attributed to ICEsIMEsStructureJustDownstreamIT in case of nested relationship
                                            currSetICEsIMEsStructure = set()
                                            currSetICEsIMEsStructure.add(keyICEsIMEsStructureIT)
                                            currSetICEsIMEsStructure.add(ICEsIMEsStructureJustDownstreamIT)
                                            integraseCouldBeAttributedToMultipleICEIMEStructures(valueIntegraseUpstreamToRegisterIT, currSetICEsIMEsStructure, True, locusTagIntegrase2Comment) #, ICEsIMEsStructure2IntegraseManuallyCheck2comment
                                            setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.add(keyICEsIMEsStructureIT)

    #print("\n****** check for possible nested element between the conj module and the inegrase : downstream integrase")
    setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister = set()
    for keyICEsIMEsStructureIT, valueSetIntegraseDownstreamToRegisterIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister.items() :
        #print("aICEsIMEsStructure2IntegraseDownstreamToRegister: keyICEsIMEsStructureIT={} ({}) -> {}".format(keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(keyICEsIMEsStructureIT.listOrderedSPs), hit.ListSPs.GetListProtIdsFromSetSP(valueSetIntegraseDownstreamToRegisterIT)))
        #get the ICEsIMEsStructureJustUpstreamIT to our currICEsIMEsStructure
        ICEsIMEsStructureJustUpstreamIT = EMStructure.ICEsIMEsStructure.getICEsIMEsStructureUpstreamOfICEsIMEsStructure(keyICEsIMEsStructureIT, listICEsIMEsStructures, False)
        #print("aICEsIMEsStructure2IntegraseUpstreamToRegister: ICEsIMEsStructureJustUpstreamIT={} ({})".format(ICEsIMEsStructureJustUpstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureJustUpstreamIT.listOrderedSPs)))
        if ICEsIMEsStructureJustUpstreamIT is not None :
            #print("ICEsIMEsStructureJustUpstreamIT is not None")
            #ICEsIMEsStructureJustUpstreamIT can be MasterMergeStructure, check if encompass our keyICEsIMEsStructureIT (has SP upstream and downstream of it), if so then fail test just as if ICEsIMEsStructureJustUpstreamIT was None
            mostUpstrSpNotIntegr_keyICEsIMEsStructureIT = keyICEsIMEsStructureIT.findMostUpstrSpNotIntegr()
            mostDownstrSpNotIntegr_keyICEsIMEsStructureIT = keyICEsIMEsStructureIT.findMostDownstrSpNotIntegr()
            mostUpstrSpNotIntegr_ICEsIMEsStructureJustUpstreamIT = ICEsIMEsStructureJustUpstreamIT.findMostUpstrSpNotIntegr()
            mostDownstrSpNotIntegr_ICEsIMEsStructureJustUpstreamIT = ICEsIMEsStructureJustUpstreamIT.findMostDownstrSpNotIntegr()
            keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustUpstreamIT = False
            if mostUpstrSpNotIntegr_keyICEsIMEsStructureIT.CDSPositionInGenome > mostUpstrSpNotIntegr_ICEsIMEsStructureJustUpstreamIT.CDSPositionInGenome and mostDownstrSpNotIntegr_keyICEsIMEsStructureIT.CDSPositionInGenome < mostDownstrSpNotIntegr_ICEsIMEsStructureJustUpstreamIT.CDSPositionInGenome :
                keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustUpstreamIT = True
            if not keyICEsIMEsStructureITIsNestedInICEsIMEsStructureJustUpstreamIT :
                #ICEsIMEsStructureJustDownstreamIT could have some integrase assign pre-merging
                if len(ICEsIMEsStructureJustUpstreamIT.listIntegraseUpstream) == 0 and len(ICEsIMEsStructureJustUpstreamIT.listIntegraseDownstream) == 0:
                    if ICEsIMEsStructureJustUpstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister :
                        #print("ICEsIMEsStructureJustUpstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister")
                        if ICEsIMEsStructureJustUpstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister :
                            #print("ICEsIMEsStructureJustUpstreamIT not in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister")

                            if not EMStructure.ICEsIMEsStructure.IsThereAnIntegraseBetweenThoseTwoConjModule(keyICEsIMEsStructureIT, ICEsIMEsStructureJustUpstreamIT, listSPs, setIntegraseTypeToCheck) :

                                convertListIntegraseDownstreamToRegisterIT = list(valueSetIntegraseDownstreamToRegisterIT)
                                convertListIntegraseDownstreamToRegisterIT.sort(key=lambda x: (x.genomeAccessionRank, x.start), reverse=False) #key=lambda x: x.start

                                if isIntegraseCorrectlyOrientedForICEsIMEsStructure(convertListIntegraseDownstreamToRegisterIT[0], ICEsIMEsStructureJustUpstreamIT, False) == 2 :

                                    # test the possibility of another integrase (or integrases that follow up each others -> make a method for that by re-using code already made) that are downstream of convertListIntegraseDownstreamToRegisterIT and that is not associated with any other structure
                                    listIntegraseGroupDownstreamOfSP = hit.ListSPs.getListIntegraseGroupJustDownstreamOfThisSP(convertListIntegraseDownstreamToRegisterIT[-1], listSPs, None, setIntegraseTypeToCheck, ICEsIMEsStructureJustUpstreamIT, locusTagIntegrase2Comment)
                                    #print("TEST listIntegraseGroupDownstreamOfSP = {}".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP)))

                                    booleanAnotherSingleIntegraseFoundDownstream = False
                                    if len(listIntegraseGroupDownstreamOfSP) > 0 :
                                        idxStructureDownstreamToCheck = keyICEsIMEsStructureIT.idxInSeedList + 1
                                        booleanHasStructureWithoutRegistredIntergaseDownstreamIT = False
                                        if idxStructureDownstreamToCheck < len(listICEsIMEsStructures):
                                            listStucturesDownstreamOfKeyICEsIMEsStructureIT = listICEsIMEsStructures[idxStructureDownstreamToCheck:]
                                            if not listStucturesDownstreamOfKeyICEsIMEsStructureIT[0].structureHasThoseIntegraseRegistred(listIntegraseGroupDownstreamOfSP) :
                                                booleanHasStructureWithoutRegistredIntergaseDownstreamIT = EMStructure.BasicEMStructure.listStucturesHasAtLeastOneStructureWithoutRegistredIntergase(listStucturesDownstreamOfKeyICEsIMEsStructureIT)

                                        if not booleanHasStructureWithoutRegistredIntergaseDownstreamIT :
                                            
                                            ## add the integrase as sure. If want to add the integrase as to check instead, see the commented code below
                                            setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister.add(keyICEsIMEsStructureIT)
                                            keyICEsIMEsStructureIT.addListIntegraseDownstream(convertListIntegraseDownstreamToRegisterIT)
                                            ICEsIMEsStructureJustUpstreamIT.addListIntegraseDownstream(listIntegraseGroupDownstreamOfSP)
                                            # add to ICEsIMEsStructure2setIntegraseMoreSure ?
                                            for currIntegraseDownstreamIT in convertListIntegraseDownstreamToRegisterIT:
                                                if keyICEsIMEsStructureIT in ICEsIMEsStructure2setIntegraseMoreSure:
                                                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[keyICEsIMEsStructureIT]
                                                    setIntegraseMoreSure.add(currIntegraseDownstreamIT)
                                                else:
                                                    setIntegraseMoreSure = set()
                                                    setIntegraseMoreSure.add(currIntegraseDownstreamIT)
                                                    ICEsIMEsStructure2setIntegraseMoreSure[keyICEsIMEsStructureIT] = setIntegraseMoreSure
                                            for currIntegraseDownstreamIT in listIntegraseGroupDownstreamOfSP:
                                                if ICEsIMEsStructureJustUpstreamIT in ICEsIMEsStructure2setIntegraseMoreSure:
                                                    setIntegraseMoreSure = ICEsIMEsStructure2setIntegraseMoreSure[ICEsIMEsStructureJustUpstreamIT]
                                                    setIntegraseMoreSure.add(currIntegraseDownstreamIT)
                                                else:
                                                    setIntegraseMoreSure = set()
                                                    setIntegraseMoreSure.add(currIntegraseDownstreamIT)
                                                    ICEsIMEsStructure2setIntegraseMoreSure[ICEsIMEsStructureJustUpstreamIT] = setIntegraseMoreSure
                                            # add comments
                                            commentITToAdd_keyICEsIMEsStructureIT = "Downstream integrase {} can be associated to structure {} as nested within structure {} associated to integrase {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseDownstreamToRegisterIT), keyICEsIMEsStructureIT.internalIdentifier, ICEsIMEsStructureJustUpstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP))
                                            keyICEsIMEsStructureIT.comment += commentITToAdd_keyICEsIMEsStructureIT
                                            for currIntegraseDownstreamIT in convertListIntegraseDownstreamToRegisterIT:
                                                icescreen_OO.addCommentToLocusTag2Comment(
                                                    currIntegraseDownstreamIT.locusTag,
                                                    commentITToAdd_keyICEsIMEsStructureIT,
                                                    locusTagIntegrase2Comment)
                                            commentITToAdd_ICEsIMEsStructureJustUpstreamIT = "Downstream integrase {} can be associated to structure {} and form a host structure of structure {} associated to integrase {}. ".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP), ICEsIMEsStructureJustUpstreamIT.internalIdentifier, keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseDownstreamToRegisterIT))
                                            ICEsIMEsStructureJustUpstreamIT.comment += commentITToAdd_ICEsIMEsStructureJustUpstreamIT
                                            for currIntegraseDownstreamIT in listIntegraseGroupDownstreamOfSP:
                                                icescreen_OO.addCommentToLocusTag2Comment(
                                                    currIntegraseDownstreamIT.locusTag,
                                                    commentITToAdd_ICEsIMEsStructureJustUpstreamIT,
                                                    locusTagIntegrase2Comment)
                                            ## If I want to add the integrase as "to manually check" instead, uncomment bolow and comment up
                                            # setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister.add(keyICEsIMEsStructureIT)
                                            # add comments
                                            # setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT = {}
                                            # commentITToAdd_keyICEsIMEsStructureIT = "Downstream integrase {} could be associated to structure {} as nested within structure {} associated to integrase {}, please check. ".format(hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseDownstreamToRegisterIT), keyICEsIMEsStructureIT.internalIdentifier, ICEsIMEsStructureJustUpstreamIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP))
                                            # #reprSetIntegraseBothUpstreamAndDownstream = hit.ListSPs.GetListProtIdsFromSetSP(setIntegraseBothUpstreamAndDownstream)
                                            # setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT[hit.ListSPs.GetListProtIdsFromSetSP(valueSetIntegraseDownstreamToRegisterIT)] = commentITToAdd_keyICEsIMEsStructureIT
                                            # ICEsIMEsStructure2IntegraseManuallyCheck2comment[keyICEsIMEsStructureIT] = setIntegraseManuallyCheck2comment_keyICEsIMEsStructureIT
                                            # setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustUpstreamIT = {}
                                            # commentITToAdd_ICEsIMEsStructureJustUpstreamIT = "Downstream integrase {} could be associated to structure {} and form a host structure of structure {} associated to integrase {}, please check. ".format(hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP), ICEsIMEsStructureJustUpstreamIT.internalIdentifier, keyICEsIMEsStructureIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(convertListIntegraseDownstreamToRegisterIT))
                                            # setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustUpstreamIT[hit.ListSPs.GetListProtIdsFromListSP(listIntegraseGroupDownstreamOfSP)] = commentITToAdd_ICEsIMEsStructureJustUpstreamIT
                                            # ICEsIMEsStructure2IntegraseManuallyCheck2comment[ICEsIMEsStructureJustUpstreamIT] = setIntegraseManuallyCheck2comment_ICEsIMEsStructureJustUpstreamIT

                                            booleanAnotherSingleIntegraseFoundDownstream = True

                                    if not booleanAnotherSingleIntegraseFoundDownstream :

                                        for valueIntegraseDownstreamToRegisterIT in convertListIntegraseDownstreamToRegisterIT:
                                            #print("valueIntegraseDownstreamToRegisterIT={}".format(valueIntegraseDownstreamToRegisterIT.locusTag))
                                            
                                            #print("isIntegraseCorrectlyOrientedForICEsIMEsStructure(valueIntegraseDownstreamToRegisterIT, ICEsIMEsStructureJustUpstreamIT, False) == 2")
                                            # this integrase could also be attributed to ICEsIMEsStructureJustUpstreamIT in case of nested relationship
                                            currSetICEsIMEsStructure = set()
                                            currSetICEsIMEsStructure.add(keyICEsIMEsStructureIT)
                                            currSetICEsIMEsStructure.add(ICEsIMEsStructureJustUpstreamIT)
                                            integraseCouldBeAttributedToMultipleICEIMEStructures(valueIntegraseDownstreamToRegisterIT, currSetICEsIMEsStructure, True, locusTagIntegrase2Comment) #, ICEsIMEsStructure2IntegraseManuallyCheck2comment
                                            setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister.add(keyICEsIMEsStructureIT)

    # delete entries with setICEsIMEsStructureToRemoveFrom
    for ICEsIMEsStructureToRemoveIt in setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister:
        del afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister[ICEsIMEsStructureToRemoveIt]
    for ICEsIMEsStructureToRemoveIt in setICEsIMEsStructureToRemoveFromafterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister:
        del afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister[ICEsIMEsStructureToRemoveIt]



    # register integrase to ICEsIMEsStructure
    for keyICEsIMEsStructureIT, valueSetIntegraseUpstreamToRegisterIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseUpstreamToRegister.items() :
        convertList_valueSetIntegraseUpstreamToRegisterIT = list(valueSetIntegraseUpstreamToRegisterIT)
        keyICEsIMEsStructureIT.addListIntegraseUpstream(convertList_valueSetIntegraseUpstreamToRegisterIT)

    for keyICEsIMEsStructureIT, valueSetIntegraseDownstreamToRegisterIT in afterMultipleICEIMEStructuresCheck_ICEsIMEsStructure2IntegraseDownstreamToRegister.items():
        convertList_valueSetIntegraseDownstreamToRegisterIT = list(valueSetIntegraseDownstreamToRegisterIT)
        keyICEsIMEsStructureIT.addListIntegraseDownstream(convertList_valueSetIntegraseDownstreamToRegisterIT)


    # was moved here
    if doNotTakeIntoAccountIfManuallyCheckButAnIntegraseMoreSureIsRegistred:
        setLocusTagToNotConsiderAseManuallyCheck = set()
        for keyICEsIMEsStructureManuallyCheck, integraseManuallyCheck2comment in ICEsIMEsStructure2IntegraseManuallyCheck2comment.items():
            for integraseManuallyCheck, valueComment in integraseManuallyCheck2comment.items():
                # keyICEsIMEsStructureManuallyCheck.comment += valueComment
                # listIntegraseLocusTag = keyReprSetIntegraseManuallyCheck.split(", ")
                # for currIntegraseLocusTagManuallyCheck in listIntegraseLocusTag:
                for keyICEsIMEsStructureMoreSure, valueSetIntegraseMoreSure in ICEsIMEsStructure2setIntegraseMoreSure.items():
                    for currIntegraseMoreSure in valueSetIntegraseMoreSure:
                        if currIntegraseMoreSure == integraseManuallyCheck:
                            # print("keyReprSetIntegraseManuallyCheck {}: collusion locus tag {} and keyICEsIMEsStructureManuallyCheck {} - keyICEsIMEsStructureMoreSure {}"\
                            #      .format(repr(keyReprSetIntegraseManuallyCheck), currIntegraseMoreSure.locusTag, keyICEsIMEsStructureManuallyCheck.internalIdentifier, keyICEsIMEsStructureMoreSure.internalIdentifier))
                            setLocusTagToNotConsiderAseManuallyCheck.add(integraseManuallyCheck.locusTag)
                            break
        setUnionLocusTagToNotConsiderAseManuallyCheck = setLocusTagToNotConsiderAseManuallyCheckSent.union(setLocusTagToNotConsiderAseManuallyCheck)
        if len(setLocusTagToNotConsiderAseManuallyCheck) >= 0 and countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure < 5 and len(setUnionLocusTagToNotConsiderAseManuallyCheck) > len(setLocusTagToNotConsiderAseManuallyCheckSent):
            if DEBUG:
                print("\n !!! rerun addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure with new setUnionLocusTagToNotConsiderAseManuallyCheck = {}\n".format(setUnionLocusTagToNotConsiderAseManuallyCheck))

            addSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure(
                    listSPs,
                    listICEsIMEsStructures,
                    maxNumberCDSForFilterIMESize,
                    groupListSPintoICEsIMEsUsingFamilyInfo,
                    setUnionLocusTagToNotConsiderAseManuallyCheck,
                    countIteraddSPIntegraseUpstreamAndDownstream_afterMergeDistantStructure,
                    allowAdjacentIntegraseOnlyForSer,
                    locusTagIntegrase2Comment,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_lowCutoffCDSDistance,
                    useCDSDistanceToChooseBetweenUpstreamAndDownstreamIntegrase_highCutoffCDSDistance)
            return


    # ICEsIMEsStructure2IntegraseManuallyCheck2comment
    for keyICEsIMEsStructure, integraseManuallyCheck2comment in ICEsIMEsStructure2IntegraseManuallyCheck2comment.items():
        for integraseManuallyCheck, valueComment in integraseManuallyCheck2comment.items():
            if valueComment not in keyICEsIMEsStructure.comment :
                keyICEsIMEsStructure.comment += valueComment
            # listIntegraseLocusTag = keyReprSetIntegraseManuallyCheck.split(", ")
            # for currIntegraseLocusTagManuallyCheck in listIntegraseLocusTag:
            # keyICEsIMEsStructure.setIntegraseToManuallyCheck.add(integraseManuallyCheck.locusTag)
            keyICEsIMEsStructure.setIntegraseToManuallyCheck.add(integraseManuallyCheck)
            icescreen_OO.addCommentToLocusTag2Comment(
                    integraseManuallyCheck.locusTag,
                    valueComment,
                    locusTagIntegrase2Comment)


    # fix if conflict regarding integrase are registred both as sure and unsure. They should all be unsure if so
    for ICEsIMEsStructuresIT in listICEsIMEsStructures :
        if ICEsIMEsStructuresIT.delMerging_idxListUpstreamStructure != -1 :
            continue
        if len(ICEsIMEsStructuresIT.listIntegraseUpstream) > 0 and len(ICEsIMEsStructuresIT.listIntegraseDownstream) > 0 :
            #print("case 1 for {} ({})".format(ICEsIMEsStructuresIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructuresIT.listOrderedSPs)))
            ICEsIMEsStructuresIT.transferAllIntegrasesToManuallyCheck()
        elif len(ICEsIMEsStructuresIT.listIntegraseUpstream) > 0 and len(ICEsIMEsStructuresIT.setIntegraseToManuallyCheck) > 0 :
            #print("case 2 for {} ({})".format(ICEsIMEsStructuresIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructuresIT.listOrderedSPs)))
            ICEsIMEsStructuresIT.transferAllIntegrasesToManuallyCheck()
        elif len(ICEsIMEsStructuresIT.listIntegraseDownstream) > 0 and len(ICEsIMEsStructuresIT.setIntegraseToManuallyCheck) > 0 :
            #print("case 3 for {} ({})".format(ICEsIMEsStructuresIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructuresIT.listOrderedSPs)))
            ICEsIMEsStructuresIT.transferAllIntegrasesToManuallyCheck()
    __integraseAttributedSureUpstream2setICEsIMEsStructures = {}
    __integraseAttributedSureDownstream2setICEsIMEsStructures = {}
    __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures = {}
    for ICEsIMEsStructuresIT in listICEsIMEsStructures :
        if ICEsIMEsStructuresIT.delMerging_idxListUpstreamStructure != -1 :
            continue
        for integraseUpstreamIT in ICEsIMEsStructuresIT.listIntegraseUpstream :
            if integraseUpstreamIT in __integraseAttributedSureUpstream2setICEsIMEsStructures:
                setICEsIMEsStructuresIT = __integraseAttributedSureUpstream2setICEsIMEsStructures[integraseUpstreamIT]
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
            else :
                setICEsIMEsStructuresIT = set()
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
                __integraseAttributedSureUpstream2setICEsIMEsStructures[integraseUpstreamIT] = setICEsIMEsStructuresIT
        for integraseDownstreamIT in ICEsIMEsStructuresIT.listIntegraseDownstream :
            if integraseDownstreamIT in __integraseAttributedSureDownstream2setICEsIMEsStructures:
                setICEsIMEsStructuresIT = __integraseAttributedSureDownstream2setICEsIMEsStructures[integraseDownstreamIT]
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
            else :
                setICEsIMEsStructuresIT = set()
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
                __integraseAttributedSureDownstream2setICEsIMEsStructures[integraseDownstreamIT] = setICEsIMEsStructuresIT

        # for integraseLocusTagsToManuallyCheckIT in ICEsIMEsStructuresIT.setIntegraseToManuallyCheck :
        #     if integraseLocusTagsToManuallyCheckIT in __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures:
        #         setICEsIMEsStructuresIT = __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures[integraseLocusTagsToManuallyCheckIT]
        #         setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
        #     else :
        #         setICEsIMEsStructuresIT = set()
        #         setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
        #         __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures[integraseLocusTagsToManuallyCheckIT] = setICEsIMEsStructuresIT
        for integraseToManuallyCheckIT in ICEsIMEsStructuresIT.setIntegraseToManuallyCheck :
            if integraseToManuallyCheckIT.locusTag in __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures:
                setICEsIMEsStructuresIT = __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures[integraseToManuallyCheckIT.locusTag]
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
            else :
                setICEsIMEsStructuresIT = set()
                setICEsIMEsStructuresIT.add(ICEsIMEsStructuresIT)
                __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures[integraseToManuallyCheckIT.locusTag] = setICEsIMEsStructuresIT


    for keyIntegraseAttributedSureUpstreamIT, valueSetICEsIMEsStructureOuter in __integraseAttributedSureUpstream2setICEsIMEsStructures.items():
        if len(valueSetICEsIMEsStructureOuter) > 1 :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 4 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()
        if keyIntegraseAttributedSureUpstreamIT in __integraseAttributedSureDownstream2setICEsIMEsStructures :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 5 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()
            for ICEsIMEsStructureInnerIT in __integraseAttributedSureDownstream2setICEsIMEsStructures[keyIntegraseAttributedSureUpstreamIT] :
                #print("case 6 for {} ({})".format(ICEsIMEsStructureInnerIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureInnerIT.listOrderedSPs)))
                ICEsIMEsStructureInnerIT.transferAllIntegrasesToManuallyCheck()
        if keyIntegraseAttributedSureUpstreamIT.locusTag in __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 7 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()
    for keyIntegraseAttributedSureDownstreamIT, valueSetICEsIMEsStructureOuter in __integraseAttributedSureDownstream2setICEsIMEsStructures.items():
        if len(valueSetICEsIMEsStructureOuter) > 1 :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 8 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()
        if keyIntegraseAttributedSureDownstreamIT in __integraseAttributedSureUpstream2setICEsIMEsStructures :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 9 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()
            for ICEsIMEsStructureInnerIT in __integraseAttributedSureUpstream2setICEsIMEsStructures[keyIntegraseAttributedSureDownstreamIT] :
                #print("case 10 for {} ({})".format(ICEsIMEsStructureInnerIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureInnerIT.listOrderedSPs)))
                ICEsIMEsStructureInnerIT.transferAllIntegrasesToManuallyCheck()
        if keyIntegraseAttributedSureDownstreamIT.locusTag in __integraseLocusTagAttributedToManuallyCheck2setICEsIMEsStructures :
            for ICEsIMEsStructureOuterIT in valueSetICEsIMEsStructureOuter :
                #print("case 11 for {} ({})".format(ICEsIMEsStructureOuterIT.internalIdentifier, hit.ListSPs.GetListProtIdsFromListSP(ICEsIMEsStructureOuterIT.listOrderedSPs)))
                ICEsIMEsStructureOuterIT.transferAllIntegrasesToManuallyCheck()


    return locusTagIntegrase2Comment
