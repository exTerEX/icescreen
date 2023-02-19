#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ICEscreen copyright Université de Lorraine - INRAE
# This file is part of ICEscreen.
# ICEscreen is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License Affero as published by the Free
# Software Foundation version 3 of the License.
# ICEscreen is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License Affero for more details.
# You should have received a copy of the GNU General Public License Affero
# along with ICEscreen.
# If not, see <https://www.gnu.org/licenses/>.

import re
import sys


def split_identifier(identifier):
    """Extract information from ICEscreen CDS identifier.

    :param identifier: ICEscreen CDS identifier
    :type identifier:  :class:`str`
    :return:           `cds_num`, `locus_tag`, `protein_id`, `genome_accession`, `cds_strand`, `cds_start` and
                       `cds_end`
    :rtype:            :class:`tuple`

    Extracted information
    ---------------------
    cds_num
      Numbering of CDS
    locus_tag
      An uniq identifier of the gene
    protein_ID
      An identifier of the protein.
    genome_accession
      An identifier of the genome record
    genome_accession_rank
      The rank of the genome record in the genbank file
    cds_strand
      Strand of gene coding for the CDS (+ or -)
    cds_start
      Start position of CDS on genome
    cds_end
      End position of CDS on genome
    """

    # Compile regex to catch all information in ICEscreen CDS identifier
    #``[CDS_number]&locus_tag=[locus_tag]&protein_id=[protein_id]&genome_accession=[ACCESSION]&genome_accession_rank=[rank]|[Strand]|[Start]..[End]|[Pseudo]``
    #reg = re.compile(r"^([0-9]*)_(.*)\|(\+|-)\|([0-9]*)\.\.([0-9]*)")
    reg = re.compile(r"^([0-9]*)&locus_tag=(.*)&protein_id=(.*)&genome_accession=(.*)&genome_accession_rank=([0-9]+)\|(\+|-)\|([0-9]*)\.\.([0-9]*)")

    CDS_numToReturn = ""
    locus_tagToReturn = ""
    protein_idToReturn = ""
    genome_accessionToReturn = ""
    genome_accession_rankToReturn = ""
    CDS_strandToReturn = ""
    CDS_startToReturn = ""
    CDS_endToReturn = ""


    # Extract information
    reg_results = re.search(reg, identifier)
    if reg_results is None :
      raise RuntimeError("Error in split_identifier: identifier {} do not match the regex as expected".format(
          str(identifier)))
    else :
      CDS_numToReturn = reg_results.group(1)
      if len(CDS_numToReturn) == 0:
        raise RuntimeError('split_identifier error: len(CDS_numToReturn) == 0 for identifier {}'.format(str(identifier)))
      locus_tagToReturn = reg_results.group(2)
      if len(locus_tagToReturn) == 0:
        locus_tagToReturn = "-"
      protein_idToReturn = reg_results.group(3)
      if len(protein_idToReturn) == 0:
        protein_idToReturn = "-"
      genome_accessionToReturn = reg_results.group(4)
      if len(genome_accessionToReturn) == 0:
        raise RuntimeError('split_identifier error: len(genome_accessionToReturn) == 0 for identifier {}'.format(str(identifier)))
      genome_accession_rankToReturn = reg_results.group(5)
      if len(genome_accession_rankToReturn) == 0:
        raise RuntimeError('split_identifier error: len(genome_accession_rankToReturn) == 0 for identifier {}'.format(str(identifier)))
      CDS_strandToReturn = reg_results.group(6)
      if len(CDS_strandToReturn) == 0:
        raise RuntimeError('split_identifier error: len(CDS_strandToReturn) == 0 for identifier {}'.format(str(identifier)))
      CDS_startToReturn = reg_results.group(7)
      if len(CDS_startToReturn) == 0:
        raise RuntimeError('split_identifier error: len(CDS_startToReturn) == 0 for identifier {}'.format(str(identifier)))
      CDS_endToReturn = reg_results.group(8)
      if len(CDS_endToReturn) == 0:
        raise RuntimeError('split_identifier error: len(CDS_endToReturn) == 0 for identifier {}'.format(str(identifier)))

    #print ("reg_results.groups(0) = {}".format(str(reg_results.groups(0))))
    #print ("locus_tagToReturn = {}".format(str(locus_tagToReturn)))
    #print ("locus_tagToReturn = {}, protein_idToReturn = {}, genome_accessionToReturn = {}, genome_accession_rankToReturn = {}".format(str(locus_tagToReturn), str(protein_idToReturn), str(genome_accessionToReturn), str(genome_accession_rankToReturn)))

    return (CDS_numToReturn, locus_tagToReturn, protein_idToReturn, genome_accessionToReturn, genome_accession_rankToReturn, CDS_strandToReturn, CDS_startToReturn, CDS_endToReturn)
    #return(reg_results.groups(0))


def breakdown_filtering_rule(rule, header_ref):
    """Break down a rule into two information.

    modulation
      Can be "min" or "max"
    field
      Can be one or two fields

    :param rule:       Rule to break down
    :param header_ref: Association table between rules and columns names
    :type rule:        :class:`str`
    :return:           Rule's modulation (:class:`str`) and fields names
                       (:class:`str` or :class:`tuple`)
    :rtype:            :class:`tuple`
    """

    # Modulation and field name are separated by "_"
    try:
        modulation, field_name = rule.split("_", maxsplit=1)
    except ValueError:
        print('ERROR: Filter rule for hit must have this format!: '
              '[min or max]_[parameter name]', file=sys.stderr)
        sys.exit(1)

    if modulation not in ["min", "max"]:
        print('ERROR: Filter rule modulator must be "min" or "max" only!',
              file=sys.stderr)
        sys.exit(1)

    try:
        field = header_ref[field_name]
    except KeyError:
        print(f'ERROR: Parameter "{field_name}" does not exist in '
              f'the header references!', file=sys.stderr)
        sys.exit(1)

    return(modulation, field)
