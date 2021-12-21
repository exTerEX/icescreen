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
    :return:           `cds_num`, `cds`, `cds_strand`, `cds_start` and
                       `cds_end`
    :rtype:            :class:`tuple`

    Extracted information
    ---------------------
    cds_num
      Numbering of CDS
    cds
      An identifier of the CDS from genbank file.
      It might be: `protein_id`, `locus_tag`, ...
    cds_strand
      Strand of gene coding for the CDS (+ or -)
    cds_start
      Start position of CDS on genome
    cds_end
      End position of CDS on genome
    """

    # Compile regex to catch all information in ICEscreen CDS identifier
    reg = re.compile(r"^([0-9]*)_(.*)\|(\+|-)\|([0-9]*)\.\.([0-9]*)")

    # Extract information
    reg_results = re.search(reg, identifier)

    return(reg_results.groups(0))


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
