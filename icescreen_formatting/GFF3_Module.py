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
from requests.utils import quote


# https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/
# GFF3 File structure is as follow:
# 9 fields tab separated. All fields must be filled out ('.' if nothing).
# (1) seq id: (mandatory) Name of reference sequence
# (2) source: (mandatory) Source of annotation
# (3) type  : (mandatory) Type of annotation
# (4) start : (mandatory) The 1-based begin position of the annotation
# (5) end   : (mandatory) The 1-based end position of the annotation
# (6) score : (optional ) Score of the annotation (floating point value)
# (7) strand: (mandatory) Strand of the annotation '+' and '-' for forward
#                         and reverse strand, '.' for features not stranded
# (8) phase : (optional ) Shift of feature regarding to the reading frame,
#                         one of "0","1","2" and "." for missing/don't care
# (9) attributes: (opt  ) A list of key/value attributes (key=value)
#                         separated by semicolons.


## global methods

def sanitize_seqid(mystring):
    """Sanitize a string for GFF3 file 9th field.

    :param mystring: String to sanitize
    :type mystring:  str
    :return:         Sanitized string
    :rtype:          str
    """
    reencoded_str = ""

    for c in mystring:
        # [a-zA-Z0-9.:^*$@!+_?-|] characters can be in seqid
        # Else have to escape it by reencoding in percent-encoding
        if re.match("[a-zA-Z0-9.:^*$@!+_?-|]", c) is not None:
            reencoded_str = reencoded_str + c
        else:
            reencoded_str = reencoded_str + quote(c)

    return(reencoded_str)


def write_GFF3_header(file_out_stream):
    file_out_stream.write("##gff-version 3\n")

def write_SeqRecord(file_out_stream, myrecord, sourceSent):

    file_out_stream.write("##sequence-region {} 1 {}\n".format(myrecord.id,
                                                    len(myrecord.seq)
                                                    ))
    
    # Build all columns first
    # 1. seqid
    seqid = sanitize_seqid(myrecord.id)

    # 2. source
    source = sourceSent

    # Other columns
    strand_dict = {1: "+", -1: "-"}
    predefined_attributes = ["id", "name", "alias", "parent", "target", "gap",
                            "derives_from", "note", "dbxref", "ontology_term",
                            "is_circular"]

    for feat in myrecord.features:
        feat_type = feat.type
        #some locus tag have fuzzy start like DS990138.1 VchoM_02810 and get a <
        # if feat.type == "CDS" and "locus_tag" in feat.qualifiers.keys() and feat.qualifiers["locus_tag"][0] == "VchoM_02810":

        temp = 0
        tmpLocationStart = str(feat.location.start)
        for chr in tmpLocationStart:
            if chr.isdigit():
                temp = tmpLocationStart.index(chr)
                break
        prefixStart = tmpLocationStart[0 : temp]
        start = int(feat.location.start) + 1 # feat.location.nofuzzy_start + 1
        start = str(prefixStart) + str(start)
        temp = 0
        tmpLocationEnd = str(feat.location.end)
        for chr in tmpLocationEnd:
            if chr.isdigit():
                temp = tmpLocationEnd.index(chr)
                break
        prefixEnd = tmpLocationEnd[0 : temp]
        end = int(feat.location.end) # feat.location.nofuzzy_end
        end = str(prefixEnd) + str(end)

        score = "."
        strand = strand_dict[feat.location.strand]

        if feat.type == "CDS":
            if "codon_start" in feat.qualifiers.keys():
                # phase = feat.qualifiers["codon_start"]
                try:
                    codon_startIT = int(feat.qualifiers["codon_start"][0])
                    phase = (codon_startIT - 1) % 3
                except Exception:
                    phase = "."
            else:
                phase = "."
        else:
            phase = "."

        attribute_field = []
        for key, value in feat.qualifiers.items():
            # There is some predefined attributes in 9th column of GFF3 file
            # these attributes names are capitalized
            if key.lower() in predefined_attributes:
                attribute_key = key.capitalize()
            else:
                attribute_key = key
            # attribute_value = []
            attribute_value = ""
            if isinstance(value, list):
                attribute_value = ", ".join(value)
            else:
                attribute_value = str(value)

            # There are some qualifiers without values (i.e. "pseudo")
            # if the qualifier exist -> True
            if value == "":
                attribute_value = "True"
            else:
                # Spaces can be in 9th column of GFF3 file
                # So encode each chunk of string in %-encoding then concatenate
                # for val in str(value).split(" "):
                #     attribute_value_part = ""
                #     # attribute_value_part = attribute_value_part + quote(val)
                #     attribute_value_part = attribute_value_part + val
                #     attribute_value.append(attribute_value_part)
                # attribute_value = " ".join(attribute_value)
                list_attribute_value = []
                for val in attribute_value.split(" "):
                    list_attribute_value.append(quote(val))
                attribute_value = " ".join(list_attribute_value)

            attribute_field.append(attribute_key + "=" + attribute_value)

        attribute_field = ";".join(attribute_field)

        file_out_stream.write(f'{seqid}\t{source}\t{feat_type}\t{start}\t{end}\t{score}'
                    f'\t{strand}\t{phase}\t{attribute_field}\n')
        