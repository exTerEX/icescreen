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

import sys
import argparse
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import Seq
from requests.utils import quote


def parse_arguments():
    """Parse generate_annotation_files.py script arguments.

    :return: List of string:

      - (input)  Path of signature proteins data file
      - (input)  Path if mobile elements data file
      - (input)  Path of genome genbank file
      - (output) Name of output annotation files

    :rtype: :class:`list`
    """
    # Parse script arguments
    parser = argparse.ArgumentParser(description="Generate annotation files "
                                                 "for the visualization of "
                                                 "ICEscreen results")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-s', '--spdata', help="SP data (.tsv)",
                          required=True)
    required.add_argument('-m', '--medata', help="Mobile element data (.tsv)",
                          required=True)
    required.add_argument('-g', '--gbfile', help="Genbank file",
                          required=True)
    required.add_argument('-o', '--output', help="Output filename",
                          required=True)

    args = parser.parse_args()

    return(args.spdata, args.medata, args.gbfile, args.output)


def get_CDS_dict(record):
    """Extract CDS features from Bio.SeqRecord object and return a dictionary
    with extracted Bio.SeqFeature as values and CDS start position as key.

    :param record: Genbank record
    :type record:  :class:`Bio.SeqRecord`
    :return:       dict of CDS features (:class:`Bio.SeqFeature.SeqFeature`)
    :rtype:        dict
    """
    # Extract CDS features
    # Key: CDS start position, Value: SeqFeature
    CDS_dict = {}

    for feature in record.features:
        if feature.type == "CDS":
            CDS_dict[feature.location.start.real] = feature

    return(CDS_dict)


def SP_to_SeqFeature(data, CDS_dict):
    """Create a `Bio.SeqFeature.BioSeqFeature` object from a signature protein
    detected by ICEscreen.

    :param data:     Signature protein detected by ICEscreen
    :param CDS_dict: dict of CDS features (:class:`Bio.SeqFeature.SeqFeature`)
    :type data:      :class:`pandas.Series`
    :type CDS_dict:  dict
    :return:         Signatures protein by ICEscreen as feature
    :rtype:          :class:`Bio.SeqFeature.SeqFeature`
    """

    colors_OK = {'Coupling protein': 6,
                 'Relaxase': 7,
                 'Tyrosine integrase': 8,
                 'Serine integrase': 8,
                 'DDE transposase': 8,
                 'VirB4': '184 134 11'}

    colors_pseudo = {'Coupling protein': '139 0 139',
                     'Relaxase': 10,
                     'Tyrosine integrase': 3,
                     'Serine integrase': 3,
                     'DDE transposase': 3,
                     'VirB4': '210 105 30'}

    functions = {'Coupling protein': 'Coupling Protein',
                 'Relaxase': 'Relaxase',
                 'Tyrosine integrase': 'Tyrosine Integrase',
                 'Serine integrase': 'Serine Integrase',
                 'DDE transposase': 'DDE transposase',
                 'VirB4': 'VirB4'}

    featureType = "CDS"

    startPos = int(data["CDS_start"]) - 1
    endPos = int(data["CDS_end"]) - 1

    if data["CDS_strand"] == "+":
        featureLocation = FeatureLocation(startPos, endPos, strand=1)
    else:
        featureLocation = FeatureLocation(startPos, endPos, strand=-1)

    # Parsing of qualifiers of the feature
    myQualifiers = {}

    # ID of CDS
    myQualifiers['origid'] = str(data['CDS']).strip()
    # Type of SP
    myQualifiers['function'] = functions[data['CDS_Protein_type']]

    all_keys = CDS_dict[startPos].qualifiers.keys()

    qualifiers_to_test = ["protein_id", "locus_tag", "gene", "codon_start"]
    for qual in qualifiers_to_test:
        if qual in all_keys:
            myQualifiers[qual] = CDS_dict[startPos].qualifiers[qual][0]

    # Check if gene is pseudogene
    ispseudo = False
    if (data["Is_pseudo"] == "Pseudo") or \
       (data["SP_blast_validation"] == "possible_pseudogene"):
        ispseudo = True
        myQualifiers["pseudo"] = ""
        # Color of the feature
        myQualifiers['color'] = colors_pseudo[data['CDS_Protein_type']]
    else:
        myQualifiers['color'] = colors_OK[data['CDS_Protein_type']]

    # Note of the feature
    # Merge results of Blast and HMM modes
    # Initialization of string "note"
    note = "ICEscreen prediction:"
    # Status of the gene
    if ispseudo is True:
        note = f'{note} Pseudo'

    if data["SP_blast_validation"] in ["manual_check_needed",
                                       "possible_new_protein"]:
        status = " (manual check required)"
    else:
        status = ""

    # Add type of SP
    note = f'{note} {data["CDS_Protein_type"]}{status}'
    # Blast result
    if data["hit_blast"] == 1:
        if data["Use_annotation"] == "yes":
            acs = "HIGH"
        else:
            acs = "LOW"
        note = (f'{note}; BlastP result (Annotation confidence: {acs}): '
                f'{data["Blast_description"]} '
                f'[Hit with {data["Query_blast"]}; ' + "Identity:{0:.2f}%; ".format(data["Ali_Identity_perc"]) + "E-value:{0:.3g}; ".format(data["E-value_blast"]) + "Query coverage:" + "{0:.2f}%]".format(data["Query_blast_coverage"]))
    # HMM result
    if data["hit_HMM"] == 1:
        note = (f'{note}; Hmmscan result: {data["Profile_description"]} '
                f'[Hit with {data["Profile_name"]} HMM profile; ' + f'E-value:{0:.3g}; '.format(data["E-value_hmm"]) + f'i-Evalue:{0:.3g}]'.format(data["i-Evalue_hmm"]))
    myQualifiers['note'] = note

    # Create feature
    myFeature = SeqFeature(featureLocation, type=featureType,
                           qualifiers=myQualifiers)

    return(myFeature)


def get_element_type(elm_info):
    """Get mobile element type based on results of ICEscreen ME detection step.

    :param elm_info: Element information from ICEscreen prediction
    :type elm_info:  :class:`pandas.Series`
    :return:         The element type
    :rtype:          str
    """
    elm_type = ""

    if elm_info["category of element"] == "Complete ICE (4 types of SP)":
        elm_type = f'Putative ICE ({elm_info["category of integrase"]})'
    elif elm_info["category of element"] == "Complete IME (R+I, R+C+I)":
        elm_type = f'Putative IME ({elm_info["category of integrase"]})'
    elif elm_info["category of element"] == "Conjugation module (R+C+V)":
        elm_type = "Complete conjugation module"
    elif "Mobilizable element" in elm_info["category of element"]:
        elm_type = "Complete mobilization module"
    elif elm_info["category of element"] == "Partial ICE (at least V)":
        elm_type = f'Partial ICE '\
                   f'({elm_info["category of integrase"]})'
    elif "Other partial element" in elm_info["category of element"]:
        elm_type = f'Partial mobile element '\
                   f'({elm_info["category of integrase"]})'

    return(elm_type)


def get_element_structure(elm_info):
    """Specify element structure based on element information from ICEscreen
    mobile element detection step results.

    :param elm_info: Element information from ICEscreen prediction
    :type elm_info:  :class:`pandas.Series`
    :return:         The element structure (single or nested)
    :rtype:          str

    Possible element structures are:
      - Single:
          Element is neither nested by another element nor guest of
          another element
      - Nested (host):
          Element hosting another element
      - Nested (guest):
          Element is hosted by another element
      - Nested (host and guest):
          Element is hosting another element
          and is hosted by another element
    """
    element_structure = ""

    if (pd.isna(elm_info["host ICE / IME"])) &\
       (pd.isna(elm_info["guest ICE / IME"])):
        element_structure = "Single"
    elif (pd.isna(elm_info["host ICE / IME"])) &\
         (not pd.isna(elm_info["guest ICE / IME"])):
        element_structure = "Nested (host)"
    elif (not pd.isna(elm_info["host ICE / IME"])) &\
         (pd.isna(elm_info["guest ICE / IME"])):
        element_structure = "Nested (guest)"
    elif (not pd.isna(elm_info["host ICE / IME"])) &\
         (not pd.isna(elm_info["guest ICE / IME"])):
        element_structure = "Nested (host and guest)"

    return(element_structure)


def get_SP_features(CDS_dict, elm_info):
    """Extract the list of signature proteins features of an element predicted
    by ICEscreen.

    :param CDS_dict: All CDS features from one genome
    :param elm_info: Mobile element information from ICEscreen prediction
    :type CDS_dict:  dict of :class:`Bio.SeqFeature.SeqFeature` objects
    :type elm_info:  :class:`pandas.Series`
    :return:         List of signature proteins features of given element
    :rtype:          list of :class:`Bio.SeqFeature.SeqFeature` objects
    """
    if pd.isna(elm_info["list ordered SPs"]):
        return([])

    # SPs of a mobile element are stored in "list ordered SPs" column
    # each SP is written as following: [protein id]-[start position]
    # and each SP are separated by a comma
    ordered_SPs = elm_info["list ordered SPs"].split(", ")
    ordered_SPs_start = []
    for SP_start in ordered_SPs:
        res = re.match("^([^,]*)-([0-9]*)$", SP_start)
        ordered_SPs_start.append(int(res.group(2)) - 1)

    return([CDS_dict[x] for x in ordered_SPs_start])


def get_element_color(elm_info):
    """Specify a color coding based on element information.

    :param elm_info: Mobile element information from ICEscreen prediction
    :type elm_info:  :class:`pandas.Series`
    :return:         The element color coding
    :rtype:          int

    +-------------------------+-------------------+------------+
    | Element type            | Element structure | Color code |
    +=========================+===================+============+
    | - Putative ICE          | - Single          | 15         |
    | - Putative IME          | - Nested (host)   |            |
    | - Complete conjugation  +-------------------+------------+
    |   module                | Nested (guest)    | 4          |
    | - Complete mobilization +-------------------+------------+
    |   module                | Nested (host and  | 1          |
    | - Partial ICE           | guest)            |            |
    +-------------------------+-------------------+------------+
    | Partial mobile element  | /                 | 13         |
    | without integrase       |                   |            |
    +-------------------------+-------------------+------------+
    | Partial mobile element  | /                 | 9          |
    | with integrase          |                   |            |
    +-------------------------+-------------------+------------+
    """

    elm_annot = ["Putative ICE",
                 "Putative IME",
                 "Complete conjugation module",
                 "Complete mobilization module",
                 "Partial ICE"]

    if any([x in elm_info["element_type"] for x in elm_annot]):
        if elm_info["element_structure"] == "Single":
            color = 15
        elif elm_info["element_structure"] == "Nested (host)":
            color = 15
        elif elm_info["element_structure"] == "Nested (guest)":
            color = 4
        elif elm_info["element_structure"] == ("Nested (host and guest)"):
            color = 1
    elif "Partial mobile element" in elm_info["element_type"]:
        if "no integrase" in elm_info["element_type"]:
            color = 13
        else:
            color = 9

    return(color)


def get_element_SP_families(elm_info):
    """Get information about transfer module superfamily and family of a
    mobile element. Also get families of signature proteins of the element.
    Return a dictionary with the following information:
        - blastSuperfamily:
            Mobile element transfer module superfamily based on BlastP results
        - blastFamily:
            Mobile element transfer module family based on BlastP results
        - RFamily:
            Superfamily/Family of relaxase
        - CFamily:
            Superfamily/Family of coupling protein
        - VFamily:
            Superfamily/Family of VirB4
    If an information is missing the value is set as "ND".

    :param elm_info: Mobile element information from ICEscreen prediction
    :type elm_info:  :class:`pandas.Series`
    :return:         Dictionary with information on element family
    :rtype:          dict
    """
    def get_family(family):
        """Check if the information exist. If not return "ND".

        :param family: Superfamily or family information
        :type family:  str (if exist) or float (if missing -> NaN)
        :return:       Superfamily or family information or status
        :rtype:        str
        """
        if isinstance(family, str):
            return(family)
        else:
            return("ND")

    cellName = "ICE SuperFamily From Blast Of SP Conj Module"
    blastSuperfamily = get_family(elm_info[cellName])
    cellName = "ICE Family From Blast Of SP Conj Module"
    blastFamily = get_family(elm_info[cellName])

    hmmFamilies = get_family(elm_info["family SP conj module HMM"])

    RFamily = get_family(elm_info["IME Family From Blast Of SP Conj Module"])

    CFamily = "ND"
    VFamily = "ND"

    if hmmFamilies != "ND":
        hmmFamiliesList = hmmFamilies.split(", ")
        for family in hmmFamiliesList:
            family_info = family.split(":")
            if family_info[0] == "R":
                if RFamily == "ND":
                    RFamily = family_info[1]
                else:
                    RFamily = f'{RFamily} ({family_info[1]})'
            elif family_info[0] == "C":
                CFamily = family_info[1]
            elif family_info[0] == "V":
                VFamily = family_info[1]

    return({"blastSuperfamily": blastSuperfamily,
            "blastFamily": blastFamily,
            "RFamily": RFamily,
            "CFamily": CFamily,
            "VFamily": VFamily})


def make_element_features(CDS_dict, elm_info):
    """Create `Bio.SeqFeature.BioSeqFeature` objects from mobile elements
    detected by ICEscreen.

    :param CDS_dict: dict of CDS features (:class:`Bio.SeqFeature.SeqFeature`)
    :param elm_info: Mobile element information from ICEscreen prediction
    :type CDS_dict:  dict
    :type elm_info:  :class:`pandas.Series`
    :return:         Mobile element detected as feature
    :rtype:          :class:`Bio.SeqFeature.SeqFeature`
    """
    # Get ordered list of SP features of given element
    sp_feats = get_SP_features(CDS_dict, elm_info)

    if len(sp_feats) == 0:
        return(None)

    # Get 1st SP's start index and last SP's end index to get range of element
    element_start_idx = sp_feats[0].location.nofuzzy_start
    element_end_idx = sp_feats[-1].location.nofuzzy_end

    # Then create Bio.SeqFeature.FeatureLocation object
    element_location = FeatureLocation(element_start_idx,
                                       element_end_idx,
                                       strand=1)  # Annotated on + strand
    # Create qualifiers of feature
    myQualifiers = {}

    elm_type = elm_info["element_type"]

    # Create mobile_element_type MANDATORY qualifier
    if "Putative ICE" in elm_type:
        metype = "other: integrative and conjugative element"
        myQualifiers["mobile_element_type"] = metype
        featType = "mobile_element"
    elif "Putative IME" in elm_type:
        metype = "other: integrative and mobilizable element"
        myQualifiers["mobile_element_type"] = metype
        featType = "mobile_element"
    else:
        featType = "misc_feature"

    # Create note
    elm_structure = elm_info["element_structure"]
    elm_ID = elm_info["ICE IME Number"]

    SP_families = get_element_SP_families(elm_info)

    if any([x in elm_type for x in ["Putative ICE", "Partial ICE",
                                    "Complete conjugation module"]]):
        elm_print = "ICE"
        if elm_type == "Complete conjugation module":
            elm_print = "Conjugation module"

        superfamily = SP_families["blastSuperfamily"]
        family = SP_families["blastFamily"]

        note = (f'ICEscreen prediction: {elm_type} [Element structure: '
                f'{elm_structure}; {elm_print} superfamily: {superfamily}; '
                f'{elm_print} family: {family}')

        if SP_families["RFamily"] != "ND":
            note = f'{note}; Relaxase family: {SP_families["RFamily"]}'
        if SP_families["CFamily"] != "ND":
            note = f'{note}; Coupling protein family: {SP_families["CFamily"]}'

    elif any([x in elm_type for x in ["Complete mobilization module",
                                      "Putative IME"]]):
        relfamily = SP_families["RFamily"]

        note = (f'ICEscreen prediction: {elm_type} [Element structure: '
                f'{elm_structure}; Relaxase family: {relfamily}')

        if SP_families["CFamily"] != "ND":
            note = f'{note}; Coupling protein family: {SP_families["CFamily"]}'

    else:
        note = (f'ICEscreen prediction: {elm_type} [Element structure: '
                f'{elm_structure}')

        superfamily = SP_families["blastSuperfamily"]
        family = SP_families["blastFamily"]

        if ((superfamily != "ND" and family != "ND") or (superfamily != "ND" and family == "ND") or (superfamily == "ND" and family != "ND")):
            note = (f'ICEscreen prediction: {elm_type} '
                    f'[Element structure: {elm_structure}; '
                    f'Transfer module superfamily: {superfamily}; '
                    f'Transfer module family: {family}')

        if SP_families["RFamily"] != "ND":
            note = f'{note}; Relaxase family: {SP_families["RFamily"]}'
        if SP_families["CFamily"] != "ND":
            note = f'{note}; Coupling protein family: {SP_families["CFamily"]}'

    note = f'{note}] (ICEscreen ID: {elm_ID})'

    myQualifiers["note"] = note

    # Create color qualifier (for visualization with Artemis)
    myQualifiers["color"] = get_element_color(elm_info)

    # Create feature
    myFeature = SeqFeature(element_location,
                           type=featType,
                           qualifiers=myQualifiers)

    return(myFeature)


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


def write_GFF3(myrecord, seq_length, outfile):
    """Write ICEscreen results as annotation in GFF3 file.

    :param myrecord:   Genbank information as Bio.SeqRecord.SeqRecord object
    :param seq_length: Genome length
    :param outfile:    Path of output GFF3 file
    :type myrecord:    :class:`Bio.SeqRecord.SeqRecord`
    :type seq_length:  int
    :type outfile:     str
    :return:           None
    """
    filout = open(outfile, "w")

    # GFF3 header
    filout.write("##gff-version 3\n")
    filout.write("##sequence-region {} 1 {}\n".format(myrecord.id,
                                                      seq_length))

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

    # Build all columns first
    # 1. seqid
    seqid = sanitize_seqid(myrecord.id)

    # 2. source
    source = "ICEscreen"

    # Other columns
    strand_dict = {1: "+", -1: "-"}
    predefined_attributes = ["id", "name", "alias", "parent", "target", "gap",
                             "derives_from", "note", "dbxref", "ontology_term",
                             "is_circular"]

    for feat in myrecord.features:
        feat_type = feat.type
        start = feat.location.nofuzzy_start + 1
        end = feat.location.nofuzzy_end
        score = "."
        strand = strand_dict[feat.location.strand]

        if feat.type == "CDS":
            if "codon_start" in feat.qualifiers.keys():
                phase = feat.qualifiers["codon_start"]
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
            attribute_value = []
            # There are some qualifiers without values (i.e. "pseudo")
            # if the qualifier exist -> True
            if value == "":
                attribute_value = "True"
            else:
                # Spaces can be in 9th column of GFF3 file
                # So encode each chunk of string in %-encoding then concatenate
                for val in str(value).split(" "):
                    attribute_value_part = ""
                    attribute_value_part = attribute_value_part + quote(val)
                    attribute_value.append(attribute_value_part)
                attribute_value = " ".join(attribute_value)

            attribute_field.append(attribute_key + "=" + attribute_value)

        attribute_field = ";".join(attribute_field)

        filout.write(f'{seqid}\t{source}\t{feat_type}\t{start}\t{end}\t{score}'
                     f'{strand}\t{phase}\t{attribute_field}')

    filout.close()


if __name__ == "__main__":
    spPath, mePath, gbPath, outPath = parse_arguments()

    # Get SeqFeature of all SPs detected
    spdata = pd.read_csv(spPath, sep="\t")

    if len(spdata.index) == 0:
        open(outPath + ".embl", "a").close()
        open(outPath + ".gff", "a").close()
        open(outPath + ".gb", "a").close()
        sys.exit(0)

    # Get SeqFeature of all CDS of the genbank
    gbdata = SeqIO.read(gbPath, 'gb')
    CDS_dict = get_CDS_dict(gbdata)
    sp_feats = spdata.apply(lambda x: SP_to_SeqFeature(x, CDS_dict),
                            axis=1).tolist()

    # Get annotation of mobile elements that have all SPs as "assigned"
    data_tsv = pd.read_csv(mePath, sep="\t")

    if len(data_tsv.index) > 0:
        data_tsv["element_type"] = data_tsv.apply(get_element_type, axis=1)
        data_tsv["element_structure"] = data_tsv.apply(get_element_structure,
                                                       axis=1)

        # Get SeqFeature of detected mobile elements
        me_feats = data_tsv.apply(lambda x: make_element_features(CDS_dict, x),
                                  axis=1).to_list()
        me_feats = [x for x in me_feats if x is not None]

        # Create SeqRecord with "light" annotation (without FASTA sequence)
        myrecord = SeqRecord(seq=Seq.Seq(''),
                             id=gbdata.id,
                             name=gbdata.name,
                             description=gbdata.description,
                             dbxrefs=gbdata.dbxrefs,
                             features=sp_feats + me_feats,
                             annotations={"molecule_type": "DNA"})
    else:
        # Create SeqRecord with "light" annotation (without FASTA sequence)
        myrecord = SeqRecord(seq=Seq.Seq(''),
                             id=gbdata.id,
                             name=gbdata.name,
                             description=gbdata.description,
                             dbxrefs=gbdata.dbxrefs,
                             features=sp_feats,
                             annotations={"molecule_type": "DNA"})

    # Save as EMBL
    with open(outPath + ".embl", "w") as output_handle:
        SeqIO.write(myrecord, output_handle, "embl")

    # And save as GFF3
    write_GFF3(myrecord, len(gbdata.seq), outPath + ".gff")

    if data_tsv.empty is False:
        # Create SeqRecord with "heavy" annotation (with FASTA sequence)
        # Sequence of the SeqRecord is extracted from en genbank
        myrecord = SeqRecord(seq=gbdata.seq,
                             id=gbdata.id,
                             name=gbdata.name,
                             description=gbdata.description,
                             dbxrefs=gbdata.dbxrefs,
                             features=gbdata.features + sp_feats + me_feats,
                             annotations={"molecule_type": "DNA"})
    else:
        # Create SeqRecord with "heavy" annotation (with FASTA sequence)
        # Sequence of the SeqRecord is extracted from en genbank
        myrecord = SeqRecord(seq=gbdata.seq,
                             id=gbdata.id,
                             name=gbdata.name,
                             description=gbdata.description,
                             dbxrefs=gbdata.dbxrefs,
                             features=gbdata.features + sp_feats,
                             annotations={"molecule_type": "DNA"})

    # Save as genbank
    with open(outPath + ".gb", "w") as output_handle:
        SeqIO.write(myrecord, output_handle, "gb")
