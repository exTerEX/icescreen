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

import argparse
from Bio import SeqIO
import yaml
import sys
import Bio.Data.CodonTable as CodonTable


def parse_arguments():
    """Parse gb_to_faa.py script arguments.

    :return: List of paths of following files:

      - (input)  Genbank file
      - (input)  YAML configuration file
      - (output) Multifasta file

    :rtype:  :class:`list`
    """

    parser = argparse.ArgumentParser(description="Create multifasta protein "
                                                 "from a genbank file")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', help="Genbank file", required=True)
    required.add_argument('-o', '--output', help="Output file", required=True)
    required.add_argument('-c', '--config',
                          help="Path to YAML configuration file",
                          required=True)
    args = parser.parse_args()

    return(args.input, args.config, args.output)


def get_codon_table(yamlpath):
    """Get codon table number for translation of coding sequences.

    :param yamlpath: YAML configuration file
    :type yamlpath:  :class:`str`
    :return:         Codon table number
    :rtype:          :class:`int`
    """

    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    try:
        codon_table = params_dict["SP_detection"]["cds_extraction"]
    except KeyError:
        print(f'ERROR: "cds_extraction" parameters are missing in "{yamlpath}"'
              f' config file!',
              file=sys.stderr)
        sys.exit(1)
    try:
        codon_table = codon_table["codon_table"]
    except KeyError:
        print(f'ERROR: "codon_table" is missing in "cds_extraction" parameters'
              f' in "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    return(codon_table)


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
    res = res.replace(',', '')
    res = res.replace(' ', '_')
    res = res.replace('|', '_')

    return(res)


def gb_to_faa(gbfile, faafile, codontable):
    """Extract CDS from genbank file and save into target multifasta file.

    :param gbfile: Path of genbank file
    :param faafile: Path of multifasta file
    :type gbfile: :class:`str`
    :type faafile: :class:`str`
    :return: None

    ===================================
    Sequence identifier (#ICEscreen_ID)
    ===================================
    An identifier is generated for each CDS with 6 information:

    Format
    ------
      ``[CDS_number]_[Protein ID]|[Strand]|[Start]..[End]|[Pseudo]``

    CDS number
      The number of the CDS in the genome
    Protein ID
      An identifier of the protein.

      By default it should be the `protein_id` qualifier but some genbank
      might not have one.

      So 4 qualifiers are searched successively:

      1. `protein_id`
      2. `label`
      3. `locus_tag`
      4. `db_xref`

    Strand
      Strand of the CDS
    Start
      Start position of the CDS
    End
      End position of the CDS
    Pseudo
      If it is a pseudogene, else this field does not exist

    Example
    -------
    For CDS with functional gene
      ``1_WP_000138202.1|+|176..1537``

    For CDS with pseudogene
      ``463_GBS_RS02765|-|491777..492813|Pseudo``
    """

    filin = open(gbfile, "r")

    # Check if there is multiple records (multigenbank)
    try:
        record = SeqIO.read(filin, "genbank")
        filin.close()
    except ValueError as e:
        filin.close()
        if "More than one record found in handle" == str(e):
            print("ERROR: Genbank file should contain only ONE record "
                  "(multigenbank files are not allowed): {gbfile}",
                  file=sys.stderr)
            sys.exit(1)
        else:
            print(f"ERROR: Cannot read genbank file: {gbfile}",
                  file=sys.stderr)
            sys.exit(1)

    # CDS counter for numbering features along the genome
    # Assuming features of the seqRecord are ordered by location on the genome
    counter = 1

    # File out: multifasta
    filout = open(faafile, "w")

    # Parsing sequences with SeqIO
    for feature in record.features:
        # Get location (start, end and strand)
        # start + 1 to transform index into location on the chromosome
        start = feature.location.start + 1
        end = feature.location.end
        strand = feature.location.strand

        # Initialization of cdsSeq and protId variables
        cdsSeq = 'Genbank file is not a Genbank full'
        protId = ''

        # Format information of location
        if strand == 1:
            strandInfo = "+|{:d}..{:d}".format(start, end)
        if strand == -1:
            strandInfo = "-|{:d}..{:d}".format(start, end)

        # Get only CDS
        if feature.type == "CDS":
            # If pseudogene might have "pseudo" or "pseudogene" qualifier
            if any(x in feature.qualifiers for x in ["pseudo", "pseudogene"]):
                pseudo = "|Pseudo"
            else:
                pseudo = ""

            # Get CDS protein sequence: If there is already a "translation"
            # qualifier, we have the CDS protein sequence
            if 'translation' in feature.qualifiers:
                cdsSeq = feature.qualifiers['translation'][0]
            else:
                nucSeq = feature.location.extract(record).seq
                # Check if CDS or pseudogene
                try:
                    cdsSeq = nucSeq.translate(table=codon_table,
                                              cds=True)
                # If it's an Translation Error then it's a pseudogene
                except Exception as e:
                    if isinstance(e, CodonTable.TranslationError):
                        pseudo = "|Pseudo"
                        # Translate sequence without CDS parameter
                        cdsSeq = nucSeq.translate(table=codon_table,
                                                  cds=False)

            # Search protein ID in 4 qualifiers
            if 'protein_id' in feature.qualifiers:
                protId = search_feature('protein_id', feature)
            elif 'label' in feature.qualifiers:
                protId = search_feature('label', feature)
            elif 'locus_tag' in feature.qualifiers:
                protId = search_feature('locus_tag', feature)
            elif 'db_xref' in feature.qualifiers:
                protId = search_feature('db_xref', feature)

            # Write the fasta in multifasta file
            filout.write(f'>{counter}_{protId}|{strandInfo}{pseudo}\n'
                         f'{cdsSeq}\n')

            # Increment counter
            counter += 1

    filout.close()


if __name__ == '__main__':
    inputfile, conffile, outputfile = parse_arguments()

    codon_table = get_codon_table(conffile)

    gb_to_faa(inputfile, outputfile, codon_table)
