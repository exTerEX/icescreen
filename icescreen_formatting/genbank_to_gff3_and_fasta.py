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
# import BCBio.GFF as GFF
from icescreen_formatting import GFF3_Module


def parse_arguments():
    """Parse genbank_to_gff3_and_fasta.py script arguments.

    :return: List of paths of following files:

      - (input)  Genbank file
      - (output) GFF annotation file
      - (output) FASTA file

    :rtype: :class:`list`
    """
    # Parse script arguments
    parser = argparse.ArgumentParser(description="Convert Genbank to GFF3 "
                                                 "(without FASTA) and "
                                                 "extract sequence to FASTA.")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', help="Genbank file", required=True)
    required.add_argument('--gff', help="GFF3 output", required=True)
    required.add_argument('--fasta', help="Fasta output", required=True)

    args = parser.parse_args()

    return(args.input, args.gff, args.fasta)


def main():
    gb_file, gff_file, fasta_file = parse_arguments()

    with open(gff_file, "w") as gff_out:
        # filin = SeqIO.parse(gb_file, "genbank")
        # GFF.write(filin, gff_out, include_fasta=False)
        GFF3_Module.write_GFF3_header(gff_out)
        record_iterator = SeqIO.parse(gb_file, "genbank")
        for gbdata in record_iterator:
            GFF3_Module.write_SeqRecord(gff_out, gbdata, "feature")

    with open(fasta_file, "w") as fasta_out:
        #filin = SeqIO.read(gb_file, "genbank")
        #SeqIO.write(filin, fasta_out, "fasta")
        record_iterator = SeqIO.parse(gb_file, "genbank")
        listRecordsToWrite = []
        for gbdata in record_iterator:
            listRecordsToWrite.append(gbdata)
        SeqIO.write(listRecordsToWrite, fasta_out, "fasta")


if __name__ == "__main__":
    main()
