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
import pandas as pd


def parse_arguments():
    """Parse process_blastp_results.py script arguments.

    :return: List of paths of following files:

      - (input)  ICEscreen BlastP results
      - (input)  YAML configuration
      - (input)  ICEscreen annotation database
      - (output) ICEscreen BlastP results annotated
      - (output) All possible Signature Proteins identified
      - (output) Best Signature Proteins-like identified

    :rtype: :class:`list`
    """

    # Create argparse parser
    parser = argparse.ArgumentParser(description="Process ICEscreen BlastP "
                                                 "results")

    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input',
                          help="ICEscreen BlastP results (.tsv) after Concatenate all results and Reorder rows in the snakemake rule gather_blastP_SP_results.",
                          required=True)
    required.add_argument('-o', '--output',
                          help="ICEscreen BlastP results (.tsv) with only the best hits for each locus tag.",
                          required=True)

    # Parse arguments
    args = parser.parse_args()

    return(args.input, args.output)

def main():
    # Parse script arguments
    (input_file, output_file) = parse_arguments()

    blast_df = pd.read_csv(input_file, sep="\t")

    blast_df = blast_df.sort_values(by=["Genome_accession", "CDS_num", "CDS_locus_tag", "Blast_ali_E-value", "Blast_ali_bitscore", "CDS_coverage", "Blast_ali_coverage_most_similar_ref_SP", "Blast_ali_identity_perc"], ascending=[True, True, True, True, False, False, False, False])

    blast_df=blast_df.drop_duplicates(subset=["Genome_accession", "CDS_num", "CDS_locus_tag"], keep='first')

    blast_df.to_csv(output_file, index=False, sep="\t", decimal=".", na_rep="NA")

if __name__ == "__main__":
    main()