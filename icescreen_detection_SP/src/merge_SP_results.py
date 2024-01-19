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
import numpy as np


def parse_arguments():
    """Parse merge_SP_results.py script arguments.

    :return: List of paths of following files:

      - (input)  ICEscreen BlastP results
      - (input)  ICEscreen hmmscan results
      - (output) ICEscreen merged results

    :rtype: :class:`list`
    """

    parser = argparse.ArgumentParser(description="Merge results of ICEscreen "
                                                 "SP scan")

    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('--blastres',
                          help="SP detected with BlastP (.tsv)",
                          required=True)
    required.add_argument('--hmmres',
                          help="SP detected with hmmscan (.tsv)",
                          required=True)
    required.add_argument('-o', '--outtsv',
                          help="SP detected with BlastP and hmmscan (.tsv)",
                          required=True)

    args = parser.parse_args()

    return(args.blastres, args.hmmres, args.outtsv)


def merge_results(blpath, hmmpath):
    """Order signature proteins by position on the genome and sort each result
    from best to worst. Sanitize string values if needed and format numeric
    values. Rename columns to more explicit terms.

    :param df:       Dataframe to work with
    :type df:        :class:`pandas.DataFrame`
    :return:         Dataframe prettified
    :rtype:          :class:`pandas.DataFrame`
    """

    # Open results
    blast_df = pd.read_csv(blpath, sep="\t")
    hmm_df = pd.read_csv(hmmpath, sep="\t")

    # Add 'Is_hit_blast' and 'Is_hit_HMM' columns
    blast_df["Is_hit_blast"] = 1
    hmm_df["Is_hit_HMM"] = 1

    merged_df = pd.merge(left=blast_df, right=hmm_df, how='outer',
                         on=["#ICEscreen_ID", "CDS_num",
                         #"CDS",
                         "Genome_accession", "Genome_accession_rank", "CDS_locus_tag", "CDS_protein_id",
                         "CDS_strand", "CDS_start", "CDS_end", "CDS_length"],
                         suffixes=('_blast', '_hmm'))

    # Fill out 'Is_hit_blast' and 'Is_hit_HMM' columns
    merged_df["Is_hit_blast"] = merged_df["Is_hit_blast"].replace({None: 0})
    merged_df["Is_hit_HMM"] = merged_df["Is_hit_HMM"].replace({None: 0})
    merged_df["Is_hit_blast"] = merged_df["Is_hit_blast"].fillna(0)
    merged_df["Is_hit_HMM"] = merged_df["Is_hit_HMM"].fillna(0)

    # Add 'Is_pseudo' column
    is_pseudo = merged_df["#ICEscreen_ID"].str.endswith("Pseudo")\
                                          .replace({True: "Pseudo",
                                                    False: "-"})
    merged_df["Is_pseudo"] = is_pseudo

    # Merge CDS Protein type prediction
    cds_protein_type = merged_df[["CDS_Protein_type_blast",
                                  "CDS_Protein_type_hmm"]].bfill(axis=1)\
                                                          .iloc[:, 0]
    merged_df["CDS_Protein_type"] = cds_protein_type

    # Rename some columns to more specific name
    merged_df.rename(columns={"CDS_Protein_type_hmm": "Protein_type_of_matching_HMM_profile"},
                     inplace=True)

    # Drop unnecessary columns
    # Drop raw results of BlastP alignments
    # cols_to_rm = ["Mismatch", "Gapopen", "Blast_ali_start_CDS", "Blast_ali_end_CDS",
    #               'Blast_ali_start_Query_blast', 'Blast_ali_end_Query_blast',
    #               'CDS_length/Length_of_blast_most_similar_ref_SP']
    cols_to_rm = ["Mismatch", "Gapopen", 'CDS_length/Length_of_blast_most_similar_ref_SP']
    # Drop raw results of hmmscan alignments and Profile ID
    cols_to_rm = cols_to_rm + ["Profile_ID", "Ali_len_HMM",
                               "Ali_len_CDS", "hmm_coord_from", "hmm_coord_to",
                               "ali_coord_from", "ali_coord_to", "#_domain",
                               "of_domain"]
    # Drop cutoff columns of SP detected with hmmscan
    cols_to_rm = cols_to_rm + list(hmm_df.columns[hmm_df.columns.str
                                                        .endswith("_cutoff")])
    cols_to_rm = cols_to_rm + list(hmm_df.columns[hmm_df.columns.str
                                                        .endswith("_OK")])
    # Drop ICEscreen ID for CDS (redundant with other fields)
    cols_to_rm = cols_to_rm + ["#ICEscreen_ID"]

    merged_df = merged_df.loc[:, ~merged_df.columns.isin(cols_to_rm)]

    # Reorder columns: First "core" common columns then other columns
    # Common columns between blastP results and hmmscan results
    common_cols = ["Is_hit_blast", "Is_hit_HMM", "CDS_num",
                    #"CDS",
                    "Genome_accession", "Genome_accession_rank", "CDS_locus_tag", "CDS_protein_id",
                    "CDS_strand", "CDS_start", "CDS_end", "CDS_length",
                    "Is_pseudo", "CDS_Protein_type"]
    # Other columns
    other_cols = [x for x in merged_df.columns if x not in common_cols]

    merged_df = merged_df[common_cols + other_cols]

    merged_df = merged_df.astype({'Is_hit_blast': np.uint32,
                                    'Is_hit_HMM': np.uint32})

    # Format numeric columns
    merged_df['HMM_ali_E-value'] = merged_df['HMM_ali_E-value'].map(
            lambda x: '{0:.3g}'.format(x))
    merged_df['HMM_ali_i-Evalue'] = merged_df['HMM_ali_i-Evalue'].map(
            lambda x: '{0:.3g}'.format(x))

    return(merged_df)


if __name__ == "__main__":

    # Parse script arguments
    blastres, hmmres, outtsv = parse_arguments()

    merged_results = merge_results(blastres, hmmres)

    # Order SP by position on the genome
    #merged_results = merged_results.sort_values(by="CDS_num", ascending=True)
    merged_results = merged_results.sort_values(by=["Genome_accession", "CDS_num"], ascending=[True, True])

    merged_results.to_csv(outtsv, index=False, sep="\t", decimal=".",
                          na_rep="NA")
