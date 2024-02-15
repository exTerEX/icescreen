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

from process_hmmscan_results import get_hmmscan_df, get_params_set_hmmscan
from process_hmmscan_results import filter_df, pretty_df, reorder_columns
import argparse
import yaml
import pandas as pd
import sys


def parse_arguments_fp():
    """Parse process_hmmscan_fp.py script arguments.

    :return: List of paths of following files:

      - (input)  ICEscreen SP results
      - (input)  hmmscan hits of false positives
      - (input)  YAML configuration
      - (output) hmmscan false positives annotated
      - (output) ICEscreen SP without detected False positives

    :rtype: :class:`list`
    """

    parser = argparse.ArgumentParser(description="Process ICEscreen hmmscan "
                                                 "results")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('--insp', help="TSV file with SP hits (.tsv)",
                          required=True)
    required.add_argument('--infp', help="hmmscan results for false positives "
                          "(.tsv)", required=True)
    required.add_argument('-c', '--config',
                          help="Path to YAML configuration file",
                          required=True)
    required.add_argument('--outfp',
                          help="False positive hits detected (.tsv)",
                          required=True),
    required.add_argument('--outfiltered',
                          help="SP hits without detected False positives "
                          "(.tsv)",
                          required=True)
    args = parser.parse_args()

    return(args.insp, args.infp, args.config, args.outfp, args.outfiltered)


def pretty_df_sp(df):
    """Order signature proteins by position on the genome and sort each result
    from best to worst. Sanitize string values if needed and format numeric
    values. Rename columns to more explicit terms.

    :param df:       Dataframe to work with
    :type df:        :class:`pandas.DataFrame`
    :return:         Dataframe prettified
    :rtype:          :class:`pandas.DataFrame`
    """

    # Format numeric columns
    df['HMM_ali_E-value'] = df['HMM_ali_E-value'].map(lambda x: '{0:.3g}'.format(x))
    df['HMM_ali_i-Evalue'] = df['HMM_ali_i-Evalue'].map(lambda x: '{0:.3g}'.format(x))
    df['HMM_coverage'] = df['HMM_coverage'].map(lambda x: '{0:.2f}'.format(x))
    df['CDS_coverage_hmm'] = df['CDS_coverage_hmm'].map(
            lambda x: '{0:.2f}'.format(x))
    for col in df.columns[df.columns.str.endswith("_cutoff")]:
        df[col] = df[col].map(lambda x: '{0:.1f}'.format(x))

    return(df)


if __name__ == "__main__":

    # Parse script arguments
    insp, infp, conffile, outfp, outfiltered = parse_arguments_fp()

    # Get hmmscan results
    data_sp = pd.read_csv(insp, sep="\t")
    data_sp = pretty_df_sp(data_sp)

    data_fp = get_hmmscan_df(infp)

    # If no results: create empty tsv files
    if len(data_fp.index) == 0:
        cols = ["#ICEscreen_ID", "CDS_num",
                #"CDS",
                "Genome_accession", "Genome_accession_rank", "CDS_locus_tag", "CDS_protein_id",
                "CDS_strand", "CDS_start", "CDS_end", "CDS_length", "CDS_Protein_type",
                "Description_of_matching_HMM_profile", "Profile_ID", "Profile_name",
                "Length_of_matching_HMM_profile", "Ali_len_HMM", "Ali_len_CDS",
                "hmm_coord_from", "hmm_coord_to", "ali_coord_from",
                "ali_coord_to", "HMM_ali_i-Evalue", "HMM_ali_E-value", "HMM_ali_Score",
                "HMM_ali_Bias", "#_domain", "of_domain", "HMM_ali_Global_score",
                "HMM_ali_Global_bias", "HMM_coverage", "CDS_coverage",
                "i-Evalue_cutoff", "i-Evalue_OK", "CDS_coverage_cutoff",
                "CDS_coverage_OK", "HMM_coverage_cutoff", "HMM_coverage_OK"]

        data_fp = pd.DataFrame(columns=cols)

        data_fp.to_csv(outfp, index=False, sep="\t", decimal=".", na_rep="NA")
        data_sp.to_csv(outfiltered, index=False, sep="\t", decimal=".",
                       na_rep="NA")

        sys.exit(0)

    with open(conffile, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    # Parameters for hmmscan results
    params_hmmscan = params_dict["SP_detection"]["hmmscan"]

    # Get dictionnary to make correspondance between fields in header and
    # rule name in configuration file
    header_ref = {"query_length": "CDS_length",
                  "ievalue": "i-Evalue",
                  "cds_coverage": "CDS_coverage",
                  "profile_coverage": "HMM_coverage"}

    # Assess hmmscan results for each HMM profile and SP type
    results = []

    for group in data_fp.groupby(by="Profile_name"):
        # To filter out non significant results need to know:
        # To get rules for filtering should have
        # and hmmprofile used
        hmmprofile = group[1]["Profile_name"].iloc[0]

        # Get parameters for filtering
        params_set = get_params_set_hmmscan(conffile,
                                            "params_profiles",
                                            hmmprofile)

        # Add results of filtering for each group
        results.append(filter_df(params_set, header_ref, group[1].copy()))

    # Concatenate all results
    data_fp = pd.concat(results, ignore_index=True)

    # Get protein type
    data_fp["CDS_Protein_type"] = "False Positive"
    # Remove "prot_type" column
    data_fp = data_fp.loc[:, data_fp.columns != "prot_type"]

    data_fp = pretty_df(data_fp)

    # All results without filtering
    data_fp = reorder_columns(data_fp.copy(deep=True))

    # All validated results are false positives
    data_fp = data_fp[data_fp["Possible_SP"] == "yes"]
    data_fp.loc[:, data_fp.columns != "Possible_SP"]\
           .to_csv(outfp, index=False, sep="\t", decimal=".", na_rep="NA")

    # Reduce number of lines
    data_fp = data_fp[~data_fp.CDS_num.duplicated()]

    if data_fp.empty is False:
        # Remove false positives from input TSV file
        data_sp = data_sp[~data_sp["CDS_num"].isin(
                data_fp["CDS_num"].values.tolist())]

    # Order SP by position on the genome and
    # sort each result from best to worst
    data_sp = data_sp.sort_values(by=["Genome_accession_rank", "CDS_num"], ascending=[True, True])

    data_sp.to_csv(outfiltered, index=False, sep="\t", decimal=".",
                   na_rep="NA")
