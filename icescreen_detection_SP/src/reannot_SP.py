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
import sys
import yaml


def parse_arguments():
    """Parse reannot_SP.py script arguments.

    :return: List of paths of following files:

      - (input)  ICEscreen SP results
      - (input)  Tyrosine integrase BlastP hits
      - (input)  YAML configuration
      - (output) ICEscreen SP results with XerS annotated
      - (output) ICEscreen SP results filtered

    :rtype: :class:`list`
    """

    # Create argparse parser
    parser = argparse.ArgumentParser(description="Reannote XerS proteins hits")

    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', help="ICEscreen SP results (.tsv)",
                          required=True)
    required.add_argument('-b', help="Tyrosine integrase hits (.tsv)",
                          required=True)
    required.add_argument('-c', '--config',
                          help="Path to YAML configuration file",
                          required=True)
    required.add_argument('--outall',
                          help="Signature proteins with FP (.tsv)",
                          required=True),
    required.add_argument('--outfiltered',
                          help="Signature proteins without FP",
                          required=True)

    # Parse arguments
    args = parser.parse_args()

    return(args.a, args.b, args.config, args.outall, args.outfiltered)


def get_XerS_annotation(data, tyr):
    """Annotate potential XerS proteins of Firmicutes genomes.

    :param data: SP to evaluate
    :param tyr:  BlastP hits with Tyrosine integrase
    :type data:  :class:`pandas.DataFrame`
    :return:     BlastP hits annotated with potential XerS proteins
    :rtype:      :class:`pandas.DataFrame`

    Each CDS are evaluated and annotated following this procedure:

      - If the CDS has a BlastP hit with at least 70% identity with
        "ACO17137" XerS then it is annotated as "Streptococcal XerS"
      - If the CDS has a BlastP hit with "WP_011835230" XerS protein and is not
        annotated as "Streptococcal XerS" then it might be a XerS protein of
        Firmicutes and will be annotated as "-"
      - If the CDS has a BlastP hit with less than 70% identity with "ACO17137"
        XerS then it might be an invertase and will be annotated
        "Possible invertase"
    """

    # Isolate hits of potential XerS proteins
    candidates = tyr[tyr.Query_blast.str.contains("WP_011835230|ACO17137")]
    # If there is potential XerS proteins
    if not candidates.empty:
        # Initialize list of verified CDS number
        xers_streptococcal = []  # XerS streptococcal proteins
        xers_firmicutes = []     # Other XerS proteins of Firmicutes genomes
        fp = []                  # False positives that will be removed

        # For each CDS, classify into 3 categories
        for candidate in candidates.groupby("CDS_num"):
            cds_num = candidate[0]
            info = candidate[1]

            xers_strep = info[info["Query_blast"] == "ACO17137"]
            xers_firmi = info[info["Query_blast"] == "WP_011835230"]
            if not xers_strep.empty:
                # If there is an hit with ACO17137 and the %id >= 70
                if any(xers_strep.Ali_Identity_perc.values >= 70):
                    xers_streptococcal.append(cds_num)
                # If %id < 70, it is not a streptococcal XerS
                elif any(xers_strep.Ali_Identity_perc.values < 70):
                    # If there is a hit with WP_011835230
                    if not xers_firmi.empty:
                        xers_firmicutes.append(cds_num)
                    else:
                        fp.append(cds_num)
            elif xers_strep.empty:
                if not xers_firmi.empty:
                    xers_firmicutes.append(cds_num)

        # Reannotate data
        for cds_num in xers_streptococcal:
            idx = data[data['CDS_num'] == cds_num].index
            data.loc[idx, "False_positives"] = "Streptococcal XerS"
            data.loc[idx, "Possible_SP"] = "no"

        for cds_num in xers_firmicutes:
            idx = data[data['CDS_num'] == cds_num].index
            #data.shape
            #data.ndim
            data.loc[idx, "False_positives"] = "-"
            data.loc[idx, "Possible_SP"] = "yes"

        for cds_num in fp:
            idx = data[data['CDS_num'] == cds_num].index
            data.loc[idx, "False_positives"] = "Possible invertase"
            data.loc[idx, "Possible_SP"] = "no"

    return(data)


def check_tyrosine_integrase(yamlpath, df, tyr):
    """This function do two things:
    1. Search streptococcal XerS in signature proteins detected and remove it.
    2. Modify the best "Tyrosine integrase" hit if the query used is a XerS

    :param yamlpath: YAML configuration file
    :param df:       Signature proteins detected by BlastP
    :param tyr:      Tyrosine integrase BlastP hits filtered
    :type yamlpath:  :class:`str`
    :type df:        :class:`pandas.DataFrame`
    :type tyr:       :class:`pandas.DataFrame`
    :return:         Input signature proteins without streptococcal XerS
    :rtype:          :class:`pandas.DataFrame`
    """

    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)
    try:
        mode = params_dict["mode"]
    except KeyError:
        print(f'ERROR: "mode" key is missing from "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    # Check tyrosine integrase only in firmicutes genomes
    if mode == "firmicutes":
        # Get annotation of potential XerS proteins
        df = get_XerS_annotation(df, tyr)

        # Reannotate XerS proteins
        df = reannot_XerS(df, tyr)

    return(df)


def reannot_XerS(df, tyr):
    """Reannot Tyrosine integrase annotation which have best hit with ACO17137.
    The function takes the annotation of second best hit possible.

    :param df:  Signature proteins detected by BlastP
    :param tyr: Tyrosine integrase BlastP hits filtered
    :type df:   :class:`pandas.DataFrame`
    :type tyr:  :class:`pandas.DataFrame`
    :return:    Input signature proteins with streptococcal XerS reannotated
    :rtype:     :class:`pandas.DataFrame`
    """
    # Check if there is XerS as best hit
    # Get index of blast hit with ACO17137
    idx = df.loc[df["Query_blast"] == "ACO17137"].index

    if len(idx) == 0:
        return(df)

    # Get 2nd hit of CDS which have ACO17137 as 1st Blast hit
    sec_annot = tyr.loc[tyr["CDS_num"].isin(
            data.iloc[idx]["CDS_num"].values)]
    sec_annot = sec_annot.groupby("CDS_num").nth(1)

    # Replace annotation of ACO17137 blast hit by 2nd best hit if possible
    # Get names of columns to replace
    cols = list(set(sec_annot.columns).intersection(set(df.columns)))

    for i in sec_annot.index:
        for c in cols:
            df.loc[df["CDS_num"] == i, c] = sec_annot.loc[i, c]

    return(df)


if __name__ == "__main__":
    # Parse script arguments
    sppath, tyrpath, conffile, outall, outfiltered = parse_arguments()

    # Import detected SP results
    data = pd.read_csv(sppath, sep="\t")

    # Import Tyrosine Integrase hits results
    tyrint_df = pd.read_csv(tyrpath, sep="\t")

    # If there is no hits
    if len(tyrint_df.index) == 0:
        data["Possible_SP"] = "yes"
        data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")
        data.loc[:, data.columns != "Possible_SP"]\
            .to_csv(outfiltered, index=False, sep="\t", decimal=".",
                    na_rep="NA")
    else:
        data["Possible_SP"] = "yes"
        data = check_tyrosine_integrase(conffile, data, tyrint_df)
        data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")

        # Remove FP
        data = data[data["Possible_SP"] == "yes"]
        data.loc[:, data.columns != "Possible_SP"]\
            .to_csv(outfiltered, index=False, sep="\t", decimal=".",
                    na_rep="NA")
