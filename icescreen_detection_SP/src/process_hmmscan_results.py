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
import yaml
import numpy as np
import pandas as pd
from process_SP_hits import split_identifier, breakdown_filtering_rule


def parse_arguments():
    """Parse process_hmmscan_results.py script arguments.

    :return: List of paths of following files:

      - (input)  ICEscreen hmmscan results
      - (input)  YAML configuration
      - (output) ICEscreen hmmscan results annotated
      - (output) All possible Signature Proteins identified
      - (output) Best Signature Proteins-like identified

    :rtype: :class:`list`
    """

    parser = argparse.ArgumentParser(description="Process ICEscreen hmmscan "
                                                 "results")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', help="tsv file (.tsv)",
                          required=True)
    required.add_argument('-c', '--config',
                          help="Path to YAML configuration file",
                          required=True)
    required.add_argument('--outall',
                          help="All hmmscan results with annotation (.tsv)",
                          required=True),
    required.add_argument('--outfiltered',
                          help="All possible Signature Proteins identified"
                               " (.tsv)",
                          required=True)
    required.add_argument('--outbest',
                          help="Best Signature Proteins-like identified"
                          " (.tsv)",
                          required=True)
    args = parser.parse_args()

    return(args.input, args.config, args.outall, args.outfiltered,
           args.outbest)


def split_target_name(target_name):
    """Extract the protein type and family of the target name in hmmscan
    results. The target name is the name of the HMM profile used in hmmscan
    result. HMM profile name is formatted like this:
        [protein type]_[family]

    :param target_name: Name of HMM profile used in hmmscan result
    :type target_name:  :class:`str`
    :return:            Type of protein and its family
    :rtype:             :class:`tuple`
    """

    protein_type, family = target_name.split("_", maxsplit=1)

    return(protein_type, family)


def get_hmmscan_df(tsvpath):
    """Parse hmmscan results from TSV file and convert into `pandas.DataFrame`.
    TSV file must have these 23 following columns:
        - target_name
        - accession
        - tlen
        - query_name
        - seq_accession
        - qlen
        - seq_E-value
        - seq_score
        - seq_bias
        - #_domain
        - of_domain
        - domain_c-Evalue
        - domain_i-Evalue
        - domain_score
        - domain_bias
        - hmm_coord_from
        - hmm_coord_to
        - ali_coord_from
        - ali_coord_to
        - env_coord_from
        - env_coord_to
        - env_coord_acc
        - description_of_target

    :param tsvpath: TSV file path
    :type tsvpath:  :class:`str`
    :return:        TSV information formatted as dataframe
    :rtype:         :class:`pandas.DataFrame`
    """

    columns_type = {"target_name": str,
                    "accession": str,
                    "tlen": np.uint32,
                    "query_name": str,
                    "seq_accession": str,
                    "qlen": np.uint32,
                    "seq_E-value": np.float32,
                    "seq_score": np.float32,
                    "seq_bias": np.float32,
                    "#_domain": np.uint32,
                    "of_domain": np.uint32,
                    "domain_c-Evalue": np.float32,
                    "domain_i-Evalue": np.float32,
                    "domain_score": np.float32,
                    "domain_bias": np.float32,
                    "hmm_coord_from": np.uint32,
                    "hmm_coord_to": np.uint32,
                    "ali_coord_from": np.uint32,
                    "ali_coord_to": np.uint32,
                    "env_coord_from": np.uint32,
                    "env_coord_to": np.uint32,
                    "env_coord_acc": np.uint32,
                    "description_of_target": str}

    try:
        data = pd.read_csv(tsvpath,
                           sep="\t",
                           dtype=columns_type,
                           header=0,
                           # Comment if want to allow additional columns
                           names=columns_type.keys()
                           )
    except ValueError:
        data = pd.read_csv(tsvpath, "\t")
    except KeyError:
        print('ERROR: There is a problem in TSV file header. '
              '(Missing / Misspelled / Additional columns).',
              file=sys.stderr)
        sys.exit(1)

    # If no results
    if len(data.index) == 0:
        return(data)

    # Split HMM profile name into protein_type and family
    columns = ["prot_type", "Profile_name"]
    targets_info = np.vectorize(split_target_name)(data["target_name"])
    targets_info = np.transpose(targets_info)
    data[columns] = pd.DataFrame(targets_info, columns=columns)

    # Split ICEscreen CDS identifier into multiple columns
    columns = ["CDS_num", "CDS", "CDS_strand", "CDS_start", "CDS_end"]
    cds_identifiers = np.vectorize(split_identifier)(data["query_name"])
    cds_identifiers = np.transpose(cds_identifiers)
    data[columns] = pd.DataFrame(cds_identifiers, columns=columns)

    # Compute alignment length for each hit of hmmscan (hits are DOMAINS)
    data["Ali_len_HMM"] = data["hmm_coord_to"] - data["hmm_coord_from"] + 1
    data["Ali_len_CDS"] = data["ali_coord_to"] - data["ali_coord_from"] + 1

    # Compute coverages
    data["HMM_coverage"] = data["Ali_len_HMM"] / data["tlen"] * 100
    data["CDS_coverage"] = data["Ali_len_CDS"] / data["qlen"] * 100

    # Drop unused columns
    data = data.loc[:, ~data.columns.isin(["domain_c-Evalue", "accession",
                                           "seq_accession"])]
    data = data.loc[:, ~data.columns.str.startswith("env_coord_")]

    # Rename some columns
    data = data.rename(columns={"qlen": "CDS_length",
                                "domain_i-Evalue": "i-Evalue",
                                "seq_E-value": "E-value"})

    return(data)


def get_params_set_hmmscan(yamlpath, mode, value):
    """Get parameters settings for filtering hmmscan results. The results are
    filtered according to the HMM profile used and the type of signature
    protein associated with the HMM profile.

    :param yamlpath: YAML configuration file
    :param mode:     Which set of parameters should be extracted
                     (`params_usage` or `params_profiles`)
    :param value:    Name of HMM profile used or type of protein associated
                     with the HMM profile
    :type yamlpath:  :class:`str`
    :type mode:      :class:`str`
    :type value:     :class:`str`
    :return:         Set of parameters for filtering hmmscan results of a
                     given HMM profile
    :rtype:          :class:`dict`
    """
    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    params_set_name = params_dict["SP_detection"]["hmmscan"]

    # Search parameters set for prot_type
    try:
        params_set_name = params_set_name[mode]
    except KeyError:
        print(f'ERROR: Set of parameters for hmmscan are "params_usage" and '
              f'"params_profiles" in "{yamlpath}" config file. "{mode}" is '
              'not one of them!',
              file=sys.stderr)
        sys.exit(1)
    try:
        params_set_name = params_set_name[value]
    except KeyError:
        print(f'ERROR: There is no set of parameters for "{value}" in the '
              f'subsection "{mode}" of "hmmscan" of "{yamlpath}" '
              f'config file!',
              file=sys.stderr)
        sys.exit(1)

    try:
        params_set = params_dict["SP_detection"]["hmmscan"]
        params_set = params_set["params_filtering"][params_set_name]
    except KeyError:
        print(f'ERROR: The parameter set "{params_set_name}" is not defined '
              f'in the subsection "hmmscan" of "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    return(params_set)


def filter_df(params_set, header_ref, df):
    """Convert a set of rules for filtering hmmscan hits (`params_set`)
    into usable command for pandas.DataFrame filtering.
    Association between columns names and rules are described in `header_ref`.

    :param params_set: Set of rules for filtering hmmscan hits
    :param header_ref: Association table between rules and column names
    :param df:         Dataframe to filter out hmmscan false positives
    :type params_set:  :class:`dict`
    :type header_ref:  :class:`dict`
    :type df:          :class:`pandas.DataFrame`
    :return:           Filtered hmmscan hits
    :rtype:            :class:`pandas.DataFrame`
    """

    # By default set as possible SP, if fail any test set to "no"
    df["Possible_SP"] = "yes"

    for rule in list(params_set.keys()):
        # Decomposition of the rule
        modulation, field = breakdown_filtering_rule(rule, header_ref)

        cutoff_col = f"{field}_cutoff"
        validation_col = f"{field}_OK"
        df[cutoff_col] = float(params_set[rule])
        df[validation_col] = True

        # Check if field exist
        if field not in df.columns:
            print(f'ERROR: Column "{field}" '
                  f'does not exist in hmmscan output!', file=sys.stderr)
            sys.exit(1)
        # Apply filter
        elif modulation == "min":
            # Valid values must be superior or equal to cutoff
            slicer = df[field] < float(params_set[rule])
            df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            df.loc[list(df[slicer].index), validation_col] = False
        elif modulation == "max":
            # Valid values must be inferior or equal to cutoff
            slicer = df[field] > float(params_set[rule])
            df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            df.loc[list(df[slicer].index), validation_col] = False
        # Modulation is incorrect
        else:
            print('ERROR: Incorrect rule for filtering. '
                  'You should not see this warning!', file=sys.stderr)
            sys.exit(1)

    return(df)


def get_protein_type_hash(yamlpath):
    """Get Protein type based on which HMM profile was used.

    :param yamlpath: YAML configuration file
    :type yamlpath:  :class:`str`
    :return:         Sanitized name of the protein type
    :rtype:          :class:`dict`
    """

    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    try:
        protein_type = params_dict["SP_detection"]["hmmscan"]
        protein_type = protein_type["params_SP_type"]
    except KeyError:
        print(f'ERROR: params_SP_type is missing in "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    return(protein_type)


def pretty_df(df):
    """Order signature proteins by position on the genome and sort each result
    from best to worst. Sanitize string values if needed and format numeric
    values. Rename columns to more explicit terms.

    :param df:       Dataframe to work with
    :type df:        :class:`pandas.DataFrame`
    :return:         Dataframe prettified
    :rtype:          :class:`pandas.DataFrame`
    """

    # Order SP by position on the genome and
    # sort each result from best to worst
    df = df.sort_values(by=["CDS_num", "i-Evalue"],
                        ascending=[True, True])

    # Sanitize annotations (remove leading and trailing whitespace)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Format numeric columns
    df['E-value'] = df['E-value'].map(lambda x: '{0:.3g}'.format(x))
    df['i-Evalue'] = df['i-Evalue'].map(lambda x: '{0:.3g}'.format(x))
    df['HMM_coverage'] = df['HMM_coverage'].map(lambda x: '{0:.2f}'.format(x))
    df['CDS_coverage'] = df['CDS_coverage'].map(lambda x: '{0:.2f}'.format(x))
    for col in df.columns[df.columns.str.endswith("_cutoff")]:
        df[col] = df[col].map(lambda x: '{0:.1f}'.format(x))

    # Rename columns
    cols = {"target_name": "Profile_ID",
            "tlen": "Profile_length",
            "query_name": "#ICEscreen_ID",
            "qlen": "CDS_length",
            "E-value": "E-value_hmm",
            "i-Evalue": "i-Evalue_hmm",
            "domain_score": "Score_hmm",
            "domain_bias": "Bias_hmm",
            "seq_score": "Global_score",
            "seq_bias": "Global_bias",
            "description_of_target": "Profile_description"}

    df = df.rename(columns=cols)

    return(df)


def reorder_columns(df):
    """Reorder columns.

    :param df:       Dataframe to work with
    :type df:        :class:`pandas.DataFrame`
    :return:         Dataframe with columns reordered
    :rtype:          :class:`pandas.DataFrame`
    """

    # Reorder columns
    cols_common = ['#ICEscreen_ID', 'CDS_num', 'CDS', 'CDS_strand',
                   'CDS_start', 'CDS_end', 'CDS_length', 'CDS_Protein_type',
                   'Profile_description', 'Profile_ID', 'Profile_name',
                   'Profile_length', 'Ali_len_HMM', 'Ali_len_CDS',
                   'hmm_coord_from', 'hmm_coord_to', 'ali_coord_from',
                   'ali_coord_to', 'i-Evalue_hmm', 'E-value_hmm', 'Score_hmm',
                   'Bias_hmm', '#_domain', 'of_domain', 'Global_score',
                   'Global_bias', 'HMM_coverage', 'CDS_coverage',
                   'Possible_SP']

    cols_validation = list(df.columns)

    for col in cols_common:
        cols_validation.remove(col)

    # If no results
    if len(df.index) == 0:
        df = pd.DataFrame(columns=cols_common + cols_validation)
    else:
        df = df[cols_common + cols_validation]

    return(df)


# %%
if __name__ == "__main__":

    # Parse script arguments
    tsvpath, conffile, outall, outfiltered, outbest = parse_arguments()

    # Get hmmscan results
    data = get_hmmscan_df(tsvpath)

    # If no results
    if len(data.index) == 0:
        cols = ["#ICEscreen_ID", "CDS_num", "CDS", "CDS_strand", "CDS_start",
                "CDS_end", "CDS_length", "CDS_Protein_type",
                "Profile_description", "Profile_ID", "Profile_name",
                "Profile_length", "Ali_len_HMM", "Ali_len_CDS",
                "hmm_coord_from", "hmm_coord_to", "ali_coord_from",
                "ali_coord_to", "i-Evalue_hmm", "E-value_hmm", "Score_hmm",
                "Bias_hmm", "#_domain", "of_domain", "Global_score",
                "Global_bias", "HMM_coverage", "CDS_coverage", "Possible_SP",
                "i-Evalue_cutoff", "i-Evalue_OK", "CDS_coverage_cutoff",
                "CDS_coverage_OK", "HMM_coverage_cutoff", "HMM_coverage_OK"]
        data = pd.DataFrame(columns=cols)

        data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")

        data = data.drop(columns=["Possible_SP"])
        data.to_csv(outfiltered, index=False, sep="\t", decimal=".",
                    na_rep="NA")
        data.to_csv(outbest, index=False, sep="\t", decimal=".", na_rep="NA")

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

    for group in data.groupby(by=["prot_type", "Profile_name"]):
        # To filter out non significant results need to know:
        # "prot_type"
        prot_type = group[1]["prot_type"].iloc[0]
        # To get rules for filtering should have
        # and hmmprofile used
        hmmprofile = group[1]["Profile_name"].iloc[0]

        # Get parameters for filtering
        params_set = get_params_set_hmmscan(conffile,
                                            "params_usage",
                                            prot_type)
        params_set.update(get_params_set_hmmscan(conffile,
                                                 "params_profiles",
                                                 hmmprofile))

        # Add results of filtering for each group
        results.append(filter_df(params_set, header_ref, group[1].copy()))

    # Concatenate all results
    data = pd.concat(results, ignore_index=True)

    # Get protein type
    prot_type = get_protein_type_hash(conffile)
    data["CDS_Protein_type"] = [prot_type[x] for x in data["prot_type"]]
    # Remove "prot_type" column
    data = data.loc[:, data.columns != "prot_type"]

    data = pretty_df(data)

    # All results without filtering
    data = reorder_columns(data)
    data = data.sort_values(by="CDS_num", ascending=False)
    data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")

    # All validated results
    data = data[data["Possible_SP"] == "yes"]
    data.loc[:, data.columns != "Possible_SP"]\
        .to_csv(outfiltered, index=False, sep="\t", decimal=".", na_rep="NA")

    # Pick best result for each SP
    data = data[~data.CDS_num.duplicated()]
    data = reorder_columns(data)
    data.loc[:, data.columns != "Possible_SP"]\
        .to_csv(outbest, index=False, sep="\t", decimal=".", na_rep="NA")
