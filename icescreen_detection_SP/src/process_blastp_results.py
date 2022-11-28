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

import os
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
import sqlite3
from process_SP_hits import split_identifier, breakdown_filtering_rule


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
                          help="ICEscreen BlastP results (.tsv)",
                          required=True)
    required.add_argument('-c', '--config',
                          help="Path to YAML configuration file",
                          required=True)
    required.add_argument('-d', '--db',
                          help="Path to ICEscreen annotation database",
                          required=True)
    required.add_argument('--outall',
                          help="All BlastP results with annotation (.tsv)",
                          required=True),
    required.add_argument('--outfiltered',
                          help="All possible Signature Proteins identified"
                               " (.tsv)",
                          required=True)
    required.add_argument('--outbest',
                          help="Best Signature Proteins-like identified"
                          " (.tsv)",
                          required=True)

    # Group of optional arguments
    optional = parser.add_argument_group('optional argument')
    optional.add_argument('--detailed', dest='detailed_mode',
                          action='store_true',
                          default=False,
                          help="Fetch detailed ICEscreen annotation",
                          required=False)

    # Parse arguments
    args = parser.parse_args()

    return(args.input, args.config, args.db, args.outall, args.outfiltered,
           args.outbest, args.detailed_mode)


def get_blastp_df(tsvpath):
    """Parse BlastP results from TSV file and convert into
    :class:`pandas.DataFrame`.

    :param tsvpath: TSV file path
    :type tsvpath:  :class:`str`
    :return:        TSV information formatted as dataframe
    :rtype:         :class:`pandas.DataFrame`

    TSV file must have these 14 following columns:

    1.  sseqid
          subject sequence identifier
    2.  qseqid
          query sequence identifier
    3.  pident
          alignment identity percentage
    4.  length
          alignment length
    5.  mismatch
          number of mismatch in alignment
    6.  gapopen
          number of gap open
    7.  qstart
          start position of alignment in query sequence
    8.  qend
          end position of alignment in query sequence
    9.  sstart
          start position of alignment in subject sequence
    10. send
          end position of alignment in subject sequence
    11. evalue
          e-value of hit
    12. bitscore
          bitscore of hit
    13. qlen
          length of query sequence
    14. slen
          lenght of subject sequence

    Three values are computed based on BlastP results:

    qcovs
      query protein coverage
    scovs
      subject protein (hit) coverage
    qlen/slen
      ratio of query protein to subject protein length
    """

    columns_type = {"qseqid": str,
                    "sseqid": str,
                    "pident": np.float32,
                    "length": np.uint32,
                    "mismatch": np.uint32,
                    "gapopen": np.uint32,
                    "qstart": np.uint32,
                    "qend": np.uint32,
                    "sstart": np.uint32,
                    "send": np.uint32,
                    "evalue": np.float64, # "evalue": np.float32,
                    "bitscore": np.float32,
                    "qlen": np.uint32,
                    "slen": np.uint32}

    try:
        data = pd.read_csv(tsvpath,
                           sep="\t",
                           dtype=columns_type,
                           header=0,
                           # Comment if want to allow additional columns
                           names=columns_type.keys()
                           )
    except ValueError:
        data = pd.read_csv(tsvpath, sep="\t")
    except KeyError:
        print('ERROR: There is a problem in TSV file header. '
              '(Missing / Misspelled / Additional columns).',
              file=sys.stderr)
        sys.exit(1)

    # If no results
    if len(data.index) == 0:
        return(data)

    # Split ICEscreen CDS identifier into multiple columns
    #columns = ["CDS_num", "CDS", "CDS_strand", "CDS_start", "CDS_end"]
    columns = ["CDS_num", "CDS_locus_tag", "CDS_protein_id", "Genome_accession", "Genome_accession_rank", "CDS_strand", "CDS_start", "CDS_end"]
    cds_identifiers = np.vectorize(split_identifier)(data["qseqid"])
    cds_identifiers = np.transpose(cds_identifiers)
    data[columns] = pd.DataFrame(cds_identifiers, columns=columns)

    # Compute coverages
    data["qcovs"] = (data["qend"] - data["qstart"] + 1) / data["qlen"] * 100
    data["scovs"] = (data["send"] - data["sstart"] + 1) / data["slen"] * 100

    # Compute lengths ratio
    data["qlen/slen"] = data["qlen"] / data["slen"]

    return(data)


def get_params_set_blastp(yamlpath, mode, dbname):
    """Get parameters settings for filtering/validate BlastP results by type of
    protein.

    :param yamlpath: YAML configuration file
    :param mode:     Which set of parameters should be extracted
                     (`params_filtering` or `params_validation`)
    :param dbname:   BlastP database used
    :type yamlpath:  :class:`str`
    :type mode:      :class:`str`
    :type dbname:    :class:`str`
    :return:         Set of parameters for filtering or validate (`mode`)
                     BlastP results of a given type of protein (`dbname`)
    :rtype:          :class:`dict`
    """
    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    try:
        params_set_name = params_dict["SP_detection"]["blastp"]
        params_set_name = params_set_name["params_usage"][dbname]
    except KeyError:
        print(f'ERROR: There is no set of parameters for the BLAST database '
              f'"{dbname}" in "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    try:
        params_set = params_dict["SP_detection"]["blastp"]
        params_set = params_set[mode][params_set_name]
    except KeyError:
        print(f'ERROR: The parameter set "{params_set_name}" is not defined '
              f'in "{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    return(params_set)


def filter_df(params_set, header_ref, df):
    """Convert a set of rules for filtering BlastP hits (`params_set`)
    into usable command for pandas.DataFrame filtering.
    Association between columns names and rules are described in `header_ref`.

    :param params_set: Set of rules for filtering BlastP hits
    :param header_ref: Association table between rules and column names
    :param df:         Dataframe to filter out BlastP false positives
    :type params_set:  :class:`dict`
    :type header_ref:  :class:`dict`
    :type df:          :class:`pandas.DataFrame`
    :return:           Nothing, `df` is modified in place
    """

    # By default set as possible SP, if fail any test set to "no"
    df["Possible_SP"] = "yes"

    for rule in list(params_set.keys()):
        # Decomposition of the rule
        modulation, field = breakdown_filtering_rule(rule, header_ref)

        if isinstance(field, tuple):
            # Check if field exist
            for f in field:
                if f not in df.columns:
                    print(f'ERROR: Column "{f}" '
                          f'does not exist in BlastP output!', file=sys.stderr)
                    sys.exit(1)
            # Apply filter
            if modulation == "min":
                # If valid: First value must be superior or equal to a given
                # percentage of the second value
                slicer = df[field[0]] < \
                         (df[field[1]] * float(params_set[rule]) / 100)
                df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            elif modulation == "max":
                # If valid: First value must be inferior or equal to a given
                # percentage of the second value
                slicer = df[field[0]] > \
                         (df[field[1]] * float(params_set[rule]) / 100)
                df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            # Modulation is incorrect
            else:
                print('ERROR: Incorrect rule for filtering. '
                      'You should not see this warning!', file=sys.stderr)
                sys.exit(1)
        else:
            # Check if field exist
            if field not in df.columns:
                print(f'ERROR: Column "{field}" '
                      f'does not exist in BlastP output!', file=sys.stderr)
                sys.exit(1)
            # Apply filter
            elif modulation == "min":
                # Valid values must be superior or equal to cutoff
                slicer = df[field] < float(params_set[rule])
                df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            elif modulation == "max":
                # Valid values must be inferior or equal to cutoff
                slicer = df[field] > float(params_set[rule])
                df.loc[list(df[slicer].index), "Possible_SP"] = "no"
            # Modulation is incorrect
            else:
                print('ERROR: Incorrect rule for filtering. '
                      'You should not see this warning!', file=sys.stderr)
                sys.exit(1)


def get_protein_type(yamlpath, dbname):
    """Get Protein type based on which BlastP database was used.

    :param yamlpath: YAML configuration file
    :param dbname:   BlastP database used
    :type yamlpath:  :class:`str`
    :type dbname:    :class:`str`
    :return:         Protein type of proteins in BlastP database used
    :rtype:          :class:`str`
    """

    with open(yamlpath, 'r') as filin:
        params_dict = yaml.safe_load(filin)

    try:
        protein_type = params_dict["SP_detection"]["blastp"]
        protein_type = protein_type["params_SP_type"][dbname]
    except KeyError:
        print(f'ERROR: There is no protein type set for "{dbname}" in '
              f'"{yamlpath}" config file!',
              file=sys.stderr)
        sys.exit(1)

    return(protein_type)


def get_annotation(dbpath, dbtable, data):
    """Fetch annotation of query protein in BlastP ICEscreen
    results from ICEscreen annotation database. Then add annotation to
    given ICEscreen result.

    :param dbpath:  ICEscreen database path
    :param dbtable: Name of table to fetch annotation from
    :param data:    ICEscreen BlastP results
    :type dbpath:   :class:`str`
    :type dbtable:  :class:`str`
    :type data:     :class:`pandas.DataFrame`
    :return:        ICEscreen BlastP results annotated
    :rtype:         :class:`pandas.DataFrame`
    """

    # Connect to annotation database
    db = sqlite3.connect(dbpath)

    # Initalization of cursor
    dbCursor = db.cursor()

    if dbtable == "Firmicutes_SP_detailed":
        # Dictionary with:
        # keys: Columns of ICEscreen annotation database table to be retrieved
        # values: Associated columns names in dataframe for retrieved info
        cols = {"Prot_ID": "sseqid",
                "Prot_Type": "Protein_type_of_blast_most_similar_ref_SP",
                "Associated_element_DynAMic_name": "Associated_element",
                "Associated_element_type": "Associated_element_type_of_blast_most_similar_ref_SP",
                "Genome_ID": "Genome_Query",
                "Associated_integrase": "Associated_integrase",
                "integrase_position_vs_integration_site": "Int_position_vs"
                                                          "_insertion_gene",
                "complement_vs_": "Int_strand_vs_insertion_gene_strand",
                "Insertion_site": "Query_Insertion_site",
                "Element_type_family": "Query_Element_type_family",
                "ICE_superfamily": "ICE_superfamily_of_most_similar_ref_SP",
                "ICE_family": "ICE_family_of_most_similar_ref_SP",
                "IME_superfamily": "IME_superfamily_of_most_similar_ref_SP",
                "relaxase_family_domain": "Relaxase_family_domain_of_most_similar_ref_SP",
                "relaxase_family_MOB": "Relaxase_family_MOB_of_most_similar_ref_SP",
                "coupling_type": "Coupling_type_of_most_similar_ref_SP",
                "False_positives": "False_positives"
                }

    elif dbtable == "Firmicutes_SP":
        # Dictionary with:
        # keys: Columns of ICEscreen annotation database table to be retrieved
        # values: Associated columns names in dataframe for retrieved info
        cols = {"Prot_ID": "sseqid",
                "Prot_Type": "Protein_type_of_blast_most_similar_ref_SP",
                "Associated_element_type": "Associated_element_type_of_blast_most_similar_ref_SP",
                "ICE_superfamily": "ICE_superfamily_of_most_similar_ref_SP",
                "ICE_family": "ICE_family_of_most_similar_ref_SP",
                "IME_superfamily": "IME_superfamily_of_most_similar_ref_SP",
                "relaxase_family_domain": "Relaxase_family_domain_of_most_similar_ref_SP",
                "relaxase_family_MOB": "Relaxase_family_MOB_of_most_similar_ref_SP",
                "coupling_type": "Coupling_type_of_most_similar_ref_SP",
                "False_positives": "False_positives"
                }

    # SQL command to fetch columns in dbtable table for given IDs
    dbCursor.execute(f"""SELECT {", ".join(f'"{c}"' for c in cols.keys())} """
                     f"""FROM {dbtable} WHERE Prot_id IN """
                     f"""('{"', '".join(data['sseqid'].unique())}')""")

    fetched = dbCursor.fetchall()

    # Formatting fetched information into pandas.DataFrame
    info = pd.DataFrame.from_records(fetched, columns=list(cols.values()))

    # Cast ID column into str
    info['sseqid'] = info['sseqid'].astype(str)
    data['sseqid'] = data['sseqid'].astype(str)

    # Add fetched annotation into data
    return(data.merge(info, on="sseqid", how="left"))


def breakdown_validation_rule(rule, header_ref):
    """Break down a rule into two information.

    1. modulation: Can be "gt", "ge", "eq", "lt" or "le"
    2. field: Can be one field only

    :param rule:       Rule to break down
    :param header_ref: Association table between rules and columns names
    :type rule:        :class:`str`
    :return:           Rule's modulation (:class:`str`) and fields names
                       (:class:`str` or :class:`tuple`)
    :rtype:            :class:`tuple`
    """

    try:
        modulation, field_name = rule.split("_", maxsplit=1)
    except ValueError:
        print("ERROR: Validation rule for BlastP hit "
              "must have this format!: [gt/ge/eq/lt/le]_[parameter name]",
              file=sys.stderr)
        sys.exit(1)

    if modulation not in ["gt", "ge", "eq", "lt", "le"]:
        print('ERROR: Validation rule modulator must be one of these: '
              '"gt", "ge", "eq", "lt", "le"!', file=sys.stderr)
        sys.exit(1)

    try:
        field = header_ref[field_name]
    except KeyError:
        print(f'ERROR: Parameter "{field_name}" does not exist '
              f'in the header references!', file=sys.stderr)
        sys.exit(1)

    return(modulation, field)


def result_validation(variable, value, modulation):
    """Test if result is valid.

    :param variable:   Which variable should be tested
    :param value:      Cutoff for validation
    :param modulation: Which test should be performed for validation
    :type variable:    :class:`float`
    :type value:       :class:`float`
    :type modulation:  :class:`str`
    :return:           Result of test
    :rtype:            :class:`bool`
    """
    if modulation == "gt":
        if variable > float(value):
            return(True)
        else:
            return(False)
    if modulation == "ge":
        if variable >= float(value):
            return(True)
        else:
            return(False)
    if modulation == "eq":
        if variable == float(value):
            return(True)
        else:
            return(False)
    if modulation == "lt":
        if variable < float(value):
            return(True)
        else:
            return(False)
    if modulation == "le":
        if variable <= float(value):
            return(True)
        else:
            return(False)
    else:
        print("ERROR: Incorrect rule for validation. "
              "You should not see this warning!", file=sys.stderr)
        sys.exit(1)


def get_SP_validation(params_set, header_ref, result):
    """Evaluate the annotation validity of each possible signature protein.
    The set of rules for validating BlastP hits is in `params_set`.
    Association between columns names and rules are described in `header_ref`.

    :param params_set: Set of rules for validating BlastP hits
    :param header_ref: Association table between rules and column names
    :param result:     BlastP hit to evaluate
    :type params_set:  :class:`dict`
    :type header_ref:  :class:`dict`
    :type result:      :class:`pandas.Series`
    :return:           The result of the evaluation
    :rtype:            :class:`str`
    """

    # Find the status of the annotation of the possible protein signature
    for status in list(params_set.keys()):
        # Initialize list to get result of each test
        validity = []

        for rule in list(params_set[status].keys()):
            # Decomposition of the rule
            modulation, field = breakdown_validation_rule(rule, header_ref)

            # Check if field exist
            if not isinstance(field, str):
                print(f'ERROR: Column "{field}" does not exist '
                      f'in BlastP output!', file=sys.stderr)
                sys.exit(1)
            else:
                if field not in result.index:
                    print(f'ERROR: Column "{field}" does not exist '
                          f'in BlastP output!', file=sys.stderr)
                    sys.exit(1)
                else:
                    # Carry out the test and add result to list
                    validation = result_validation(float(result[field]),
                                                   float(params_set[status]
                                                   [rule]),
                                                   modulation)
                    validity.append(validation)

        # Return validation status if the annotation passed all tests
        if all(validity):
            return(status)

    # If failed all tests for all validation status
    return("No validation")


def get_annotation_validation(df, column):
    """Check if ICEscreen annotation should be use or not based on identity
    percentage.

    :param df:     Signature proteins found by ICEscreen by BlastP
    :param column: Name of column to stock validation status
    :type df:      :class:`pandas.DataFrame`
    :type column:  :class:`str`
    :return:       Nothing, modify in place
    """

    df[column] = "no"
    mask = df["pident"] >= 30
    df.loc[mask, column] = "yes"


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
    df = df.groupby("CDS_num", as_index=False, group_keys=True)\
           .apply(lambda x: x.sort_values(["evalue", "bitscore", "qcovs",
                                           "scovs", "pident"],
                                          ascending=[True, False, False,
                                                     False, False]))

    # Sanitize annotations (remove leading and trailing whitespace)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    # Format numeric columns
    df['pident'] = df['pident'].map(lambda x: '{0:.2f}'.format(x))
    df['evalue'] = df['evalue'].map(lambda x: '{0:.3g}'.format(x))
    df['bitscore'] = df['bitscore'].map(lambda x: '{0:.1f}'.format(x))
    df['qcovs'] = df['qcovs'].map(lambda x: '{0:.2f}'.format(x))
    df['scovs'] = df['scovs'].map(lambda x: '{0:.2f}'.format(x))
    df['qlen/slen'] = df['qlen/slen'].map(lambda x: '{0:.2f}'.format(x))

    # Rename columns
    cols = {"qseqid": "#ICEscreen_ID",                     # ICEscreen_ID
            "sseqid": "Id_of_blast_most_similar_ref_SP",                       # qblast
            "pident": "Blast_ali_identity_perc",                 # ali_pident
            "length": "Blast_ali_length",                        # ali_len
            "mismatch": "Mismatch",                        # mismatch
            "gapopen": "Gapopen",                          # gapopen
            "qstart": "Blast_ali_start_CDS",                     # ali_start_cds
            "qend": "Blast_ali_end_CDS",                         # ali_end_cds
            "sstart": "Blast_ali_start_Query_blast",             # ali_start_qblast
            "send": "Blast_ali_end_Query_blast",                 # ali_end_qblast
            "evalue": "Blast_ali_E-value",                     # evalue_blast
            "bitscore": "Blast_ali_bitscore",                  # bitscore_blast
            "qlen": "CDS_length",                          # cds_len
            "slen": "Length_of_blast_most_similar_ref_SP",                  # qblast_len
            "qcovs": "CDS_coverage",                       # cds_cov
            "scovs": "Blast_ali_coverage_most_similar_ref_SP",               # qblast_cov
            "qlen/slen": "CDS_length/Length_of_blast_most_similar_ref_SP"}  # cds_len/qblast_len

    df = df.rename(columns=cols)

    return(df)


def reorder_columns(df, detailed):
    """Reorder columns.

    :param df:       Dataframe to work with
    :param detailed: If annotation in dataframe is detailed or not
    :type df:        :class:`pandas.DataFrame`
    :type detailed:  :class:`bool`
    :return:         Dataframe with columns reordered
    :rtype:          :class:`pandas.DataFrame`
    """

    # Reorder columns
    cols_common = [["#ICEscreen_ID", "CDS_num",
                    #"CDS", # "CDS" has been changed to  "CDS_locus_tag", "CDS_protein_id", "Genome_accession", "Genome_accession_rank" -> change it everywhere
                    "Genome_accession", "Genome_accession_rank", "CDS_locus_tag", "CDS_protein_id",
                    "CDS_strand", "CDS_start", "CDS_end", "CDS_length", "CDS_Protein_type",
                    "Description_of_blast_most_similar_ref_SP",
                    "Id_of_blast_most_similar_ref_SP", "Length_of_blast_most_similar_ref_SP", "Blast_ali_length",
                    "Mismatch", "Gapopen", "Blast_ali_start_CDS", "Blast_ali_end_CDS",
                    "Blast_ali_start_Query_blast", "Blast_ali_end_Query_blast",
                    "Blast_ali_identity_perc", "Blast_ali_E-value", "Blast_ali_bitscore",
                    "CDS_coverage", "Blast_ali_coverage_most_similar_ref_SP",
                    "CDS_length/Length_of_blast_most_similar_ref_SP",
                    "Protein_type_of_blast_most_similar_ref_SP"],
                   ["ICE_superfamily_of_most_similar_ref_SP", "ICE_family_of_most_similar_ref_SP", "IME_superfamily_of_most_similar_ref_SP",
                    "Relaxase_family_domain_of_most_similar_ref_SP", "Relaxase_family_MOB_of_most_similar_ref_SP",
                    "Coupling_type_of_most_similar_ref_SP", "False_positives", "Possible_SP",
                    "SP_blast_validation", "Use_annotation"]]

    # Not detailed results columns
    cols_not_detailed = cols_common[0] + ["Associated_element_type_of_blast_most_similar_ref_SP"] +\
        cols_common[1]
    # Detailed results columns
    cols_detailed = cols_common[0] + ["Associated_element",
                                      "Associated_element_type_of_blast_most_similar_ref_SP",
                                      "Genome_Query",
                                      "Associated_integrase",
                                      "Int_position_vs_insertion_gene",
                                      "Int_strand_vs_insertion_gene_strand",
                                      "Query_Element_type_family"] +\
        cols_common[1]

    # If no results
    if len(df.index) == 0:
        if detailed is False:
            df = pd.DataFrame(columns=cols_not_detailed)
        elif detailed is True:
            df = pd.DataFrame(columns=cols_detailed)
    else:
        if detailed is False:
            df = df[cols_not_detailed]
        elif detailed is True:
            df = df[cols_detailed]

    return(df)

def remove_family_info_if_perc_identity_less_40(df):
    """Remove the family information for SP whose percentage identity is less than 40 percent.

    :param df:       dataframe of BlastP hits that might be a signature protein
    :return:         Nothing, modify in place
    """
    df['ICE_family_of_most_similar_ref_SP'] = np.where((df['pident'] < 40),
                           "-",
                           df['ICE_family_of_most_similar_ref_SP'])
                           

def remove_family_info_from_type_integrase(df):
    """Hide the family information for SP of type integrase.

    :param df:       dataframe of BlastP hits that might be a signature protein
    :return:         Nothing, modify in place
    """
    df['ICE_superfamily_of_most_similar_ref_SP'] = np.where((df['Protein_type_of_blast_most_similar_ref_SP'] == 'Integrase'),
                           "-",
                           df['ICE_superfamily_of_most_similar_ref_SP'])

    df['ICE_family_of_most_similar_ref_SP'] = np.where((df['Protein_type_of_blast_most_similar_ref_SP'] == 'Integrase'),
                           "-",
                           df['ICE_family_of_most_similar_ref_SP'])
    df['IME_superfamily_of_most_similar_ref_SP'] = np.where((df['Protein_type_of_blast_most_similar_ref_SP'] == 'Integrase'),
                           "-",
                           df['IME_superfamily_of_most_similar_ref_SP'])



def get_description(df, cds_type, col):
    """Generate a clearer description of BlastP hit result.

    :param df:       BlastP hits that might be a signature protein
    :param cds_type: Type of signature protein
    :param col:      Name of column for description
    :type df:        :class:`pandas.DataFrame`
    :type cds_type:  str
    :type col:       str
    :return:         Nothing, modify in place
    """

    # If signature protein type is any integrase or VirB4, the description is
    # the type of the protein
    df[col] = df["CDS_Protein_type"]

    # If the CDS is a false positive, the description is information about
    # the type of false positive
    df[col] = np.where(df['False_positives'] != '-',
                       df['False_positives'],
                       df[col])

    # If signature protein type is Coupling protein, the description is the
    # type of protein and superfamily/family of the protein
    if cds_type == "Coupling protein":
        df[col] = np.where((df['False_positives'] == '-') & (df["Coupling_type_of_most_similar_ref_SP"] != "-"),
                           df[col] + " " + df["Coupling_type_of_most_similar_ref_SP"],
                           df[col])
    # If signature protein type is Relaxase, the description is the
    # type of protein and superfamily/family of the protein
    elif cds_type == "Relaxase":
        df[col] = np.where((df['False_positives'] == '-') & (df["Relaxase_family_MOB_of_most_similar_ref_SP"] != "-") & (df["Relaxase_family_domain_of_most_similar_ref_SP"] != "-"),
                           df[col] + " " + df["Relaxase_family_MOB_of_most_similar_ref_SP"] + " (" + df["Relaxase_family_domain_of_most_similar_ref_SP"] + ")",
                           df[col])
        df[col] = np.where((df['False_positives'] == '-') & (df["Relaxase_family_MOB_of_most_similar_ref_SP"] != "-") & (df["Relaxase_family_domain_of_most_similar_ref_SP"] == "-"),
                           df[col] + " " + df["Relaxase_family_MOB_of_most_similar_ref_SP"],
                           df[col])
        df[col] = np.where((df['False_positives'] == '-') & (df["Relaxase_family_MOB_of_most_similar_ref_SP"] == "-") & (df["Relaxase_family_domain_of_most_similar_ref_SP"] != "-"),
                           df[col] + " " + df["Relaxase_family_domain_of_most_similar_ref_SP"],
                           df[col])


if __name__ == "__main__":
    # Parse script arguments
    tsvpath, conffile, annotdb, outall, outfiltered,\
        outbest, detailed = parse_arguments()

    # Get BlastP results
    data = get_blastp_df(tsvpath)

    # If no results
    if len(data.index) == 0:
        data = reorder_columns(data, detailed)
        data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")

        data = data.drop(columns=["Possible_SP"])

        data.to_csv(outfiltered, index=False, sep="\t", decimal=".",
                    na_rep="NA")
        data.to_csv(outbest, index=False, sep="\t", decimal=".", na_rep="NA")

        sys.exit(0)

    # Filter out non significant hits from dataframe
    # Get Blastp database name used for these results (in filename)
    dbname = os.path.splitext(os.path.basename(tsvpath))[0].split("_")[-1]

    # Get parameters for filtering BlastP hits
    params_set = get_params_set_blastp(conffile, "params_filtering", dbname)

    # Get dictionnary to make correspondance between fields in header and
    # rule name in configuration file
    header_ref = {"query_length": "qlen",
                  "percentage_identity": "pident",
                  "alignment_coverage_to_query": "qcovs",
                  "alignment_coverage_to_subject": "scovs",
                  "evalue": "evalue",
                  "length_ratio": "qlen/slen"}

    # Filter BlastP hits to get Signature Proteins (SPs)
    filter_df(params_set, header_ref, data)

    prot_type = get_protein_type(conffile, dbname)
    data["CDS_Protein_type"] = prot_type

    # Add ICEscreen annotation to each possible SPs
    if detailed is True:
        dbtable = "Firmicutes_SP_detailed"
    elif detailed is False:
        dbtable = "Firmicutes_SP"

    data = get_annotation(annotdb, dbtable, data)

    # Add description of each result
    get_description(data, prot_type, "Description_of_blast_most_similar_ref_SP")

    # Add validation of SP:
    # If the gene coding for the signature protein is functional or not
    # Get parameters for validating Signature Protein
    params_set = get_params_set_blastp(conffile, "params_validation", dbname)
    data["SP_blast_validation"] = data.apply(
        lambda x: get_SP_validation(params_set, header_ref, x),
        axis=1)

    # Add validation of SP enriched annotation:
    # If % identity is too low, the enriched annotation should not be used
    get_annotation_validation(data, "Use_annotation")

    # Do not take into account family for integrase
    remove_family_info_from_type_integrase(data)
    remove_family_info_if_perc_identity_less_40(data)

    # Save results
    data = pretty_df(data)

    # All results without filtering
    data = reorder_columns(data, detailed)

    data.to_csv(outall, index=False, sep="\t", decimal=".", na_rep="NA")

    # All validated results
    data = data[data["Possible_SP"] == "yes"]
    data.loc[:, data.columns != "Possible_SP"]\
        .to_csv(outfiltered, index=False, sep="\t", decimal=".", na_rep="NA")

    # Pick best result for each SP
    data = data[~data.CDS_num.duplicated()]

    data = reorder_columns(data, detailed)
    data.loc[:, data.columns != "Possible_SP"]\
        .to_csv(outbest, index=False, sep="\t", decimal=".", na_rep="NA")
