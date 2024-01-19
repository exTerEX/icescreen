#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
from Bio import SeqIO
import sys
import Bio.Data.CodonTable as CodonTable
import os
from pathlib import Path
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import yaml
import commonMethods


def translate_sequence(input_seq, translation_tableSend, isCDSTrueSend, locus_tagIT):
    """Wrapper for Biopython translate function.  Bio.Seq.translate will complain if input sequence is 
    not a mulitple of 3.  This wrapper function passes an acceptable input to Bio.Seq.translate in order to
    avoid this warning."""
    trailing_bases = len(input_seq) % 3
    # print("BEFORE : translate_sequence on locus_tagIT = {} ; len(input_seq) = {} ; trailing_bases = {} ; translation_tableSend = {} ; isCDSTrueSend = {}".format(locus_tagIT, str(len(input_seq)), str(trailing_bases), str(translation_tableSend), str(isCDSTrueSend)))
    if trailing_bases == 1:
        input_seq = input_seq + "NN"
        # input_seq = ''.join([input_seq, 'NN'])
    if trailing_bases == 2:
        input_seq = input_seq + "N"
        # input_seq = ''.join([input_seq, 'N'])
    # print("AFTER : translate_sequence on locus_tagIT = {} ; len(input_seq) = {} ; trailing_bases = {}".format(locus_tagIT, str(len(input_seq)), str(len(input_seq) % 3)))
    # if trailing_bases:
    #     input_seq = ''.join([input_seq, 'NN']) if trailing_bases == 1 else ''.join([input_seq, 'N'])

    output_seq = Seq.translate(input_seq, table=translation_tableSend, cds=isCDSTrueSend)

    if trailing_bases == 1 or trailing_bases == 2 :
        #remove last residue if input needed to be extended because of trailing bases
        output_seq = output_seq[:-1]

    return output_seq


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



def gb_to_faa(record, sequence_type, feature_type, codon_table_if_unspecified_in_input_file, counter, dictLocusTags, dictGenomeIds, in_file, genome_accession_rank_sent):
    
    """Extract CDS from genbank file and save into target multifasta file.

    ===================================
    Sequence identifier (#ICEscreen_ID)
    ===================================
    An identifier is generated for each CDS with 6 information:

    Format
    ------
     ``[CDS_number]&locus_tag=[locus_tag]&protein_id=[protein_id]&genome_accession=[ACCESSION]&genome_accession_rank=[rank]|[Strand]|[Start]..[End]|[Pseudo]``

    CDS number
      The number of the CDS in the genome
    Locus_tag
      An uniq identifier of the gene
    Protein_ID
      An identifier of the protein.
    Genome_accession
      An identifier of the genome record
    Genome_accession_rank
      The rank of the genome record in the genbank file
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

    genome_accession = record.id
    if len(genome_accession) == 0:
        raise RuntimeError("Error in gb_to_faa: empty genome_accession in file {}".format(
            str(in_file)))
    if genome_accession in dictGenomeIds :
        raise RuntimeError("Error in gb_to_faa: duplicate genome_accession {} in file {}".format(
            str(genome_accession),
            str(in_file)))
    else :
        dictGenomeIds[genome_accession] = 1

    # First order features by location on genome
    record.features.sort(key=lambda x: int(x.location.start))

    new_records = []
    # Parsing sequences with SeqIO
    for feature in record.features:
        # Get location (start, end and strand)
        # start + 1 to transform index into location on the chromosome
        start = feature.location.start + 1
        end = feature.location.end
        strand = feature.location.strand

        # Initialization of temp_seq and protId variables
        temp_record = ''
        temp_seq = ''
        locus_tagIT = ""
        protein_idIT = ""

        # Format information of location
        if strand == 1:
            strandInfo = "+|{:d}..{:d}".format(start, end)
        if strand == -1:
            strandInfo = "-|{:d}..{:d}".format(start, end)

        if feature.type == feature_type: # What kind of feature to extract. Usually CDS or tRFLP
            # If pseudogene might have "pseudo" or "pseudogene" qualifier
            if any(x in feature.qualifiers for x in ["pseudo", "pseudogene"]):
                pseudo = "|Pseudo"
            else:
                pseudo = ""

                
            # Search protein ID in 4 qualifiers
            if 'protein_id' in feature.qualifiers:
                protein_idIT = commonMethods.search_feature('protein_id', feature)
            # elif 'label' in feature.qualifiers:
            #     protId = search_feature('label', feature)
            # elif 'locus_tag' in feature.qualifiers:
            #     protId = search_feature('locus_tag', feature)
            # elif 'db_xref' in feature.qualifiers:
            #     protId = search_feature('db_xref', feature)
            if 'locus_tag' in feature.qualifiers:
                locus_tagIT = commonMethods.search_feature('locus_tag', feature)
                if len(locus_tagIT) != 0:
                    if locus_tagIT in dictLocusTags:
                        raise RuntimeError("Error in gb_to_faa: duplicate locus_tag {} in file {}".format(
                            str(locus_tagIT),
                            str(in_file)))
                    else :
                        dictLocusTags[locus_tagIT] = 1


            if sequence_type == 'nt':
                #header=f'{counter}_{protId}|{strandInfo}'
                header=f'{counter}&locus_tag={locus_tagIT}&protein_id={protein_idIT}&genome_accession={record.id}&genome_accession_rank={genome_accession_rank_sent}|{strandInfo}'
                temp_record = SeqRecord(feature.extract(record.seq), id = header,\
                description = '', annotations={"molecule_type": "DNA"})
            elif sequence_type == 'taa':
                translation_table = None
                if "transl_table" in feature.qualifiers:
                    translation_table = feature.qualifiers["transl_table"][0]
                else:
                    translation_table = codon_table_if_unspecified_in_input_file
                #temp_record = SeqRecord(feature.extract(record.seq).translate(table = translation_table),\
                #id = header)
                nucSeq = feature.location.extract(record).seq
                #nucSeq = pad_seq(nucSeq_no_pad)
                # Check if CDS or pseudogene
                try:
                    isCDSTrue = False
                    if feature_type == "CDS" :
                        isCDSTrue = True
                    temp_seq = translate_sequence(nucSeq, translation_table, isCDSTrue, locus_tagIT) # was Seq.translate
                    # temp_seq = nucSeq.translate(table=translation_table, cds=isCDSTrue)
                    # temp_seq = Seq.translate(nucSeq, table=translation_table, cds=isCDSTrue)
                # If it's an Translation Error then it's a pseudogene
                except Exception as e:
                    if isinstance(e, CodonTable.TranslationError):
                        pseudo = "|Pseudo"
                        # Translate sequence without CDS parameter
                        temp_seq = translate_sequence(nucSeq, translation_table, False, locus_tagIT)
                        # temp_seq = nucSeq.translate(table=translation_table, cds=False)
                        # temp_seq = Seq.translate(nucSeq, table=translation_table, cds=False)

                #header=f'{counter}_{protId}|{strandInfo}{pseudo}'
                header=f'{counter}&locus_tag={locus_tagIT}&protein_id={protein_idIT}&genome_accession={record.id}&genome_accession_rank={genome_accession_rank_sent}|{strandInfo}{pseudo}'
                temp_record = SeqRecord(temp_seq, id = header, \
                description = '', annotations={"molecule_type": "protein"})
            elif sequence_type == 'aa':
                if "translation" in feature.qualifiers:
                    temp_seq = Seq(feature.qualifiers["translation"][0]) #, IUPAC.protein
                else:           
                    translation_table = None
                    if "transl_table" in feature.qualifiers:
                        translation_table = feature.qualifiers["transl_table"][0]
                    else:
                        translation_table = codon_table_if_unspecified_in_input_file
                    nucSeq = feature.location.extract(record).seq
                    # Check if CDS or pseudogene
                    try:
                        isCDSTrue = False
                        if feature_type == "CDS" :
                            isCDSTrue = True
                        temp_seq = translate_sequence(nucSeq, translation_table, isCDSTrue, locus_tagIT)
                        # temp_seq = nucSeq.translate(table=translation_table, cds=isCDSTrue)
                        # temp_seq = Seq.translate(nucSeq, table=translation_table, cds=isCDSTrue)
                    # If it's an Translation Error then it's a pseudogene
                    except Exception as e:
                        if isinstance(e, CodonTable.TranslationError):
                            pseudo = "|Pseudo"
                            # Translate sequence without CDS parameter
                            temp_seq = translate_sequence(nucSeq, translation_table, False, locus_tagIT)
                            # temp_seq = nucSeq.translate(table=translation_table, cds=False)
                            # temp_seq = Seq.translate(nucSeq, table=translation_table, cds=False)
                    

                #header=f'{counter}_{protId}|{strandInfo}{pseudo}'
                header=f'{counter}&locus_tag={locus_tagIT}&protein_id={protein_idIT}&genome_accession={record.id}&genome_accession_rank={genome_accession_rank_sent}|{strandInfo}{pseudo}'
                temp_record = SeqRecord(temp_seq, id = header, \
                description = '', annotations={"molecule_type": "protein"})
                
            if len(str(temp_seq)) > 0:
                new_records.append(temp_record)
                # Increment counter
                counter += 1
    return (new_records, counter, dictLocusTags, dictGenomeIds)


# def gb_to_faa_whole(record, annotation_list):
#     header = []
#     for item in annotation_list:
#         if item in record.annotations:
#             header_part = record.annotations[item]
#             if type(header_part) == type([]): #Some attributes are lists. Must turn into string
#                 header_part = ' : '.join(header_part)
#             header_part = header_part.replace("\n"," ") #Catch improper newline character
#             header_part = header_part.replace("                      ", " ")#Catch inproper spaces
#             header.append(header_part)
#         else:
#             header_part = 'missing_%s_annotation' % item
#             header.append(header_part)
#     header = "|".join(header)
#     temp_record = SeqRecord(record.seq, id = header, description = '')
#     return [temp_record] #Return a list because will be used in a context requiring a list


##############
#### Main ####
##############
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Create multifasta protein "
                                                 "from a genbank file")
    # Group of mandatory arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--in_file', help="Specify the input file that you wish to convert", required=True)
    required.add_argument('-o', '--out_file', help="Specify the path and name of the output fasta file you wish to create.", required=True)
    facultative = parser.add_argument_group('facultative arguments')
    facultative.add_argument('-m', '--file_format', help="Specify the input file format. Specify 'genbank' or 'embl'. Default is genbank.", required=False)
    facultative.add_argument('-d', '--codon_table_if_unspecified_in_input_file', help="codon table must be one of the following -> https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. Default is 11.", required=False)
    facultative.add_argument('-force_delete_output_file', help="Flag to force delete the output fasta file ; default=False", action='store_true')
    facultative.add_argument('-s', '--sequence_type', help="Specify the kind of sequence you would like to extract. Options are 'aa' "
                "(feature amino acids), 'nt' (feature nucleotides), 'whole' (the entire "
                "sequence, not just sequence corresponding to features) and 'taa'    (amino acids "
                "translated on the fly, which generates amino acid sequence by translating the "
                "nucleotide sequence rather than extracting from the feature table)."
                "Default is 'aa'.", required=False)
    facultative.add_argument('-f', '--feature_type', help="Specify the type of feature that you would like to extract. This option "
                "accepts arbitrary text, and will fail if you input a non-existent feature name. "
                "Common options are 'CDS', 'rRNA', 'tRNA', or 'gene'. Default is 'CDS'.", required=False)
    facultative.add_argument('-c', '--config',
                          help="Path to YAML configuration file. An example file is Projects/icescreen/pipelines/mode/bacillota.yml",
                          required=False)
    
    args = parser.parse_args()
    
    # see gb_to_faa.py from the BYG_ASS_CAT project for a version with directory as input and conversion of input files recursively
    in_file = args.in_file
    in_file = os.path.abspath(in_file)
    out_file = args.out_file
    out_file = os.path.abspath(out_file)
    #print(f"{in_file.name}:\n{in_file.read_text()}\n")
    #print("File '%s'" % in_file)
    (in_filePath, in_fileWholeName) = os.path.split(in_file)
    #relativePathAfterInDir = Path(in_filePath).relative_to(in_dir)
    (in_fileBase, in_fileExt) = os.path.splitext(in_fileWholeName)
    #Figure out what our out_file is.
    #out_file = os.path.join(out_dir, relativePathAfterInDir, in_fileBase + '.fasta')
    #out_file = os.path.abspath(out_file)
    (out_filePath, out_fileWholeName) = os.path.split(out_file)
    Path(out_filePath).mkdir(parents=True, exist_ok=True)
    
    force_delete_output_file = args.force_delete_output_file
    
    if (os.path.exists(out_file)):
        if force_delete_output_file is True :
            os.remove(out_file)
        else:
            print ("The out_file you specified already exists and -force_delete_output_file was not used. Use '-h' for help.")
            sys.exit() 
    Path(out_filePath).mkdir(parents=True, exist_ok=True)

    sequence_type = args.sequence_type
    if not sequence_type:
        sequence_type = "aa"
    feature_type = args.feature_type
    if not feature_type:
        feature_type = "CDS"
    file_format = args.file_format
    if not file_format:
        file_format = "genbank"
        
        
    codon_table_if_unspecified_in_input_file = args.codon_table_if_unspecified_in_input_file
    conffile = args.config
    if conffile:
        codon_table_from_conffile = get_codon_table(conffile)
        if codon_table_if_unspecified_in_input_file :
            if codon_table_if_unspecified_in_input_file != codon_table_from_conffile:
                raise RuntimeError ("Error: The codon table specified in -codon_table_if_unspecified_in_input_file is different from the one specified in -config.")
                #sys.exit()
        else :
            codon_table_if_unspecified_in_input_file = codon_table_from_conffile
    if not codon_table_if_unspecified_in_input_file:
        codon_table_if_unspecified_in_input_file = 11

    
    # for in_file in glob.iglob(in_dir + '**/**', recursive=True):
    #     in_file = os.path.abspath(in_file)
    
    #if in_file.is_file():
    if os.path.isfile(in_file):
        
        in_file_handle = open (in_file, 'r', encoding='utf-8')
        out_file_handle = open (out_file, 'w')
        
        record_iterator = SeqIO.parse(in_file_handle, file_format)
        
        # CDS counter for numbering features along the genome
        # Assuming features of the seqRecord are ordered by location on the genome
        counter = 1
        
        dictLocusTags = {}
        dictGenomeIds = {}
        genome_accession_rank = 0
        for record in record_iterator:
            genome_accession_rank += 1
            #print ("\tConverting '%s' to fasta ..." % record.description)
            if sequence_type in ['nt', 'aa', 'taa']:
                (fasta_records, counter, dictLocusTags, dictGenomeIds) = gb_to_faa(record, sequence_type, feature_type, codon_table_if_unspecified_in_input_file, counter, dictLocusTags, dictGenomeIds, in_file, genome_accession_rank)
            elif sequence_type == 'whole': #whole records are handled separetly
                raise RuntimeError ("whole sequence_type is not supported for now.")
                #fasta_records = gb_to_faa_whole(record, "organism,accessions")
            else:
                raise RuntimeError ("Unrecognized sequence_type. Use '-h' for help.")
                #sys.exit()
            SeqIO.write(fasta_records, out_file_handle, 'fasta')
        
        in_file_handle.close()
        out_file_handle.close()
