import os
import glob

################################################################################
#                          INITIALIZATION OF VARIABLES                         #
################################################################################

#------------------------------ COMMON VARIABLES ------------------------------#
# Parameters of 'config':
# - outdir : Path to results folder 'ICEscreen_results/'
# - gbdir  : Path to data folder 'genbank' of 'outdir'
# - rootdir: Path to the 'icescreen' folder
# - mode   : Path to the YAML config file for a search mode

outdir = config["outdir"]
gbdir = config["gbdir"]
rootdir = config["rootdir"]
mode = config["mode"]

# Get genbank names (Must be without spaces)
gbnames = []
""" multilines comment
for gb in glob.glob(os.path.join(gbdir, "*")):
    gbname = os.path.basename(gb)
    gbname = os.path.splitext(gbname)[0]
    gbnames.append(gbname)

"""
for gb in glob.iglob(os.path.join(gbdir, "*"), recursive=False):
	if os.path.isfile(gb):
		gbname = os.path.basename(gb)
		gbname = os.path.splitext(gbname)[0]
		gbnames.append(gbname)
#		print("The file "+ str(gb) +" is registred to be analysed.")
#	else :
#		print("Warning : Subdirectory "+ str(gb) +" is ignored, only genbank files at the root of the input genbank directory are analysed.")


#--------------------------------- BLAST MODE ---------------------------------#

# Get Blast databases names
blast_dbnames = ["Relaxase", "Couplingprot", "Couplingprot2", "VirB4", "IntTyr", "IntSer", "DDE"]

#---------------------------------- MODE HMM ----------------------------------#

# Path to HMM profile file (.hmm)
db_sp_hmmprofiles = [hmmfile for hmmfile in glob.glob("{}/icescreen_detection_SP/database/hmmdb/SP_profiles/*hmm".format(rootdir))]

# Dictionary with key as path to hmm file and value as hmm profile name
PATH_SP_HMM_PROFILES = {}
for gbname in gbnames:
    for sp_hmmprofile in db_sp_hmmprofiles:
        PATH_SP_HMM_PROFILES[os.path.splitext(os.path.basename(sp_hmmprofile))[0]] = sp_hmmprofile

# Path to HMM profile file (.hmm)
db_fp_hmmprofiles = [hmmfile for hmmfile in glob.glob("{}/icescreen_detection_SP/database/hmmdb/FP_profiles/*hmm".format(rootdir))]

# Dictionary with key as path to hmm file and value as hmm profile name
PATH_FP_HMM_PROFILES = {}
for gbname in gbnames:
    for fp_hmmprofile in db_fp_hmmprofiles:
        PATH_FP_HMM_PROFILES[os.path.splitext(os.path.basename(fp_hmmprofile))[0]] = fp_hmmprofile

################################################################################
#                                     RULES                                    #
################################################################################

#--------------------------------- FINAL RULES --------------------------------#

# Rule all
rule all:
    input:
        expand(os.path.join(outdir, "results/{gbname}/visualization_files/{gbname}_icescreen.gff"), gbname=gbnames),
        expand(os.path.join(outdir, "results/{gbname}/icescreen.conf"), gbname=gbnames)


# Create annotation files for the visualization of ICEscreen results
rule create_annotation_files:
    input:
        gbfile = os.path.join(gbdir, "{gbname}.gb"),
        spfile = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/{gbname}_detected_SP.tsv"),
        mefile = os.path.join(outdir, "results/{gbname}/icescreen_detection_ME/{gbname}_detected_ME.tsv")
    params:
        outname = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_icescreen"),
        conffile = os.path.join(outdir, "results/{gbname}/icescreen.conf")
    output:
        srcfasta = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_source.fa"),
        srcgff = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_source.gff"),
        resgff = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_icescreen.gff"),
        resembl = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_icescreen.embl"),
        resgb = os.path.join(outdir, "results", "{gbname}", "visualization_files", "{gbname}_icescreen.gb")
    shell:
        """
        # Extract genome sequence from genbank into fasta file
        python3 {rootdir}/icescreen_formatting/genbank_to_gff3_and_fasta.py -i {input.gbfile} \
                                                                    --gff {output.srcgff} \
                                                                    --fasta {output.srcfasta}
        # Generate annotation files with ICEscreen results (GFF3, EMBL and Genbank files)
        python3 {rootdir}/icescreen_formatting/generate_annotation_files.py -s {input.spfile} -m {input.mefile} -g {input.gbfile} -o {params.outname}
        # Add output paths to ICEscreen config file
        bash {rootdir}/icescreen_pipelines/edit_configfile.sh --conffile {params.conffile} \
                                                    --srcfastafile {output.srcfasta} \
                                                    --srcgfffile {output.srcgff} \
                                                    --gffannotfile {output.resgff} \
                                                    --emblannotfile {output.resembl} \
                                                    --gbannotfile {output.resgb}
        """


#------------------ DETECTION ME (MOBILE ELEMENTS) PART RULES -----------------#

# Run icescreen_OO.py to search mobile elements
rule detect_mobile_elements:
    input:
        spfile = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/{gbname}_detected_SP.tsv"),
        meconffile = os.path.join(rootdir, "icescreen_detection_ME/icescreen.conf"),
        gb = os.path.join(gbdir, "{gbname}.gb")
    params:
        conffile = os.path.join(outdir, "results/{gbname}/icescreen.conf")
    log:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_ME/{gbname}_detected_ME.log")
    output:
        mefile = os.path.join(outdir, "results/{gbname}/icescreen_detection_ME/{gbname}_detected_ME.tsv"),
        spidsfile = os.path.join(outdir, "results/{gbname}/icescreen_detection_ME/{gbname}_detected_SP_withMEIds.tsv")
    shell:
        """
        # Search ICE and IME
        python3 {rootdir}/icescreen_detection_ME/src/icescreen_OO.py -i {input.spfile} -c {input.meconffile} -o {output.mefile} -m {output.spidsfile} -l {log} --gb_input {input.gb} --taxo_mode_file {mode}
        # Add output paths to ICEscreen config file
        bash {rootdir}/icescreen_pipelines/edit_configfile.sh --conffile {params.conffile} \
                                                    --mefile {output.mefile} \
                                                    --spidsfile {output.spidsfile} \
                                                    --melogfile {log}
        """


# ---------------- DETECTION SP (SIGNATURE PROTEINS) PART RULES ----------------#

########----------------- PREPROCESSING SP PART RULES ------------------########



# Extract CDS protein sequence from a genbank (genome)
# Inputs:
#   - GenBank file (.gb)
# Output: Multifasta proteic file (.faa)
rule gb_to_faa:
    input:
        gb = os.path.join(gbdir, "{gbname}.gb")
    output:
        os.path.join(outdir, "faa", "{gbname}.faa")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/gb_to_faa.py -i {input.gb} -o {output} -c {mode}
        """


########------------ BLAST MODE OF DETECTION SP PART RULES -------------########

# Run blastP to search SPs
# Input: Multifasta proteic file (.faa)
# Output: TSV file (.tsv)
rule search_blastP_hits:
    input:
        os.path.join(outdir, "faa", "{gbname}.faa")
    params:
        blast_db = os.path.join(rootdir, "icescreen_detection_SP/database/blastdb", "{blast_dbname}")
    output:
        os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP", "Blast_mode", "blastp_output", "{gbname}_{blast_dbname}.tsv")
    shell:
        """
        # Run BlastP
        blastp -out {output} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
               -db {params.blast_db} \
               -evalue 0.001 \
               -query {input} \
               -task blastp \
               -max_target_seqs 10;

        # Add header to Blastp output
        blastpfields="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen"

        # If there is results (file exist), add header
        if [ -s {output} ]
        then
            sed -i "1s/^/$blastpfields\\n/" {output}
        # If there is no results (file does not exist), create file with header only
        else
            echo $blastpfields > {output}
        fi
        """

# Process BlastP results
rule process_blastp_results:
    input:
        os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode/blastp_output", "{gbname}_{blast_dbname}.tsv")
    output:
        outall = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode/unfiltered_results", "{gbname}_{blast_dbname}_unfiltered.tsv"),
        outfiltered = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode/filtered_results", "{gbname}_{blast_dbname}_filtered.tsv"),
        outbest = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode/filtered_results", "{gbname}_{blast_dbname}_filtered_best.tsv")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/process_blastp_results.py -i {input}\
                                                                       -c {mode}\
                                                                       -d {rootdir}/icescreen_detection_SP/database/blastdb/ICE_Finder.db\
                                                                       --outall {output.outall}\
                                                                       --outfiltered {output.outfiltered}\
                                                                       --outbest {output.outbest}
        """

# Gather all results for one genome into one tsv file
rule gather_blastP_SP_results:
    input:
        expand(os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode/filtered_results", "{gbname}_{blast_dbname}_filtered_best.tsv"), blast_dbname=blast_dbnames, allow_missing=True)
    output:
        temp(os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode", "{gbname}_SP.tmp.tsv")),
        temp(os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode", "{gbname}_SP.tmp2.tsv")),
        os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/Blast_mode", "{gbname}_blast_SP.tsv")
    shell:
        """
        # Concatenate all results
        awk 'FNR==1 && NR!=1 {{next;}}{{print}}' {input} > {output[0]};
        # Reorder rows
        tail -n+2 {output[0]} | sort -t$'\t' -k2n | cat <(head -1 {output[0]}) - > {output[1]}
        python3 {rootdir}/icescreen_detection_SP/src/retain_only_best_blast_hit_for_each_locus_tag.py -i {output[1]} -o {output[2]}
        """


########------------- HMM MODE OF DETECTION SP PART RULES --------------########

# Run hmmscan with CDS as query
rule search_hmmscan_hits:
    input:
        os.path.join(outdir, "faa", "{gbname}.faa")
    params:
        hmmprofile = lambda w: PATH_SP_HMM_PROFILES[w.hmmname]
    output:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/HMM_mode/hmmscan_output/{gbname}_{hmmname}.txt"),
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/HMM_mode/hmmscan_output/{gbname}_{hmmname}.tsv")
    shell:
        """
        # Run hmmscan
        hmmscan --domtblout {output[0]} {params.hmmprofile} {input} > /dev/null
        # Format to tsv
        sed '/^#/d' {output[0]} |\
        awk -F '[[:space:]]+' 'BEGIN{{ OFS="\t"; }} {{ descr=$23; for(i=24; i<=NF; i++) descr=descr" "$i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,descr; }}' |\
        awk 'BEGIN{{ print "target_name\taccession\ttlen\tquery_name\tseq_accession\tqlen\tseq_E-value\tseq_score\tseq_bias\t#_domain\tof_domain\tdomain_c-Evalue\tdomain_i-Evalue\tdomain_score\tdomain_bias\thmm_coord_from\thmm_coord_to\tali_coord_from\tali_coord_to\tenv_coord_from\tenv_coord_to\tenv_coord_acc\tdescription_of_target" }}1' > {output[1]};
        """


# Gather all results for one genome into one tsv file
rule gather_hmmscan_results:
    input:
        expand(os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/HMM_mode/hmmscan_output", "{gbname}_{hmmname}.tsv"), hmmname=PATH_SP_HMM_PROFILES.keys(), allow_missing=True)
    output:
        os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/HMM_mode/hmmscan_output", "{gbname}.tsv")
    shell:
        """
        # Concatenate all results
        awk 'FNR==1 && NR!=1 {{next;}}{{print}}' {input} > {output};
        """


# Filter out hmmscan hits
rule process_hmmscan_results:
    input:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/HMM_mode/hmmscan_output/{gbname}.tsv")
    output:
        outall = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/HMM_mode/unfiltered_results", "{gbname}_unfiltered.tsv"),
        outfiltered = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/HMM_mode/filtered_results", "{gbname}_filtered.tsv"),
        outbest = os.path.join(outdir, "results", "{gbname}", "icescreen_detection_SP/HMM_mode", "{gbname}_hmm_SP.tsv")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/process_hmmscan_results.py -i {input}\
                                                                        -c {mode}\
                                                                        --outall {output.outall}\
                                                                        --outfiltered {output.outfiltered}\
                                                                        --outbest {output.outbest}
        """


########------------ HITS CLEANING DETECTION SP PART RULES -------------########

# Merge SP results obtained with the Blast mode and HMM mode into folder for hits cleaning
rule merge_SP_results:
    input:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/Blast_mode/{gbname}_blast_SP.tsv"),
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/HMM_mode/{gbname}_hmm_SP.tsv")
    output:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_source.tsv")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/merge_SP_results.py --blastres {input[0]} --hmmres {input[1]} -o {output}
        """


# Merge SP results obtained with the Blast mode and HMM mode into folder for hits cleaning
rule extract_SP_faa:
    input:
        tsvfile = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_source.tsv"),
        faafile = os.path.join(outdir, "faa", "{gbname}.faa")
    output:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_source.faa")
    shell:
        """
        nrow=`cat {input.tsvfile} | wc -l`
        if [ $nrow -eq 1 ]
        then
            touch {output}
        else
            protid=$(tail -n +2 {input.tsvfile} | cut -f3 -d'\t' | sed 's/^/\>/' | sed 's/$/_\\\\/' | paste -s -d'|' | sed 's/^/\//' | sed 's/\\\\$/\/,+1p/')
            sed -n "$protid" {input.faafile} > {output}
        fi
        """


# Run hmmscan with CDS as query
rule search_hmmscan_FP:
    input:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_source.faa")
    params:
        hmmprofile = lambda w: PATH_FP_HMM_PROFILES[w.hmmname]
    output:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/{gbname}_{hmmname}.txt"),
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/{gbname}_{hmmname}.tsv")
    shell:
        """
        nrow=`cat {input} | wc -l`
        if [ $nrow -eq 0 ]
        then
            touch {output[0]};
            touch {output[1]};
        else
            hmmscan --domtblout {output[0]} {params.hmmprofile} {input} > /dev/null
            sed '/^#/d' {output[0]} |\
            awk -F '[[:space:]]+' 'BEGIN{{ OFS="\t"; }} {{ descr=$23; for(i=24; i<=NF; i++) descr=descr" "$i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,descr; }}' |\
            awk 'BEGIN{{ print "target_name\taccession\ttlen\tquery_name\tseq_accession\tqlen\tseq_E-value\tseq_score\tseq_bias\t#_domain\tof_domain\tdomain_c-Evalue\tdomain_i-Evalue\tdomain_score\tdomain_bias\thmm_coord_from\thmm_coord_to\tali_coord_from\tali_coord_to\tenv_coord_from\tenv_coord_to\tenv_coord_acc\tdescription_of_target" }}1' > {output[1]};
        fi
        """


# Gather all false positives results for one genome into one tsv file
rule gather_hmmscan_fp:
    input:
        expand(os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_output/{gbname}_{hmmname}.tsv"), hmmname=PATH_FP_HMM_PROFILES.keys(), allow_missing=True)
    output:
        os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_compiled_output/{gbname}.tsv")
    shell:
        """
        # Concatenate all results
        awk 'FNR==1 && NR!=1 {{next;}}{{print}}' {input} > {output};
        """


# Filter out hmmscan hits
rule process_hmmscan_fp:
    input:
        insp = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_source.tsv"),
        infp = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/hmmscan_compiled_output/{gbname}.tsv")
    output:
        outfp = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/proteins_to_remove/{gbname}_hmm_fp_hits.tsv"),
        outfiltered = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_hmm_cleaned.tsv")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/process_hmmscan_fp.py --insp {input.insp}\
                                                                   --infp {input.infp}\
                                                                   -c {mode}\
                                                                   --outfp {output.outfp}\
                                                                   --outfiltered {output.outfiltered}
        """


# Reannot SP
rule reannot_SP:
    input:
        insp = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_hmm_cleaned.tsv"),
        intyr = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/Blast_mode/filtered_results/{gbname}_IntTyr_filtered.tsv")
    params:
        conffile = os.path.join(outdir, "results/{gbname}/icescreen.conf")
    output:
        outall = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/hits_cleaning/{gbname}_detected_SP_hmm_cleaned_reannotated.tsv"),
        outfiltered = os.path.join(outdir, "results/{gbname}/icescreen_detection_SP/{gbname}_detected_SP.tsv")
    shell:
        """
        python3 {rootdir}/icescreen_detection_SP/src/reannot_SP.py -a {input.insp}\
                                                           -b {input.intyr}\
                                                           -c {mode}\
                                                           --outall {output.outall}\
                                                           --outfiltered {output.outfiltered}
        bash {rootdir}/icescreen_pipelines/edit_configfile.sh --conffile {params.conffile} \
                                                    --spfile {output.outfiltered}
        """
