

## The dev directory
 - The dev directory is not to be pushed to the public repository of ICEscreen, it contains unused SP or detailed information about the classification of the SP or the delineation of the elements that are yet to be released.


---


## The database -> blastdb directory
This directory deals with the blastp detection of SP using curated SP protein sequences in fasta format.

#### to export "Firmicute_SP" table so biologist such as Gerard Guédon can modify it
 - Fichier -> exporter -> table vers un fichier csv
 - open this file with Excel and save as .xls format
 - send the file to Gérard

#### to import a modified version of the database
 - get the excel file modified by biologist such as Gerard Guédon
 - using excel, save it as a text file (tab delimited)
 - open ICE_Finder.db with software "DB-Browser for SQLite" (https://sqlitebrowser.org/)
 - right click on the table "Metadata_SP" and delete it
 - Fichier -> importer -> table depuis un fichier csv
 - name the table "Metadata_SP", tick box header as first line, then ok
 - check that the column name is "IME_superfamily"
 - save the modifications
 
#### to modify the ICE_IME_analysis/icescreen_dev/icescreen_detection_SP/dev/database/blastdb/ICE_Finder.db
 - to open and modify the file icescreen_detection_SP/database/blastdb/ICE_Finder.db, a software such as https://sqlitebrowser.org/ can be used. Run 'sqlitebrowser' in a terminal after installing the software to run it.

#### How to add a new SP to the blastp database by creating a new fasta file 

 - create a new fasta file with the protein sequences of the SP you want to add in the directory database -> blastdb. For example, see the file Couplingprot2.faa.

 - modify the file icescreen_pipelines/icescreen.snakefile. The changes must happen in the lines where "Couplingprot2" is present and follow the same coding pattern.
blast_dbnames = ["Relaxase", "Couplingprot", "Couplingprot2", ....

 - add a new file icescreen_pipelines/mode/[__TAXONOMY_DETECTION_MODE__].yml OR modify the appropriate existing [__TAXONOMY_DETECTION_MODE__].yml. Example of such file is icescreen_pipelines/mode/bacillota.yml. You can follow the links between the different sections of the yaml file by following the example of "Couplingprot2". 

 - if needed, modify **all** the files icescreen_pipelines/mode/[__TAXONOMY_DETECTION_MODE__].yml so that the new yaml sections added at the previous step are referenced across all the [__TAXONOMY_DETECTION_MODE__].yml files. In other words, the titles for the sections should be identical and present across all the taxonomic yaml files, only their value can be different. If a type of SP referenced by a section of a specific taxonomic yaml file (i.e. taxonomy_A.yml) should not be detected in another specific taxonomic yaml file (i.e. taxonomy_B.yml), then the threshold values for the latter should be such that no match can be considered valid.

 - for each protein entry, create an entry in icescreen_detection_SP/database/blastdb/ICE_Finder.db. The fields to create a new entry are as follow :
	"Prot_ID"	TEXT,
	"Prot_Type"	TEXT,
	"Associated_element_type"	TEXT,
	"ICE_superfamily"	TEXT,
	"ICE_family"	TEXT,
	"IME_superfamily"	TEXT,
	"relaxase_family_domain"	TEXT,
	"relaxase_family_MOB"	TEXT,
	"coupling_type"	TEXT,
	"False_positives"	TEXT
The superfamily field is used to forbid association of SP with two different superfamilies within the same ICE IME structure.
The family field is used to preferentially associate SP with similar family within the same ICE IME structure.

 - Notify your changes in the file /database/blastdb/LOG/summary_log_modif_SP_database.md and add a detailed log file as well in the directory /database/blastdb/LOG/detailled_log.

 - run icescreen --index_genomic_resources in order to index the genomic resources related to the signature proteins 
'''
icescreen --index_genomic_resources
'''

 - run a serie of tests to make sure your changes affect the sofware positively and does not constitue a regression.

#### How to add a new SP to the blastp database by modifying an existing fasta file

 - Modify an existing fasta file by adding the reference protein SP in the directory database -> blastdb. For example, see the file Couplingprot2.faa.

 - for each protein entry, create an entry in icescreen_detection_SP/database/blastdb/ICE_Finder.db. The fields to create a new entry are as follow :
	"Prot_ID"	TEXT,
	"Prot_Type"	TEXT,
	"Associated_element_type"	TEXT,
	"ICE_superfamily"	TEXT,
	"ICE_family"	TEXT,
	"IME_superfamily"	TEXT,
	"relaxase_family_domain"	TEXT,
	"relaxase_family_MOB"	TEXT,
	"coupling_type"	TEXT,
	"False_positives"	TEXT
The superfamily field is used to forbid association of SP with two different superfamilies within the same ICE IME structure.
The family field is used to preferentially associate SP with similar family within the same ICE IME structure.

 - Notify your changes in the file /database/blastdb/LOG/summary_log_modif_SP_database.md and add a detailed log file as well in the directory /database/blastdb/LOG/detailled_log.

 - run icescreen --index_genomic_resources in order to index the genomic resources related to the signature proteins 
'''
icescreen --index_genomic_resources
'''

 - run a serie of tests to make sure your changes affect the sofware positively and does not constitue a regression.


---


## The database -> hmmdb directory

This directory deals with the HmmScan detection of more distantly related SP using set of HMM profile, either imported or custom made.

#### The FP_profiles directory

It includes HMM profile for the detection of the false positives protein that are somewhat similar to SP but that are not SP.

#### The SP_profiles directory

It includes HMM profile for the detection of the signature proteins (SP).

#### How to add a new SP HMM profile resource

 - add a new HMM profile file (.hmm) to the directory database -> hmmdb -> SP_profiles. See Couplingprot_T4SS_t4cp1.hmm for example. The name of the file should follow the convention of the other file in this directory. Note that the convention of the naming is separated into three sections that are handled separetly in the code, i.e. "Couplingprot", "T4SS", and "t4cp1". Pay attention to the information fields in the HMM file as well as the filename for the HMM file itself. Make sure that in the HMM profile file (.hmm), the field "NAME" (line 2) contains the same string as the name of the HMM file itself. For example the file Couplingprot_T4SS_t4cp2.hmm has a "NAME" field of "Couplingprot_t4cp2".
 
 - add a new file icescreen_pipelines/mode/[__TAXONOMY_DETECTION_MODE__].yml OR modify the appropriate existing [__TAXONOMY_DETECTION_MODE__].yml. Example of such file is icescreen_pipelines/mode/bacillota.yml. 

 - if needed, modify **all** the files icescreen_pipelines/mode/[__TAXONOMY_DETECTION_MODE__].yml so that the new yaml sections added at the previous step are referenced across all the [__TAXONOMY_DETECTION_MODE__].yml files. In other words, the titles for the sections should be identical and present across all the taxonomic yaml files, only their value can be different. If a type of SP referenced by a section of a specific taxonomic yaml file (i.e. taxonomy_A.yml) should not be detected in another specific taxonomic yaml file (i.e. taxonomy_B.yml), then the threshold values for the latter should be such that no match can be considered valid.

 - Notify your changes in the file /database/blastdb/LOG/summary_log_modif_SP_database.md and add a detailed log file as well in the directory /database/blastdb/LOG/detailled_log.

 - run icescreen --index_genomic_resources in order to index the genomic resources related to the signature proteins 
'''
icescreen --index_genomic_resources
'''

 - run a serie of tests to make sure your changes affect the sofware positively and does not constitue a regression.



#### How to modify the min_query_length for a specific HMM profile

 - The first part of the HMM profile file name (before the first underscore) must reference a distinct name, i.e. RelaxaseShorter_SXT_MOBI.hmm if we have to decrease the min_query_length for the HMM profile relaxase for SXT.
 
 - run ./icescreen --index_genomic_resources to index this new HMM file
 
 - The min_query_length is referenced in the yaml files (i.e. bacillota.yml) and can be modified directly in it. In the hmmscan section of the yaml file, add "RelaxaseShorter: Relaxase" in params_SP_type as an extension of the repeated nodes SP_types
'''
    hmmscan:
        params_SP_type:
            <<: *SP_types
            RelaxaseShorter: Relaxase #extension repeated nodes
'''
In params_profiles, add a reference to the second part of the HMM profile file name (after the first underscore), it can be default as min_query_length is not referenced in this section.
'''
        params_profiles:
        # Relaxase
            [...]
            SXT_MOBI: default
'''
In params_usage, add a reference to RelaxaseShorter like in the example below :
'''
        params_usage:
            [...]
            RelaxaseShorter: relaxase_shorter
'''
In params_filtering:, add a relaxase_shorter section that reference the min_query_length, i.e:
'''
        params_filtering:
            [...]
            relaxase_shorter:
                min_query_length: 130
'''


#### How to build a new profile HMM (from Lao J, Lacroix T, Guédon G, Coluzzi C, Payot S, Leblond-Bourget N, Chiapello H. ICEscreen: a tool to detect Firmicute ICEs and IMEs, isolated or enclosed in composite structures. NAR Genom Bioinform. 2022 Oct 21;4(4):lqac079. doi: 10.1093/nargab/lqac079. PMID: 36285285; PMCID: PMC9585547.)

To create the new HMM profiles of relaxases, one or more reference sequences were selected from the ICEscreen protein database of SPs and/or from the NCBI GenBank public database. The functional domain architecture was identified using CD-Search and SPARCLE. Coding DNA Sequences (CDSs) from GenBank with a similar functional domain architecture were retrieved. Redundant sequences were removed by a clustering at 95% identity over 100% of the length using CD-HIT. To further decrease the number of seed sequences if needed, a second clustering at 40% identity over 100% of the length was performed. A multiple alignment was constructed with MAFFT v7.407 using the FFTNS1 algorithm. Manual curation of the multiple alignment was carried out to remove too divergent proteins by using a phylogenetic tree constructed with SeaView 4.2. The Neighbor-Joining method BioNJ was used with the following parameters: gaps excluded, Poisson distance and bootstrap of 100 iterations. When a set of relaxases was still too distant, the multiple alignment was further divided into smaller but coherent sets of alignments to subsequently create HMM profiles. When the functional domain was found at the N-terminal position, Clustal Omega v1.2.4 with default parameters was used instead of MAFFT to align subsets of sequences and minimize gaps in the N-terminal part of the alignment. To prune away the ends of the multiple alignments carrying little information on conservation, BMGE was used with a BLOSUM 30 substitution matrix. The HMM profile was then built with the HMMbuild tool from the HMMER 3.2.1 suite with default parameters.



---



