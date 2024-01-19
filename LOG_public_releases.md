
To access the different releases, see https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/releases


## v1.3.0, January 2024
#### Minor changes
* The following conda dependencies were upgraded: pandas >=2.1.4, and snakemake-minimal >=8.2.1.
* The command line argument "--mode" has been renamed "--phylum" to better reflect its scope and is now case insensitive. The command line argument "--mode" is still supported but raises a deprecation warning. Absence of the argument "--phylum" (or "--mode") in the command line now raises a deprecation warning, this argument is to become mandatory in a future release. The option "Phylum of the genomes to analyze" is now visible and expanded by default on the galaxy user interface. "Firmicutes" has been replaced by "Bacillota" to keep up with the current taxonomy nomenclature. "Firmicutes" is still supported for "--phylum" but raises a deprecation warning.
* The columns "Start_of_most_upstream_SP" and "Stop_of_most_downstream_SP" were added in the _detected_ME.tsv output file. Those columns are **not to be mistaken for the start and stop of the element however**.
* Genome accession that generates no significant alignment match with any of the ICEscreen reference signature proteins now produces a complete summary output file filled with 0 values instead of an empty output file.
* New ICEscreen command line option: "--print_version_dependencies" prints the version of the dependencies.
* The ICEscreen command line option "--test_installation" now works for user without OS root privileges.
* The experimental support for "Streptomyces" has been removed for now as it needs more polishing.
* v1.3.0 (January 2024) of the reference signature proteins database, see log file icescreen_detection_SP/database/LOG/summary_log_modif_SP_database.md for more details.
#### Bugs and warnings fixes
* [Issue #15](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/15) and [Issue #16](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/16): incompatibility between bcbio-gff version 0.7.0 and biopython version 1.82 and 1.83. The version of those two conda dependencies were set as follow: bcbio-gff =0.7.0, biopython =1.81.
* If a query locus tag had multiple valid alignment matches from multiple different fasta files of ICEscreen reference signature proteins, all those matches were considered and generated duplicated lines for a query locus tag in the final output files. This has been fixed and only the overall best alignment match for a query locus tag is now considered.
* Ascending sorting by i-Evalue of multiple HMM profile matches for a query locus tag was not working properly. As a result, the best scoring HMM profile was not always selected. HMM profile matches are now sorted by domain_score (descending) to ensure the selection of the best HMM profile.
* [Issue #13](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/13): in some case, there are multiple ICE or IME element structures that could be validly attributed to a query signature protein. When that happens, this query signature protein is attributed multiple potential ICE_IME_id which should appear in the column "ICE_IME_id_need_manual_curation" of the file _detected_SP_withICEIMEIds.tsv. A bug was fixed where the multiple potential ICE_IME_id would appear in the column "ICE_IME_id" when multiple such query signature protein were adjacent to each other.
* [Issue #14](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/14): generate_annotation_files.py 'DataFrame' object has no attribute 'tolist'. This bug happens for a Genbank file with multiple accessions (multiple chromosomes or contigs) and at least one of the accession has no significant alignment match with any of the ICEscreen reference signature proteins.
* A genome accession that generates no significant alignment match with any of the ICEscreen reference signature proteins produced empty visualization files (gb, embl, gff). It now produces valid gb, embl, and gff files with no ICEscreen result.
* Fixed warning "FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.". Updated dependencies for pandas >=2.1.0.
* Fixed warning "relative file path is strongly discouraged. It can also lead to inconsistent results of the file-matching approach used by Snakemake."
#### Other source code improvements
* The module to detect the structure of the mobile element has been updated: the source code is now generic for the type of signature proteins from the conjugation module (relaxase, coupling protein, virB4). This modification will make the integration of other type of signature proteins from the conjugation module easier in the future.


## v1.2.0, April 2023
#### Major changes
* Merging of the integrase fragments that are detectable by blast (i.e. insertion of mobile element in the integrase).
#### Minor changes
* The version of the conda dependencies were upgraded: bcbio-gff >=0.7.0, biopython >=1.81. Python 3.11 is not yet supported by some other dependencies.
* Do not add family or superfamily information to the structure if the signature protein is a fragment or pseudogene.
* Backtracking adding SP in conflict toward upstream is now done recursively. Previsouly, backtracking adding SP toward upstream was only considering N-1 structure for SP in conflict.	
#### Bug fixes
* Fixed Error in rule detect_mobile_elements: RuntimeError: determineIfResultSPFileHasMultipleGenomeAccesion error: len(dictGenomeAccesion) == 0. Genbank files that yield no hit with the ICEscreen database return this error instead of an empty result file. [Link in gitlab](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/7).
* Fixed Error in IsThereAnIntegraseBetweenThoseTwoConjModule: unrecognized positioning ofICEsIMEsStructureOne. This error arise while checking for the presence of an integrase between two adjacent structure that share the same SP (SP in conflict). [Link in gitlab](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/8).
* Fixed Error in mergeWith: not greenLightAddSPConjugaisonModule and not currSP.setICEsIMEsStructureInConflict for merging locus tag. This error happen when merging a distant fragmented SP. [Link in gitlab](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/9).
* Fixed Error in icescreen_detection_ME/src/rulesMergeICEIMEStructures.py : AttributeError: type object 'BasicEMStructure' has no attribute 'getGetListInternIdFromListEMStructure'. This error happens when a fragmented SP could be attributed to multiple adjacent SP from the conj module. [Link in gitlab](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/10).
* Fixed RuntimeError: Error in IsThereAnIntegraseBetweenThoseTwoConjModule: unrecognized positioning ofICEsIMEsStructureOne. This error occur when ICEsIMEsStructureTwo is upstream of ICEsIMEsStructureOne but they share one or more SP (SP in conflict) [Link in gitlab](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/11).
* Fixed BiopythonWarning: Partial codon, len(sequence) not a multiple of three.
* Fixed BiopythonDeprecationWarning: Use int(feature.end) and int(feature.start) rather than feature.nofuzzy_end and feature.nofuzzy_start.

## v1.1.1, February 2023
#### Minor changes
* Set version of dependency biopython =1.80 instead of biopython >=1.80. The dependency biopython version 1.81 is incompatible with bcbio-gff version 0.6.9, it throws an ImportError: cannot import name 'UnknownSeq' from 'Bio.Seq'. See https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/6 for more details.


## v1.1.0, November 2022
#### Major changes
* v1.1.0 (November 2022) of the signature proteins database, see log file icescreen_detection_SP/database/blastdb/LOG/summary_log_modif_SP_database.md for more details.
* Support for multigenbank files (i.e. gbff files featuring multiple genome or contig records back to back).
* Improvements regarding the assignation of the integrases to the conjugation modules in some specific cases: (1) take into account the possibility for a multiple succession of integrases in between conjugation modules (only attribute the integrases if there are no ambiguity) ; (2) take into account the possibility for two sequential integrases to form a host/guest relationship with 2 downstream or upstream sequential conjugation modules.
* Detection of some fragmented SP of the conjugation module and merge into a structure.
* If available, the locus tag is now used as the primary identifier for the gene (was formerly protein_id + start).
* Sequential numbering of ICE and IME elements IDs found by ICEscreen (no gap in numbering).
#### Minor changes
* Changes of the column names in the output files _detected_SP_withMEIds.tsv and _detected_ME.tsv to be less ambiguous and more precise.
* Standardization of the format of empty values to "-".
* Clarification in the output files whether the detected SPs are annotated as pseudogenes or not. Do not take pseudogenes into account when assigning the element type.
* In the column "Category_of_element", display the first letter of the types of the SP for the element specifically rather than all the combinations of the types of the SP related to this category.
* In the results_EM.tsv output file, the "Category_of_element" now mentions "(no integrase assigned)" when relevant to take into account for the following possibilities: (1) the absence of integrase around the detected element or (2) the impossibility for the program to choose between multiple valid integrases.
* Display of the family of the ICE only if there is a consensus between the SPs of the conjugation module that compose it. Only consider data on ICE families for CDS with more than 60% identity with the most similar reference SP.
* The "ICE Family" is associated to SPs in the report only if their blast alignment identity is >=40%.
* "IME superfamily" is now mentioned instead of "IME family" in the reports to better reflect the scope for this field.
* The "ICE superfamily" and "ICE family" of the integrase are not mentioned in the report anymore to avoid the confusion about their modular nature and the fact that they can associate to any ICE or IME.
* The formulation of comments were changed regarding isolated SP not associated to any structure.
* The threshold for the CDS distance between 2 consecutive SP for IME is increased from 10 to 14 (based on manual curation).
* Bug fix: incorrect formatting of output file _icescreen.gff: missing line breaks (https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/4).
* Bug fix: the stops of the SPs detected by ICEscreen were off by 1 nt in the genbank, embl, and gff reports.
* Bug fix: nested directory structure for input genbank files was causing problems. Now only genbank files that are at the root of the input genbank directory are taken into account (sub-directories are ignored).
* Bug fix: very low E-values were wrongly rounded down to 0.0 in the output files, they are now correctly displayed in a scientific format.
* Warning fix: File path ...ICEscreen_results//genbank/XXX.gb contains double '/'. This is likely unintended. It can also lead to inconsistent results of the file-matching approach used by Snakemake.
* Warning fix: Not prepending group keys to the result index of transform-like apply. In the future, the group keys will be included in the index, regardless of whether the applied function returns a like-indexed object.
* Five genomes were added to test for the coherence of the pipeline.
* The documentation was updated to reflect the changes for this version.
* The version of the conda dependencies were upgraded: python >=3.9, pandas >=1.5.2, snakemake-minimal >=7.18.2, biopython >=1.80, bcbio-gff >=0.6.9.
 

## v1.0.4, March 2022
This minor update fixes issues related to the InvalidIndexError error of v1.0.3 due to pandas version 1.4.0 and up. ICEscreen can now use the latest pandas version. 


## v1.0.3, January 2022
ICEscreen v1.0.3 stabilizes Pandas version 1.3.5 to avoid a bug in version 1.4.0 and blast and HMMer version to prevent future bugs related to changes to their underlying database schema.
* Stabilize Pandas version to 1.3.5 because of the appearance of a bug with version 1.4.0: when running ICEscreen on some genomes (i.e. NZ_CP026548.1.gb) that contains a tyrosine_integrase that resemble XerS, a pandas.errors.InvalidIndexError: Int64Index([0], dtype='int64' is thrown upon `data.at[idx, "False_positives"] = "-"` which was not thrown with prior versions.
* Stabilize blast version to 2.12 and HMMer version to 3.3.2 to avoid potential issues with changes in the underlying database schema between the indexing step at compile time and the execution at run-time.


## v1.0.2, January 2022
This minor update fixes issues related to the ordering of Genbank features by genomic position and BLAST Database error pre-fetching sequence data :
* Genbank features in the Genbank genome file need to be ordered by genomic position prior to runnning ICEscreen otherwise the results are inaccurate.
* BLAST Database error pre-fetching sequence data: Blast requirement changed from version >=2.9.0 to >=2.12. The issue comes down to a change in the blast DB format. NCBI updated their default DB format from V4 to V5. blast 2.9 was the last version that read (and wrote) V4 by default (optionally it would read/write V5). From version 2.11 V4 support is dropped all together. As a consequence, some instance of ICEscreen tends to exhibit the BLAST Database error pre-fetching sequence data.


## v1.0.1, December 2021
This minor update fixes compatibility issues with Galaxy and running ICEscreen without any write privilege on the installation directory.


## v1.0.0, December 2021
The version 1.0.0 is the first public release of the tool. ICEscreen detects and annotates ICEs (Integrative and Conjugative Elements) and IMEs (Integrative and Mobilizable Elements) in Firmicutes genomes. Its main features are:
* Detection of signature proteins (SPs) of ICEs/IMEs by using blastP on a curated resource. BlastP allows for an accurate assignment of hits to a given ICE/IME superfamily or family. The curated resource was derived from an analysis of over 120 ICEs and IMEs in Streptococcus genomes by the DINAMIC lab.
* Detection of distant homologs of SPs by using HMM profiles of ICEs/IMEs protein families. The HMM profiles have been either imported from trusted resources or created and curated when needed.
* Detection of the ICE/IME structures: ICEScreen groups together SPs that belong to the same ICE/IME structures to the best of its ability.
* Delimitation of the elements at the gene or nucleotide level is not yet implemented and still needs manual curation.
