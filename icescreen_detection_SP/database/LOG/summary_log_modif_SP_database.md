
## v1.3.0, January 2024
#### Minor changes
* [Issue #12](https://forgemia.inra.fr/ices_imes_analysis/icescreen/-/issues/12): the integrase ADV69694 was wrongfully annotated as having IME_superfamily PF01076, relaxase_family_domain PF01076, and relaxase_family_MOB MOBV. Those metadata were causing an error : "Input file error: adding Relaxase_family_domain_of_most_similar_ref_SP for Type Serine integrase". Those metadata have been deleted.

## v1.1.0, November 2022
* The IME_superfamily "PF01719-like" has been renamed "PF01719" to conform to the naming convention.
* Replacement of family PF01076+PF02486 by PF02486 and of family MobV+MobT by mobT for AAM99131, ABA45702, AFS45036, AGU82204, AGU83894, BAN61292, CAD45862, CBZ47880, CCW37084.
* Replacement of family PF01076+PF02486 by PF01076 and of family MobV+MobT by mobV for ABA44955, AFS45038, AGU82202, AGU83896, BAH80804, CAD45864, CBZ47878, CCW37086.
* Removed the improper assignment of relaxase family PF01441 to two SPs that are not relaxases but belonged to an IME encoding a relaxase PF01446: coupling protein AGP68114 and integrase AGP68118.
* Fixed error domain identifier PF01441 with correct identifier PF01446 for relaxase AGP68116.
* Removal of serine integrases from MGI SCCmec (BAA86648, BAC67561, BAG06212) whose transfer mechanisms are still unknown.
* Assignment of CP BAN75392, Relaxase BAN75394, VirB4 BAN7540, Integrase BAN75415 to an IME of Tn916 family instead of the ICE_Lcas393_rpsI family. The CP, relaxase and VirB4 having respectively 59% id, 46% id and 66% id with that of Tn916.
* Replacement of IME superfamily PF05840/PHA00330 encoding integrase CCW3779, Relaxase CCW37800, CP CCW37801 as well as superfamily PF05840/PHA00330 of relaxase CCW37800 by PHA00330. According to the bank and a CDSearch, CCW37800 has a PF05840 domain and not PHA00330. However, the analyzes carried out during the work for the article on the 124 streptococcal genomes demonstrated that the relaxases having this domain are related on the region corresponding to this domain to the region of the relaxases having the PHA00330 domain. Therefore considered as a particular family of PHA00330 relaxases in the publication, nomenclature to be included in ICEScreen.
* The relaxase family "PHA00330-like" is changed to "PHA00330" to conform to the naming convention.
* Replacement of AGZ23064 (wrong start, missing 26 aa) by WP_041179472.1.
* The integrase AAO81984 has been removed.
