#!/bin/bash

#ICEscreen copyright Université de Lorraine - INRAE
#This file is part of ICEscreen. ICEscreen is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License Affero as published by the Free Software Foundation version 3 of the License. ICEscreen is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License Affero for more details. You should have received a copy of the GNU General Public License Affero along with ICEscreen.  If not, see <https://www.gnu.org/licenses/>.

function die
{
    local test_number=$1
    local my_command=${@:2}
    [ -z "$my_command" ] && message="Died"
    echo "ERROR ${BASH_SOURCE[1]}: line ${BASH_LINENO[0]}" >&2
    echo "test $test_number failed for command:" >&2
    echo "$my_command" >&2
    echo "with message:" >&2
    cat tmpTestIceScreenOO__.log >&2
    rm -f tmpTestIceScreenOO__.log
    exit 1
}

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "Starting tests:" >&2
rm -f tmpTestIceScreenOO__.log
touch tmpTestIceScreenOO__.log

: '
 - Test 1 :
       - What is tested : the file tests/test_files/test1_NC_004668.1_detected_SP.tsv contains 6 ICE canonical conjugaison module (including 2 that are guest / host) and 2 partial conjugaison module
       - Output files :
'

test_number="1"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test1_NC_004668.1_detected_SP.tsv -o $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.tsv -m $SCRIPT_DIR/results/test1_NC_004668.1_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.log"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test1_NC_004668.1_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test1_NC_004668.1_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test1_NC_004668.1_detected_SP_withICEIMEIds.tsv"
# Check columns: "ICE_IME_id","ICE_IME_id_need_manual_curation","Segment_number","Comments_ICE_IME_structure","Is_hit_blast","Is_hit_HMM","CDS_num","CDS","Id_of_blast_most_similar_ref_SP"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

#Remark: adding diff -I \"Date\" does not work for some reason...
#Remark: do not work either... my_command="diff <(sed '/Date/d' $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary) <(sed '/Date/d' $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary)"
my_command="diff $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary > $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary > $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test1_NC_004668.1_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.tsv $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.log $SCRIPT_DIR/results/test1_NC_004668.1_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test1_NC_004668.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 2 : 
	- What is tested : This test includes (among other) a host ICE without int that has 3 guests elements (3 IME) ; it also test for integrase IME that have different rules than int ICE (no strand check but distance check) ; it also test for relaxase that are not adjacent but separated by 1 CDS
	- Output files :
'
test_number="2"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test2_NZ_AP013072_detected_SP_ANNOTATED.tsv -o $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.tsv -m $SCRIPT_DIR/results/test2_NZ_AP013072_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test2_NZ_AP013072_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary $SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary > $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.summary > $SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test2_NZ_AP013072_results_EM.summary_WITHOUT_DATE"
#$my_command >> tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.tsv $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.log $SCRIPT_DIR/results/test2_NZ_AP013072_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test2_NZ_AP013072_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 3 : 
	- What is tested : This test includes (among other) an ICE with int DDE and an IME with composed Element_family (PF01719-PF00910) and single Element_family (PF01719)
	- HowOutput files :
'
test_number="3"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test3_Streptococcus_salivarius_strain_LAB813_detected_SP.tsv -o $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.tsv -m $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary $SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary > $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary > $SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.tsv $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.log $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test3_Streptococcus_salivarius_strain_LAB813_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 4 : 
	- What is tested : host / guest ICES with only the integrase being separated from the SP of the conjugaison module for the host ICE
	- How to run the test (the diff command should yield an empty line) :
'
test_number="4"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test4_NZ_AFRY01000001_detected_SP_ANNOTATED.tsv -o $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.tsv -m $SCRIPT_DIR/results/test4_NZ_AFRY01000001_detected_SP_ANNOTATED_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test4_NZ_AFRY01000001_detected_SP_ANNOTATED_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_detected_SP_ANNOTATED_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary $SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary > $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.summary > $SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test4_NZ_AFRY01000001_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.tsv $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.log $SCRIPT_DIR/results/test4_NZ_AFRY01000001_detected_SP_ANNOTATED_withICEIMEIds.tsv $SCRIPT_DIR/results/test4_NZ_AFRY01000001_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 5 : 
	- What is tested : complex zone with partial conjugaison modules
	- How to run the test (the diff command should yield an empty line) :
'
test_number="5"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test5_NC_022239_detected_SP.tsv -o $SCRIPT_DIR/results/test5_NC_022239_results_EM.tsv -m $SCRIPT_DIR/results/test5_NC_022239_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test5_NC_022239_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test5_NC_022239_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test5_NC_022239_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test5_NC_022239_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary $SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary > $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.summary > $SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test5_NC_022239_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test5_NC_022239_results_EM.tsv $SCRIPT_DIR/results/test5_NC_022239_results_EM.log $SCRIPT_DIR/results/test5_NC_022239_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test5_NC_022239_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 6 : 
	- What is tested : test case number 2 without false positive WP_080655072.1 for merged ICE #3
	- How to run the test (the diff command should yield an empty line) :
'
test_number="6"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test6_medley_detected_SP_ANNOTATED.tsv -o $SCRIPT_DIR/results/test6_medley_results_EM.tsv -m $SCRIPT_DIR/results/test6_medley_detected_SP_ANNOTATED_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test6_medley_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test6_medley_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test6_medley_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test6_medley_detected_SP_ANNOTATED_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test6_medley_detected_SP_ANNOTATED_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test6_medley_results_EM.summary $SCRIPT_DIR/tests/test_results/test6_medley_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test6_medley_results_EM.summary > $SCRIPT_DIR/results/test6_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test6_medley_results_EM.summary > $SCRIPT_DIR/tests/test_results/test6_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test6_medley_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test6_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test6_medley_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test6_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test6_medley_results_EM.tsv $SCRIPT_DIR/results/test6_medley_detected_SP_ANNOTATED_withICEIMEIds.tsv $SCRIPT_DIR/results/test6_medley_results_EM.log $SCRIPT_DIR/results/test6_medley_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2



: '
 - Test 7 : 
	- What is tested : complex case with coupling in conflict between 2 ICE IME
	- How to run the test (the diff command should yield an empty line) :
'
test_number="7"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test7_NC_015977.1_detected_SP.tsv -o $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.tsv -m $SCRIPT_DIR/results/test7_NC_015977.1_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test7_NC_015977.1_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test7_NC_015977.1_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test7_NC_015977.1_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary $SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary > $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.summary > $SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test7_NC_015977.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.tsv $SCRIPT_DIR/results/test7_NC_015977.1_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.log $SCRIPT_DIR/results/test7_NC_015977.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2



: '
 - Test 8 : 
	- What is tested : merge case can not choose at first between different structure to merge ; integrase fetched from same family but not directly upstream or downstream need to check strand if virb4
	- How to run the test (the diff command should yield an empty line) :
'
test_number="8"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test8_medley_detected_SP.tsv -o $SCRIPT_DIR/results/test8_medley_results_EM.tsv -m $SCRIPT_DIR/results/test8_medley_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test8_medley_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test8_medley_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test8_medley_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test8_medley_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test8_medley_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test8_medley_results_EM.summary $SCRIPT_DIR/tests/test_results/test8_medley_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test8_medley_results_EM.summary > $SCRIPT_DIR/results/test8_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test8_medley_results_EM.summary > $SCRIPT_DIR/tests/test_results/test8_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test8_medley_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test8_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test8_medley_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test8_medley_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test8_medley_results_EM.tsv $SCRIPT_DIR/results/test8_medley_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test8_medley_results_EM.log $SCRIPT_DIR/results/test8_medley_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2


: '
 - Test 9 : 
	- What is tested : real difficult case, curated by Gerard
	- How to run the test (the diff command should yield an empty line) :
'
test_number="9"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test9_NZ_LT635479.1_detected_SP.tsv -o $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.tsv -m $SCRIPT_DIR/results/test9_NZ_LT635479.1_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test9_NZ_LT635479.1_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary $SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary > $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.summary > $SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test9_NZ_LT635479.1_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.tsv $SCRIPT_DIR/results/test9_NZ_LT635479.1_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.log $SCRIPT_DIR/results/test9_NZ_LT635479.1_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2






: '
 - Test 10 : 
	- What is tested : test10_fragmented_SPs_detected_SP.tsv
'
test_number="10"
echo " ** Starting test #$test_number" >&2
my_command="python3 $SCRIPT_DIR/src/icescreen_OO.py -c $SCRIPT_DIR/icescreen.conf -i $SCRIPT_DIR/tests/test_files/test10_fragmented_SPs_detected_SP.tsv -o $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.tsv -m $SCRIPT_DIR/results/test10_fragmented_SPs_detected_SP_withICEIMEIds.tsv -l $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.log"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.tsv"
reference="$SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.tsv"
my_command="diff $result $reference"
#my_command="diff <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $result) <(cut -f1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20 $reference)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
result="$SCRIPT_DIR/results/test10_fragmented_SPs_detected_SP_withICEIMEIds.tsv"
reference="$SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_detected_SP_withICEIMEIds.tsv"
my_command="diff $result $reference"
#my_command="diff <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $reference) <(awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$14}' FPAT=\"([^,]+)|(\\\"[^\\\"]+\\\")\" $result)"
eval $my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi

my_command="diff $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary $SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary > $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="sed '/Date/d' $SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.summary > $SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="diff $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
#my_command="rm -f $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE $SCRIPT_DIR/tests/test_results/test10_fragmented_SPs_results_EM.summary_WITHOUT_DATE"
#$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
my_command="rm -f $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.tsv $SCRIPT_DIR/results/test10_fragmented_SPs_detected_SP_withICEIMEIds.tsv $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.log $SCRIPT_DIR/results/test10_fragmented_SPs_results_EM.summary"
$my_command > tmpTestIceScreenOO__.log 2>&1
if [ -s tmpTestIceScreenOO__.log ]
then
      die $test_number $my_command
fi
echo " ** Done with test #$test_number" >&2





# end of script
rm -f tmpTestIceScreenOO__.log
echo "Tests successful" >&2
