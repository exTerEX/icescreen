#!/bin/bash
set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

########################################
# FUNCTIONS
########################################
# die function (show message then quit).
function die() {
    printf '%s\n' "$1" >&2
    exit 1
}

# usage / help for parameters
function usage() {
cat << EOF
usage: bash edit_configfile.sh [--conffile <Path to icescreen.conf> --gbfile <Path to genome genbank file> --meconffile <Path to icescreen.conf of detection ME> --spfile <Path to detected SP tsv file> --mefile <Path to detected ME tsv file> --spidsfile <Path to detected ME SP with ME Ids tsv file> --melogfile <Path to detected ME log file> --srcfastafile <Path to genome source fasta file> --srcgfffile <Path to genome source gff3 file> --gffannotfile <Path to ICEscreen results GFF3 annotation file> --emblannotfile <Path to ICEscreen results EMBL annotation file> --gbannotfile <Path to ICEscreen results GB annotation file>][--help]
EOF
}

########################################
# VARIABLES INITIALIZATION
########################################

# Initialize all option variables.
# This ensures we are not contaminated by variables from the environment.
# Path to icescreen.conf
conffile=
# Genome genbank file
gbfile=
# Path to icescreen.conf of detection ME
meconffile=
# Path to detected SP tsv file
spfile=
# Path to detected ME tsv file
mefile=
# Path to detected ME SP with ME Ids tsv file
spidsfile=
# Path to detected ME log file
melogfile=
# Path to genome source fasta file
srcfastafile=
# Path to genome source gff3 file
srcgfffile=
# Path to ICEscreen results GFF3 annotation file
gffannotfile=
# Path to ICEscreen results EMBL annotation file
emblannotfile=
# Path to ICEscreen results GB annotation file
gbannotfile=

########################################
# ARGUMENTS PARSING
########################################

# Catch all arguments given and check if they are empty.
while :; do
    case $1 in
        -h|-\?|--help)
            usage       # Show help.
            exit 0
            ;;
        --conffile)
            if [ "$2" ]; then
                conffile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --conffile=?*)
            conffile=${1#*=}
            ;;
        --conffile=)
            die 'ERROR: "--conffile" requires a non-empty option argument.'
            ;;
        --gbfile)
            if [ "$2" ]; then
                gbfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --gbfile=?*)
            gbfile=${1#*=}
            ;;
        --gbfile=)
            die 'ERROR: "--gbfile" requires a non-empty option argument.'
            ;;
        --spfile)
            if [ "$2" ]; then
                spfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --spfile=?*)
            spfile=${1#*=}
            ;;
        --spfile=)
            die 'ERROR: "--spfile" requires a non-empty option argument.'
            ;;
        --mefile)
            if [ "$2" ]; then
                mefile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --mefile=?*)
            mefile=${1#*=}
            ;;
        --mefile=)
            die 'ERROR: "--mefile" requires a non-empty option argument.'
            ;;
        --spidsfile)
            if [ "$2" ]; then
                spidsfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --spidsfile=?*)
            spidsfile=${1#*=}
            ;;
        --mecsvfile=)
            die 'ERROR: "--spidsfile" requires a non-empty option argument.'
            ;;
        --melogfile)
            if [ "$2" ]; then
                melogfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --melogfile=?*)
            melogfile=${1#*=}
            ;;
        --melogfile=)
            die 'ERROR: "--melogfile" requires a non-empty option argument.'
            ;;
        --srcfastafile)
            if [ "$2" ]; then
                srcfastafile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --srcfastafile=?*)
            srcfastafile=${1#*=}
            ;;
        --srcfastafile=)
            die 'ERROR: "--cdscov" requires a non-empty option argument.'
            ;;
        --srcgfffile)
            if [ "$2" ]; then
                srcgfffile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --srcgfffile=?*)
            srcgfffile=${1#*=}
            ;;
        --srcgfffile=)
            die 'ERROR: "--srcgfffile" requires a non-empty option argument.'
            ;;
        --gffannotfile)
            if [ "$2" ]; then
                gffannotfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --gffannotfile=?*)
            gffannotfile=${1#*=}
            ;;
        --gffannotfile=)
            die 'ERROR: "--gffannotfile" requires a non-empty option argument.'
            ;;
        --emblannotfile)
            if [ "$2" ]; then
                emblannotfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --emblannotfile=?*)
            emblannotfile=${1#*=}
            ;;
        --emblannotfile=)
            die 'ERROR: "--emblannotfile" requires a non-empty option argument.'
            ;;
        --gbannotfile)
            if [ "$2" ]; then
                gbannotfile=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --gbannotfile=?*)
            gbannotfile=${1#*=}
            ;;
        --gbannotfile=)
            die 'ERROR: "--gbannotfile" requires a non-empty option argument.'
            ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            die
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac
    
    shift
done

if [ -z "$conffile" ]
then
    die 'ERROR: Config file (--conffile) is required'
else
    if [ -f "$conffile" ]
    then
        # Genbank orignal and sanitized names
        if [ ! -z "$gbfile" ]
        then
            sed -i 's| genome_genbank:| genome_genbank: '"$gbfile"'|' $conffile
        fi
        # Path to detected SP csv file
        if [ ! -z "$spfile" ]
        then
            sed -i 's|detected_sp:|detected_sp: '"$spfile"'|' $conffile
        fi
        # Path to detected ME tsv file
        if [ ! -z "$mefile" ]
        then
            sed -i 's|detected_me:|detected_me: '"$mefile"'|' $conffile
        fi
        # Path to detected ME SP with ME Ids csv file
        if [ ! -z "$spidsfile" ]
        then
            sed -i 's|detected_sp_with_me_ids:|detected_sp_with_me_ids: '"$spidsfile"'|' $conffile
        fi
        # Path to detected ME log file
        if [ ! -z "$melogfile" ]
        then
            sed -i 's|log_file:|log_file: '"$melogfile"'|' $conffile
        fi
        # Path to genome source fasta file
        if [ ! -z "$srcfastafile" ]
        then
            sed -i 's|source_sequence_fasta:|source_sequence_fasta: '"$srcfastafile"'|' $conffile
        fi
        # Path to genome source gff3 file
        if [ ! -z "$srcgfffile" ]
        then
            sed -i 's|source_annotation_gff3:|source_annotation_gff3: '"$srcgfffile"'|' $conffile
        fi
        # Path to ICEscreen results GFF3 annotation file
        if [ ! -z "$gffannotfile" ]
        then
            sed -i 's|icescreen_annotation_gff3:|icescreen_annotation_gff3: '"$gffannotfile"'|' $conffile
        fi
        # Path to ICEscreen results EMBL annotation file
        if [ ! -z "$emblannotfile" ]
        then
            sed -i 's|icescreen_annotation_embl:|icescreen_annotation_embl: '"$emblannotfile"'|' $conffile
        fi
        # Path to ICEscreen results GB annotation file
        if [ ! -z "$gbannotfile" ]
        then
            sed -i 's|icescreen_annotation_genbank:|icescreen_annotation_genbank: '"$gbannotfile"'|' $conffile
        fi
    else
        die 'ERROR: Config file '"$conffile"' does not exist'
    fi
fi


