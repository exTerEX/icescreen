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
usage: bash edit_configfile.sh [--resgff_ori <Path to original result gff file> --main_outdir <Path to main outdir> --gb_name <Name of the genbank file> --resgff_reorganized <Path to reorganized result gff file> --icescreen_conf_reorganized <Path to reorganized icescreen conf file>][--help]
EOF
}

########################################
# VARIABLES INITIALIZATION
########################################

# Initialize all option variables.
# This ensures we are not contaminated by variables from the environment.
# Path to original result gff file
resgff_ori=
# Path to main outdir
main_outdir=
# Name of the genbank file
gb_name=
# Path to reorganized result gff file
resgff_reorganized=
# Path to reorganized icescreen conf file
icescreen_conf_reorganized=

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
        --resgff_ori)
            if [ "$2" ]; then
                resgff_ori=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --resgff_ori=?*)
            resgff_ori=${1#*=}
            ;;
        --resgff_ori=)
            die 'ERROR: "--resgff_ori" requires a non-empty option argument.'
            ;;
        --main_outdir)
            if [ "$2" ]; then
                main_outdir=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --main_outdir=?*)
            main_outdir=${1#*=}
            ;;
        --main_outdir=)
            die 'ERROR: "--main_outdir" requires a non-empty option argument.'
            ;;
        --gb_name)
            if [ "$2" ]; then
                gb_name=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --gb_name=?*)
            gb_name=${1#*=}
            ;;
        --gb_name=)
            die 'ERROR: "--gb_name" requires a non-empty option argument.'
            ;;
        --resgff_reorganized)
            if [ "$2" ]; then
                resgff_reorganized=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --resgff_reorganized=?*)
            resgff_reorganized=${1#*=}
            ;;
        --resgff_reorganized=)
            die 'ERROR: "--resgff_reorganized" requires a non-empty option argument.'
            ;;
        --icescreen_conf_reorganized)
            if [ "$2" ]; then
                icescreen_conf_reorganized=$2
                shift
            else
                die 'ERROR: '$1' requires a non-empty option argument.'
            fi
            ;;
        --icescreen_conf_reorganized=?*)
            icescreen_conf_reorganized=${1#*=}
            ;;
        --icescreen_conf_reorganized=)
            die 'ERROR: "--icescreen_conf_reorganized" requires a non-empty option argument.'
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

mkdir -p "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats"

if [ -f "${main_outdir}/faa/${gb_name}.faa" ]
then
	mv "${main_outdir}/faa/${gb_name}.faa" "${main_outdir}/results/${gb_name}/${gb_name}.faa"
else
	die 'ERROR: '"${main_outdir}/faa/${gb_name}.faa"' does not exist'
fi

if [ -f "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.log" ]
then
	touch "${main_outdir}/${gb_name}/param.conf"
	tail -n +8 "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.log" > "${main_outdir}/${gb_name}/param.conf"
	rm -f "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.log"
else
	die 'ERROR: '"${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.log"' does not exist'
fi

if [ -f "${main_outdir}/results/${gb_name}/icescreen.conf" ]
then
	# mv "${main_outdir}/results/${gb_name}/icescreen.conf" "${main_outdir}/${gb_name}/icescreen.conf"
	tail -n +25 "${main_outdir}/results/${gb_name}/icescreen.conf" >> "${main_outdir}/${gb_name}/param.conf"
	rm -f "${main_outdir}/results/${gb_name}/icescreen.conf"
	gzip "${main_outdir}/${gb_name}/param.conf"
else
	die 'ERROR: '"${main_outdir}/results/${gb_name}/icescreen.conf"' does not exist'
fi

if [ -f "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_SP_withMEIds.tsv" ]
then
	mv "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.summary" "${main_outdir}/${gb_name}/detected_mobile_elements/${gb_name}_detected_ME.summary"
	mv "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_ME.tsv" "${main_outdir}/${gb_name}/detected_mobile_elements/${gb_name}_detected_ME.tsv"
	mv "${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_SP_withMEIds.tsv" "${main_outdir}/${gb_name}/detected_mobile_elements/${gb_name}_detected_SP_withMEIds.tsv"
	rmdir "${main_outdir}/results/${gb_name}/icescreen_detection_ME"
else
	die 'ERROR: '"${main_outdir}/results/${gb_name}/icescreen_detection_ME/${gb_name}_detected_SP_withMEIds.tsv"' does not exist'
fi

if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff" ]
then

	gzip "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.embl"
	if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.embl" ]
	then
		rm -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.embl"
	fi
	mv "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.embl.gz" "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats/${gb_name}_icescreen.embl.gz"
	
	gzip "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gb"
	if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gb" ]
	then
		rm -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gb"
	fi
	mv "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gb.gz" "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats/${gb_name}_icescreen.gb.gz"
	
	gzip "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff"
	if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff" ]
	then
		rm -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff"
	fi
	mv "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff.gz" "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats/${gb_name}_icescreen.gff.gz"

	gzip "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.fa"
	if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.fa" ]
	then
		rm -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.fa"
	fi
	mv "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.fa.gz" "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats/${gb_name}_source.fa.gz"

	gzip "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.gff"
	if [ -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.gff" ]
	then
		rm -f "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.gff"
	fi
	mv "${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_source.gff.gz" "${main_outdir}/${gb_name}/detected_mobile_elements/standard_genome_annotation_formats/${gb_name}_source.gff.gz"
	
	rmdir "${main_outdir}/results/${gb_name}/visualization_files"
else
	die 'ERROR: '"${main_outdir}/results/${gb_name}/visualization_files/${gb_name}_icescreen.gff"' does not exist'
fi

if [ -d "${main_outdir}/results/${gb_name}/icescreen_detection_SP" ]
then
	# echo tar -czf "${main_outdir}/${gb_name}/tmp_intermediate_files.tar.gz" --quiet -C "${main_outdir}" "${main_outdir}/results/${gb_name}/"
	tar -C "${main_outdir}" -czf "${main_outdir}/${gb_name}/tmp_intermediate_files.tar.gz" "results/${gb_name}/"
	rm -rf "${main_outdir}/results/${gb_name}"
else
	die 'ERROR: '"${main_outdir}/results/${gb_name}/icescreen_detection_SP/${gb_name}_icescreen.gff"' does not exist'
fi


