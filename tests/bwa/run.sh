#!/bin/bash
set -o nounset # Fail on unset variables
set -o errexit # Fail on uncaught non-zero returncodes
set -o pipefail # Fail is a command in a chain of pipes fails

echo cd $(dirname $0)

if [ ! -e "ValidateSamFile.jar" ];
then
	echo "Please place 'ValidateSamFile.jar' in the current folder."
	exit 1
fi

# Fetch and build relevant versions of BWA
function build_bwa ()
{
	local version=$1
	local bwaarch=bwa-${version}.tar.bz2
	local bwadir=bwa-${version}

	cd ./builds
	if [ ! -e "${bwadir}/bwa" ];
	then
		echo "Building BWA v${version}"

		if [ ! -e "${bwadir}" ];
		then
			if [ ! -e "${bwaarch}" ];
			then
				wget "http://sourceforge.net/projects/bio-bwa/files/bwa-${version}.tar.bz2/download" \
					-O ${bwaarch}
			fi

			tar xvjf bwa-${version}.tar.bz2
		fi

		nice -19 make -C bwa-${version}
	fi

	cd ..
}

mkdir -p builds

build_bwa 0.7.12
build_bwa 0.7.11
build_bwa 0.7.10
build_bwa 0.7.9a
build_bwa 0.7.8
build_bwa 0.7.7
build_bwa 0.7.6a # This version is broken (no aln command!)
build_bwa 0.7.5a
build_bwa 0.7.4
build_bwa 0.7.3a
build_bwa 0.7.2
build_bwa 0.7.1
# build_bwa 0.7.0 # Broken, does not compile
#build_bwa 0.6.2
#build_bwa 0.6.1
#build_bwa 0.5.9
#build_bwa 0.6.0
#build_bwa 0.5.10
#build_bwa 0.5.9rc1 # Oldest acceptable version


cd builds
if [ ! -e "bwa-git" ];
then
	git clone "https://github.com/lh3/bwa.git" "bwa-git"
fi

cd bwa-git
git pull
make
cd ../../


# Errors to ignore during valiation
IGNORE="IGNORE=RECORD_MISSING_READ_GROUP IGNORE=MISSING_READ_GROUP"

ls -d testcases/* |
while read testcase;
do
	echo "Running testcase ${testcase}: $(head -n1 ${testcase}/README)"
	ls builds/*/bwa |
	while read BWA;
	do
		echo -n "  $BWA "

		rm -rf temp
        folder="runs/$testcase/$(dirname $BWA | xargs basename)"
        rm -rf $folder
        mkdir -p $folder
        ln -s $folder temp

		msg=""
		returncode=-1

		cp ${testcase}/* temp/
		if [ -e "temp/run.sh" ];
		then
			bash "temp/run.sh" "$(pwd)/$BWA" && returncode=$? || returncode=$?

			if [ -e "temp/run.log" ];
			then
				msg="$(head -n1 temp/run.log)"
			fi
		elif [ -e "temp/reads1.fasta" ];
		then
			PREFIX=temp/prefix.fasta
			READS1=temp/reads1.fasta
			READS2=temp/reads2.fasta
			RESULTS=temp/results

            command="index"
			if $BWA index ${PREFIX} 2> ${PREFIX}.log;
			then
                command="aln #1"
				if $BWA aln ${PREFIX} ${READS1} > ${READS1}.fai 2> ${READS1}.log;
				then
                    command="aln #2"
					if $BWA aln ${PREFIX} ${READS2} > ${READS2}.fai 2> ${READS2}.log;
					then
                        command="sampe"
						if $BWA sampe ${PREFIX} ${READS1}.fai ${READS2}.fai ${READS1} ${READS2} 2> ${RESULTS}.log | \
							paleomix cleanup --paired --fasta ${PREFIX} --temp-prefix temp/cleanup 2> ${RESULTS}.cleanup.log >  ${RESULTS}.bam;
						then
							java -jar "ValidateSamFile.jar" ${IGNORE} I=${RESULTS}.bam &> ${RESULTS}.bam.validated \
								&& returncode=$? || returncode=$?

							if [ -e "${RESULTS}.bam.validated" ];
							then
								msg="$((grep ERROR ${RESULTS}.bam.validated || true) | head -n1)"
							fi
						fi
					fi
				fi
			fi
		fi

		if test $returncode -eq 0;
		then
			echo -e "\033[32m[OK]\033[0m"
		elif test $returncode -eq -1;
		then
			echo -e "\033[31m[TEST ERROR]: $command\033[0m"
		else
			echo -e "\033[33m[FAILED]\033[0m: $msg"
		fi
	done
done

rm -rf temp
