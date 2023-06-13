#!/bin/bash
TMPDIR="/tmp"
JUNCTIONFILTEREDFILE=$4
seq=$1
str=$2
end=$3
end=$((end+1))
tmpfile=$(basename `mktemp`)
awk -v s=$seq '$1==s' $JUNCTIONFILTEREDFILE | awk -v s=$str '$2==s' | awk -v e=$end '$5==e' | cut -f3,6 > ${TMPDIR}/tmpfile
awk -v s=$seq '$1==s' junction.filter1 | awk -v s=$str '$5==s' | awk -v e=$end '$2==e' | cut -f3,6 >> ${TMPDIR}/tmpfile
strand=$(sort ${TMPDIR}/tmpfile | uniq -c | sort -k1,1nr | head -n1 | cut -f2)
if [ "$strand" = "+" ]; then
	echo "$strand"
elif [ "$strand" = "-" ];then
	echo "$strand"
else
	echo "."
fi
rm -f ${TMPDIR}/tmpfile
