#!/bin/bash

set -u


SDIR="$( dirname "$( cd "$( dirname "$0" )" && pwd )" )"


while read PROJECT NAME GENOME
do
  if [[ $PROJECT =~ ^\# ]]; then
    continue
  fi
  bash "$SDIR"/lib/mapping.s1seq.nextseq.20240711.sh \
    $PROJECT $NAME $GENOME
  tput bel
done <<EOS
# PROJECT NAME GENOME
Project_XXX SAMPLE01 mm10
Project_XXX SAMPLE02 mm10
EOS
