#!/bin/bash

if [ ! -d "$1" ]
then
    echo "No target directory specified" 1>&2
    exit 1
fi

SCRIPT_PATH="$(dirname "$(readlink -f "$0")")"

find "$1" -mindepth 1 -maxdepth 1 -type d -print0 | while IFS= read -r -d '' dir; do
    "${SCRIPT_PATH}/generate-sif.sh" "${dir}"/*/*pdb.gz
    "${SCRIPT_PATH}/sif-to-database.py" "${dir}/database" "${dir}"/*/*.sif
done

"${SCRIPT_PATH}/sif-to-database.py" -c "$1/classification-database" "$1"/*/*/*.sif
