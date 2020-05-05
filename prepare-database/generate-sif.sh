#!/bin/bash
TMP_DIR_BASE="$(mktemp -d -t "prepare-data.XXXXX"|| exit 1)"
SCRIPT_PATH="$(dirname "$(realpath -s "$0")")"
RINERATOR_DIR="${SCRIPT_PATH}/../rinerator/Source/"

MESSAGES=""

SUFFIX=""

while getopts "s:" opt; do
  case $opt in
    s)
      SUFFIX="_${OPTARG}"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

shift $((OPTIND-1))

for RAW_INPUT in "$@"
do
    INPUT="$(realpath --no-symlinks "${RAW_INPUT}")"
    INPUT_DIR="$(dirname ${INPUT})"
    OUTPUT_DIR="$(dirname ${INPUT})"
    PDB_ID="$(basename ${INPUT} | sed -e 's/\([^_\.]*\).*/\1/g')"
    SEGMENTS_FILE="${OUTPUT_DIR}/${PDB_ID}.seg"

    if [ ! -f "${SEGMENTS_FILE}" ]
    then
        MESSAGES="${MESSAGES}${SEGMENTS_FILE} not found\n"
        MESSAGES="${MESSAGES}-> Skipping ${PDB_ID}\n"
        continue
    fi

    # RINerator
    TMP_OUTPUT_DIR="${TMP_DIR_BASE}/${PDB_ID}"
    mkdir -p "${TMP_OUTPUT_DIR}"


    PDB_FILE="${TMP_OUTPUT_DIR}/${PDB_ID}.pdb"
    if [[ "${INPUT}" == *.gz ]]
    then
        gunzip -c "${INPUT}" > "${PDB_FILE}"
    else
        # reduce modifies the original file, so use a copy
        cp "${INPUT}" "${PDB_FILE}"
    fi

    pushd "${RINERATOR_DIR}" > /dev/null
    ./get_segments.py "${PDB_FILE}" "${TMP_OUTPUT_DIR}" "${SEGMENTS_FILE}"
    popd > /dev/null

    mv "${TMP_OUTPUT_DIR}/${PDB_ID}_h.sif" "${OUTPUT_DIR}/${PDB_ID}${SUFFIX}.sif"
    rm -f "${PDB_FILE}"
done

rm -rf "${TMP_DIR_BASE}"

echo -e "$MESSAGES"
