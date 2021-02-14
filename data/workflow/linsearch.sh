#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}


#pre processing
[ -z "$BIOSNAKE" ] && echo "Please set the environment variable \$BIOSNAKE to your BIOSNAKE binary." && exit 1;
# check number of input variables
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exist
[ ! -f "$1.dbtype" ] && echo "$1.dbtype not found!" && exit 1;
[ ! -f "$2.dbtype" ] && echo "$2.dbtype not found!" && exit 1;
[   -f "$3.dbtype" ] && echo "$3.dbtype exists already!" && exit 1;
[ ! -d "$4" ] && echo "tmp directory $4 not found!" && mkdir -p "$4";

QUERY="$1"
TARGET="$2"
TMP_PATH="$4"

# 1. Finding exact $k$-mer matches.
if notExists "${TMP_PATH}/pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$BIOSNAKE" kmersearch "$QUERY" "$TARGET" "${TMP_PATH}/pref" ${KMERSEARCH_PAR} \
        || fail "kmermatcher died"
fi

# 1. Ungapped alignment filtering
RESULTDB="${TMP_PATH}/pref"
if [ -n "$FILTER" ]; then
    if notExists "${TMP_PATH}/reverse_ungapaln.dbtype"; then
        # shellcheck disable=SC2086
        $RUNNER "$BIOSNAKE" rescorediagonal "$TARGET" "$QUERY" "${RESULTDB}" "${TMP_PATH}/reverse_ungapaln" ${RESCORE_FILTER_PAR} \
        || fail "Rescorediagonal step died"
    fi

    if notExists "${TMP_PATH}/pref_filter.dbtype"; then
        "$BIOSNAKE" filterdb "${TMP_PATH}/pref" "${TMP_PATH}/pref_filter" --filter-file "${TMP_PATH}/reverse_ungapaln" --positive-filter 0 \
            || fail "Filterdb step died"
    fi
    RESULTDB="${TMP_PATH}/pref_filter"
fi


# 2. Local gapped sequence alignment.
if notExists "${TMP_PATH}/reverse_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$BIOSNAKE" "${ALIGN_MODULE}" "$TARGET" "$QUERY" "${RESULTDB}" "${TMP_PATH}/reverse_aln" ${ALIGNMENT_PAR} \
        || fail "Alignment step died"
fi


# 3. swap logic
if [ -n "$NUCL" ]; then

    if notExists "${TMP_PATH}/aln.dbtype"; then
        # shellcheck disable=SC2086
        "$BIOSNAKE" swapresults "$TARGET" "$QUERY" "${TMP_PATH}/reverse_aln" "${TMP_PATH}/aln" ${SWAPRESULT_PAR} \
            || fail "Alignment step died"
    fi

    if [ -n "$FILTER" ]; then
        if notExists "${TMP_PATH}/reverse_ungapaln_swap.dbtype"; then
             "$BIOSNAKE" swapresults "$TARGET" "$QUERY" "${TMP_PATH}/reverse_ungapaln" "${TMP_PATH}/ungap_aln" \
              || fail "Mergedbs died"
        fi

        if notExists "${TMP_PATH}/aln_merged.dbtype"; then
             "$BIOSNAKE" concatdbs "${TMP_PATH}/ungap_aln" "${TMP_PATH}/aln" "${TMP_PATH}/aln_merged" --preserve-keys --take-larger-entry\
              || fail "Mergedbs died"
        fi
        RESULT="${TMP_PATH}/aln_merged"
    fi


    if notExists "$3.dbtype"; then
        # shellcheck disable=SC2086
        "$BIOSNAKE" offsetalignment "$QUERY" "${QUERY}" "${TARGET}" "${TARGET}" ${RESULT} "$3" ${OFFSETALIGNMENT_PAR} \
            || fail "Offset step died"
    fi
else
    # 3. Local gapped sequence alignment.
    if notExists "$3.dbtype"; then
        # shellcheck disable=SC2086
       "$BIOSNAKE" swapresults "$TARGET" "$QUERY" "${TMP_PATH}/reverse_aln" "$3" ${SWAPRESULT_PAR} \
            || fail "Alignment step died"
    fi
fi



if [ -n "$REMOVE_TMP" ]; then
    # shellcheck disable=SC2086
    "$BIOSNAKE" rmdb "${TMP_PATH}/pref" ${VERBOSITY}
    # shellcheck disable=SC2086
    "$BIOSNAKE" rmdb "${TMP_PATH}/reverse_aln" ${VERBOSITY}
    rm -f "${TMP_PATH}/linsearch.sh"
fi
