#!/bin/sh -e
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

if notExists "${TMP_PATH}/query.dbtype"; then
    # shellcheck disable=SC2086
    "$BIOSNAKE" createdb "${QUERY}" "${TMP_PATH}/query" ${CREATEDB_QUERY_PAR} \
        || fail "query createdb died"
    QUERY="${TMP_PATH}/query"
fi

if notExists "${TARGET}.dbtype"; then
    if notExists "${TMP_PATH}/target"; then
        # shellcheck disable=SC2086
        "$BIOSNAKE" createdb "${TARGET}" "${TMP_PATH}/target" ${CREATEDB_PAR} \
            || fail "target createdb died"
    fi
    TARGET="${TMP_PATH}/target"
fi

if notExists "${INTERMEDIATE}.dbtype"; then
    # shellcheck disable=SC2086
    "$BIOSNAKE" rbh "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${TMP_PATH}/rbh_tmp" ${SEARCH_PAR} \
        || fail "Search died"
fi

if notExists "${TMP_PATH}/alis.dbtype"; then
    # shellcheck disable=SC2086
    "$BIOSNAKE" convertalis "${QUERY}" "${TARGET}" "${TMP_PATH}/result" "${RESULTS}" ${CONVERT_PAR} \
        || fail "Convert Alignments died"
fi

if [ -n "${REMOVE_TMP}" ]; then
    # shellcheck disable=SC2086
    "$BIOSNAKE" rmdb "${TMP_PATH}/result" ${VERBOSITY}
    if [ -z "${LEAVE_INPUT}" ]; then
        if [ -f "${TMP_PATH}/target" ]; then
            # shellcheck disable=SC2086
            "$BIOSNAKE" rmdb "${TMP_PATH}/target" ${VERBOSITY}
            # shellcheck disable=SC2086
            "$BIOSNAKE" rmdb "${TMP_PATH}/target_h" ${VERBOSITY}
        fi
        # shellcheck disable=SC2086
        "$BIOSNAKE" rmdb "${TMP_PATH}/query" ${VERBOSITY}
        # shellcheck disable=SC2086
        "$BIOSNAKE" rmdb "${TMP_PATH}/query_h" ${VERBOSITY}
    fi
    rm -rf "${TMP_PATH}/rbh_tmp"
    rm -f "${TMP_PATH}/easyrbh.sh"
fi
