#!/bin/bash
_biosnake() {
	local cur
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"

	if [[ ${COMP_CWORD} -eq 1 ]] ; then
		COMPREPLY=( $(LC_COLLATE=C compgen -W "$(biosnake shellcompletion 2> /dev/null)" -- "${cur}") )
		return 0
	fi

	if [[ ${COMP_CWORD} -gt 1 ]] ; then
		COMPREPLY=( $(LC_COLLATE=C compgen -f -W "$(biosnake shellcompletion "${COMP_WORDS[1]}" 2> /dev/null)" -- "${cur}") )
		return 0
	fi

}
complete -o plusdirs -F _biosnake biosnake
