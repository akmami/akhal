_akhal() {
    local bin cur cmds

    # whatever the binary is called on this line, so ./akhal and a renamed or second-version binary complete against their own command list
    bin="${COMP_WORDS[0]}"
    cur="${COMP_WORDS[COMP_CWORD]}"
    COMPREPLY=()

    # only the word right after `akhal` names a command
    if [ "$COMP_CWORD" -eq 1 ]; then
        cmds=$("$bin" -h 2>/dev/null | awk '/^  /{print $1}')
        COMPREPLY=( $(compgen -W "${cmds}" -- "$cur") )
    fi

    return 0
}

complete -o default -F _akhal akhal
