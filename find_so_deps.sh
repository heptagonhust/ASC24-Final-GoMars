#!/bin/bash

# Usage:
#   ./find_so_deps.sh path/to/execuable
# Example:
#   ./find_so_deps.sh ./gmcore/build/gmcore_driver.exe

find_dependenies() {
    name="$1"
    path="$2"
    echo "$name:"
    objdump -x "$path" 2> /dev/null | grep NEEDED | grep -oP '(  \S*\.so\S*)'
}

[ -n "$1" ] || exit 1
[ -f "$1" ] || exit 1

IFS=$'\n'

for line in $(ldd "$1" | grep -oP '.*\.so.* => \/.* \(.*\)'); do
    name="$(grep -oP '(lib.*so.*)(?= =>)' <<< "$line")"
    path="$(grep -oP '(?<= => )(.*so.*)(?= \()' <<< "$line")"
    find_dependenies "$name" "$path"
done

