#!/bin/bash

# VERSION=0.1.2

# Files look something like this after Galaxy sets permissions via external_chown_script.
# # file: test
# # owner: letaw
# # group: HPCUsers
# user::rw-
# user:galaxyuser:rw-
# group::---
# group:galaxy:r--
# mask::rw-
# other::---

echo "Invoking ACL change operation." | tee $4

if [ "$#" -lt 4 ]; then
    echo "Needs 4 arguments <user/group> <user/group_name> <path> <logfile>" | tee $4
    exit -1;
fi

modifier="$1"
if [[ $modifier != 'g' && $modifier != 'u' ]]; then
    printf "Invalid option: %s\nMust be either g or u.\n" "$modifier" | tee $4
    exit -1;
fi
name="$2"
path="$3"

echo "Current ACLs:" | tee $4
getfacl -p $path | tee $4

if [ -e "$path" ]; then
    if [ -d "$path" ]; then
	path="${path}/";
    fi
    #Recursively provide read access
    setfacl -RP -m ${modifier}:${name}:r "${path}"
    if [ -d "$path" ]; then
	#Provide execute access to top level directory
	setfacl -P -m ${modifier}:${name}:rx "${path}"
	#Provide execute access to all sub-directories
	find ${path} -type d -exec setfacl -P -m ${modifier}:${name}:rx {} \;
    fi
else
    printf "File or directory %s not found" "${path}" | tee $4
    exit -1;
fi

echo "Modified ACLs:" | tee $4
getfacl -p $path | tee $4

echo "Access change complete." | tee $4
