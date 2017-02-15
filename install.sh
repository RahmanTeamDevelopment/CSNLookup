#!/bin/bash

message () {
echo; printf '%.0s-' {1..100}; echo
echo "$1"
printf '%.0s-' {1..100}; echo; echo
}

ver="v1.0.0"

pver=$(python -c "import sys; print sys.version[:3]")
if [ $pver != '2.7' ]; then
	echo; echo "CSN Lookup requires Python 2.7.x"; echo
	exit
fi

message "Installing CSN Lookup "$ver

mkdir submissions
mkdir transdbs

tar -zxvf CAVA-1.2.0.tar.gz
mv CAVA-1.2.0 cava
cd cava
./install.sh

message "CSN Lookup "$ver" installed."
