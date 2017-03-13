#! /bin/bash

gitroot=$(dirname $0)

# use this script to compile the source code downloaded from the
# repository.

cd $gitroot

autoreconf -vfi && ./configure && make date

if [ ! -f $gitroot/src/releaseDate.hpp ]; then
    read -r -d  '' MSG <<EOF
\n
Bootstrapping failed. Please check the messages above for
errors (exabayes most likely will not compile).

Did you install autoconf and autoconf-archive (must be released later
than 2015)?

EOF

    echo -e "$MSG"
    exit

else

    read -r -d '' MSG <<EOF
\n
Bootstrap successful. You can now run configure and make.
EOF

    echo -e "$MSG"
    exit

fi
