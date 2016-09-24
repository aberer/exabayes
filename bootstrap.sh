#! /bin/bash

gitroot=$(dirname $0)

if [ ! -d "$gitroot/.git" ]; then

    echo ".git folder not found. Most likely you are attempting to"
    echo "bootstrap from a github generated tar-ball which is not"
    echo "supported. Please download the source-tarball from the"
    echo "contributor's site or try to bootstrap in a cloned repository."

    exit
fi

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
