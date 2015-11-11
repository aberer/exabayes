#! /bin/bash

d=$( date "+%F" )

echo "#define RELEASE_DATE \"$d\"" > src/releaseDate.hpp

echo -e  "\n\nupdated release date: "
cat src/releaseDate.hpp
echo -e  "\n\n"

