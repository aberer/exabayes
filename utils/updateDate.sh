#! /bin/bash

date=$( date "+%F" )
commitId=$(git log  | head -n 1  | cut -f 2  -d ' ' )

echo "#define RELEASE_DATE \"$date\"" > src/releaseDate.hpp
echo "#define GIT_COMMIT_ID \"$commitId\"" >> src/releaseDate.hpp

echo -e  "\n\nupdated release date: "
cat src/releaseDate.hpp
echo -e  "\n\n"


