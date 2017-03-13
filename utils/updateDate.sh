#! /bin/bash

root=$(dirname $(dirname $0))

echo $root


isGit=$?

if [ -d $root/.git ]; then
    commitId=$(git log  | head -n 1  | cut -f 2  -d ' ' )
    date=$(git log --date=iso -1  | grep Date  | sed 's/Date:[ ]*\(.*\)/\1/' | cut -f 1,2 -d ' ')
else
    version=$(cat $root/VERSION)
    commitId="-- not built from git (version: $version) --"
    date="$(date +'%Y-%m-%d %H:%M:%S')"
fi

echo "#define RELEASE_DATE \"$date\"" > src/releaseDate.hpp
echo "#define GIT_COMMIT_ID \"$commitId\"" >> src/releaseDate.hpp

echo -e  "\n\nupdated release date: "
cat src/releaseDate.hpp
echo -e  "\n\n"
