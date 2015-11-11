#! /bin/bash

d=$( date "+%F" )

echo "#define RELEASE_DATE \"$d\"" > src/releaseDate.hpp
