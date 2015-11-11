red='\e[0;31m'
NC='\e[0m' # No Color
green='\e[0;32m'

echo "configuring with clang"
./configure CC="ccache clang" CXX="ccache clang++" --enable-mpi  2> /dev/null  > /dev/null
make clean > /dev/null 2> /dev/null
echo "makeing" 
make -j4 > /dev/null   2> /dev/null

if [ $? != 0 ]; then
    echo -e  "[ ${red} failure ${NC} ] with clang"
else  
    echo -e  "[ ${green} OK ${NC} ] with clang"
fi



