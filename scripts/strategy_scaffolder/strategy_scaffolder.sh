
if [[ $# -eq 0 ]] ; then
    echo 'Please provide a strategy name as the first argument.'
    exit 0
fi

NAME=$1
LNAME=`echo "$NAME" | tr '[:upper:]' '[:lower:]'`

echo Making new strategy with name ${NAME} in result/
RScript -e "source('./new_strategy.R'); new_strategy('${NAME}')"

echo Strategy $NAME created
echo New files created:
echo "\tR/${LNAME}.R"
echo "\tsrc/${LNAME}_strategy.cpp"
echo "\tinst/include/plant/${LNAME}_strategy.h"
echo "\ttests/testthat/test-strategy-${LNAME}.R"
echo
echo Modified Files:
echo "\tinst/include/plant.h"
echo "\tinst/include/RccpR6_classes.yml"
echo "\tsrc/plant_plus.cpp"
echo "\tsrc/plant_tools.cpp"
echo 
echo run the following in the root directory to recompile and test plant.
echo "make clean; make; make test; make install"