
countoutput=`ls -1 *.sh.o* 2>/dev/null | wc -l`
counterror=`ls -1 *.sh.e* 2>/dev/null | wc -l`
countcsv=`ls -1 *.csv 2>/dev/null | wc -l`

if [ $# -gt 0 ] && [ "$1" = "all" ]
then 
    if [ $(($countoutput + $counterror + $countcsv)) -gt 0 ]
    then
        echo "Cleaning logs and CSV files"
        rm *.sh.o* *.sh.e* *.csv
    else
        echo "No files to clean"
    fi
else
    if [ $(($countoutput + $counterror)) -gt 0 ]
    then
        echo "Cleaning log files"
        rm *.sh.o* *.sh.e*
    else
        echo "No files to clean"
    fi

fi 

