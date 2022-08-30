
countoutput=`ls -1 *.sh.o* 2>/dev/null | wc -l`
counterror=`ls -1 *.sh.e* 2>/dev/null | wc -l`
countcsv=`ls -1 *.csv 2>/dev/null | wc -l`
countslurm=`ls -1 slurm-*.out 2>/dev/null | wc -l`

if [ $# -gt 0 ] && [ "$1" = "all" ]
then 
    if [ $(($countoutput + $counterror + $countcsv + $countslurm)) -gt 0 ]
    then
        echo "Cleaning logs and CSV files"
        
  		if [ $countoutput -gt 0 ]
  		then
  			rm *.sh.o*
  		fi
  		
  		if [ $counterror -gt 0 ]
  		then
  			rm *.sh.e*
  		fi
  		
  		if [ $countcsv -gt 0 ]
  		then
  			rm *.csv
  		fi
  		
  		if [ $countslurm -gt 0 ]
  		then
  			rm slurm-*.out
  		fi
    else
        echo "No files to clean"
    fi
else
    if [ $(($countoutput + $counterror + $countslurm)) -gt 0 ]
    then
        echo "Cleaning log files"
  		if [ $countoutput -gt 0 ]
  		then
  			rm *.sh.o*
  		fi
  		
  		if [ $counterror -gt 0 ]
  		then
  			rm *.sh.e*
  		fi
  		
  		if [ $countslurm -gt 0 ]
  		then
  			rm slurm-*.out
  		fi
    else
        echo "No files to clean"
    fi

fi 

