#!/bin/sh

#  sync.sh
#  
#
#  Created by Andrzej Chlebicki on 23.12.2018.
#  


# Sources
# rsync -r -avz --prune-empty-dirs --exclude '*.html' --exclude 'obj/' --exclude '*.csv' --exclude '*.png' --exclude '.*' --exclude 'Release' --exclude 'Debug' . ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalIntegrator


# Make the binary

options=(-avz --prune-empty-dirs )
destination="ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalIntegrator/"
results="./results"
if [ $# -gt 0 ] && [ "$1" == "mysz" ]
then 
    options+=('-e' 'ssh -A -J ac357729@kruk-host.fuw.edu.pl')
    destination="achleb@10.42.11.24:~/FunctionalIntegrator/"
    results="./results_mysz" 
fi
    
if [ $# -gt 0 ] && [ "$1" == "cluster" ]
then
    destination="achlebicki@cluster.uy:/clusteruy/home/achlebicki/FunctionalIntegrator/"
    results="./results_cluster" 
fi  

while ! [[ $confirm =~ [nNyY] ]]
    do 
	    read -p "Remake project? (y/n): " confirm
    done
    if [[ $confirm =~ [yY] ]]
    then
	    make clean
	    make static
    fi

echo "Synchronizing with ${destination}"

rsync "${options[@]}" --exclude '*.nb' --exclude '*.html' --exclude 'obj' \
    --exclude 'LINCENSE' --exclude '*.csv' \
    --exclude '*.png' --exclude '.*' --exclude 'Release' --exclude 'Debug'  \
    . $destination
#--exclude 'makefile' --exclude '*.h' --exclude '*.cpp' 

# Results
rsync "${options[@]}" --include="*.csv" --exclude="*" $destination $results

