#!/bin/sh

#  sync.sh
#  
#
#  Created by Andrzej Chlebicki on 23.12.2018.
#  


# Sources
# rsync -r -avz --prune-empty-dirs --exclude '*.html' --exclude 'obj/' --exclude '*.csv' --exclude '*.png' --exclude '.*' --exclude 'Release' --exclude 'Debug' . ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalIntegrator
make
rsync -r -avz --prune-empty-dirs --exclude '*.nb' --exclude '*.html' --exclude 'obj' --exclude 'LINCENSE' --exclude 'makefile' --exclude '*.csv' --exclude '*.h' --exclude '*.cpp' --exclude '*.png' --exclude '.*' --exclude 'Release' --exclude 'Debug'  . ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalIntegrator

# Results
rsync -avz --include="*.csv" --exclude="*" ac357729@kruk-host.fuw.edu.pl:~/Documents/FunctionalIntegrator/ ./results

