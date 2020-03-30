#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/readTier /nfs/gerda6/users/burlac/psd/analysis 6.13 98 > /nfs/gerda6/users/burlac/psd/analysis/read.out

chmod 666 *.*

