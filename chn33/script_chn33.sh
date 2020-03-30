#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/chn33 /nfs/gerda6/users/burlac/test/test.tier2.root 33 > /nfs/gerda6/users/burlac/psd/chn33/psd_chn33.out

chmod 666 *.*

