#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/chn18 /nfs/gerda6/users/burlac/test/test.tier2.root 18 > /nfs/gerda6/users/burlac/psd/chn18/psd_chn18.out

chmod 666 *.*

