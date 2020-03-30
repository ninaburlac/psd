#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/chn1 /nfs/gerda6/users/burlac/test/test.tier2.root 1 > /nfs/gerda6/users/burlac/psd/chn1/psd_chn1.out

chmod 666 *.*

