#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/chn31 /nfs/gerda6/users/burlac/test/test.tier2.root 31 > /nfs/gerda6/users/burlac/psd/chn31/psd_chn31.out

chmod 666 *.*

