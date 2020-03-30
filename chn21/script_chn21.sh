#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/chn21 /nfs/gerda6/users/burlac/test/test.tier2.root 21 > /nfs/gerda6/users/burlac/psd/chn21/psd_chn21.out

chmod 666 *.*

