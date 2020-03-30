#!/bin/sh

source /afs/lngs.infn.it/experiment/gerda/MJsetupDev.sh
#source /nfs/gerda6/users/dandrea/DPLMS/MJsetupVal.sh

/nfs/gerda6/users/burlac/psd/processPSD /nfs/gerda6/users/burlac/psd/analysis/chn0 /nfs/gerda6/users/burlac/psd/analysis/run0095-cal.tier.root 0 > /nfs/gerda6/users/burlac/psd/analysis/chn0/psd_chn0.out

chmod 666 *.*

