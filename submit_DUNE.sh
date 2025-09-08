#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00 
#SBATCH --mem=50G

# at seawulf
#source /gpfs/home/uyevarouskay/Work/env_gundam_home_v4.sh;

# at SBU NN home cluster
#source /home/uyevarou/Work/env_dune_lbl.sh
source /home/isgould/work/gundam/install/setup.sh

# at seawulf cluster
#export OA_INPUTS=/gpfs/scratch/uyevarouskay/DUNE/OA-inputs/lbl/gudi-inputs/v16/
#export NUOSCILLATOR_ROOT_LIB=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
#source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

# at SBU NN home cluster
export OA_INPUTS=/storage/shared/DUNE/OA-inputs/lbl/gudi-inputs/v5/
export NUOSCILLATOR_ROOT_LIB=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

# if you request an Asimov fit, simply disregard the `-d` option, --scan is to build llh scans
gundamFitter -a -d --scan -c ./config_DUNE.yaml -t 8
# To disable  oscillation parameters consider override the parameters config:
#gundamFitter -a -d --scan -c ./config_DUNE.yaml -of overrides/disableOscillationParameters.yaml -t 8

