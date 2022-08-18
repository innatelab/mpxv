#! /bin/sh
#SBATCH -o /gpfs/scratch/pn69ha/ge54heq2/logs/mpxv/%x_%j.log
#SBATCH -J msglm_mpxv_fp
#SBATCH --get-user-env

#SBATCH --clusters=cm2_tiny
#SBATCH --nodes=1-4

##SBATCH --clusters=cm2
##SBATCH --partition=cm2_std
##SBATCH --qos=cm2_std
##SBATCH --nodes=3-24

##SBATCH --clusters=cm2
##SBATCH --partition=cm2_large
##SBATCH --qos=cm2_large
##SBATCH --nodes=25-32

#SBATCH --ntasks-per-node=4
#SBATCH --mincpus=7
##SBATCH --mem-per-cpu=2GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=7
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=yiqi.huang@tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

module load slurm_setup
module load charliecloud

IMAGES_PATH=$(realpath $HOME)
CHUDIS_PATH=$HOME/projects/cool_chunk_dispatcher

PROJECT_ID=mpxv
DATA_VERSION=20220817
FIT_VERSION=20220817

srun --wait=0 --no-kill --distribution=block --exclusive=user \
     -o $SCRATCH/logs/$PROJECT_ID/%x_%j_%t.log \
${CHUDIS_PATH}/slurmstep_process_chunks.sh "${PROJECT_ID}_fp_${FIT_VERSION}" $USER \
"ch-run $IMAGES_PATH/archpc.msglm \
  -t --no-home --unset-env='*PATH' \
  --set-env=$HOME/projects/adhoc/$PROJECT_ID/archpc.env \
  -b $HOME/projects/adhoc:/projects/adhoc \
  -b $HOME/data:/data \
  -b $HOME/analysis:/analysis \
  -b $SCRATCH:/scratch \
  -- Rscript /projects/adhoc/$PROJECT_ID/msglm_fit_chunk_fp.R \
  $PROJECT_ID $SLURM_JOB_NAME $DATA_VERSION $FIT_VERSION"

