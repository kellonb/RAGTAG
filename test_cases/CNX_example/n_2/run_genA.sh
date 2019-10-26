#!/bin/bash

#genA=/mnt/raidc2/kbelfon/04.programs_and_scripts/genA_kb/GenA_current/genA
genA=/mnt/raidc3/kbelfon/04.programs_and_scripts/genA_kb/GenA_Aug19/genA

#submit the equilibration script
#sbatch << EOF
#!/bin/bash

#SBATCH --gres=gpu:1
#SBATCH --job-name n.${i}
#SBATCH --dependency=singleton
#SBATCH --partition="680,780,980,1080"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kellonbelfon@gmail.com
###SBATCH --nodelist="naga84"

#echo $CUDA_VISIBLE_DEVICES

$genA -i genA_input -p parmfile.4 -r CNX.restart -s CNX.scores -f CNX.frcmod -y CNX.fit -o CNX.log  

#EOF
