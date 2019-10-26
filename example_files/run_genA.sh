#!/bin/bash

genA=/mnt/raidc3/kbelfon/04.programs_and_scripts/genA_kb/GenA_Aug19/genA

#submit the equilibration script
sbatch << EOF
#!/bin/bash

#SBATCH --gres=gpu:1
#SBATCH --job-name n.${i}
#SBATCH --dependency=singleton
#SBATCH --partition="680,780,980,1080"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kellonbelfon@gmail.com

echo $CUDA_VISIBLE_DEVICES

$genA -i genA_input -p parmfile.4 -r CNX_4p.restart -s CNX_4p.scores -f CNX_4p.frcmod -y CNX_4p.fit -o CNX_4p.log  
$genA -i genA_input -p parmfile.4.1 -c CNX_4p.restart -r CNX_4p.restart1 -s CNX_4p.scores1 -f CNX_4p.frcmod1 -y CNX_4p.fit1 -o CNX_4p.log1 
$genA -i genA_input -p parmfile.4.2 -c CNX_4p.restart1 -r CNX_4p.restart2 -s CNX_4p.scores2 -f CNX_4p.frcmod2 -y CNX_4p.fit2 -o CNX_4p.log2 

EOF
