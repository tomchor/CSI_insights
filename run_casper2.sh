#!/bin/bash -l
name=TEST2
sbatch <<EOT
#!/bin/bash -l
#SBATCH --account=UMCP0012
#SBATCH --job-name=$name
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=0:00:10
#SBATCH --nodes=1
#SBATCH --gres=gpu:v100:1 # GPUs per Node
#SBATCH --mail-user=tchor@umd.edu
#SBATCH --mail-type=ALL

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda

echo SOMETHING$name
EOT
