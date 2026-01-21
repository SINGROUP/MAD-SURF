#!/bin/bash
#SBATCH --account=project_2008059
#SBATCH --partition=gpusmall #gpusmall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:15:00
#SBATCH --gres=gpu:a100:1

#module load gcc/11.2.0 openmpi/4.1.6 scicomp-python-env
#source activate boss_nequip
#source activate /scratch/work/gonzalm9/conda_envs/boss_nequip
#source /scratch/project_2012660/venv_ase_flare_mace/ase_flare_mace/bin/activate
#export PATH="/scratch/project_2012660/mace_env/bin:$PATH" 
#export PATH="/scratch/project_2012660/mace_env_cueq/bin:$PATH"
export PATH="/scratch/project_2012660/mace_faiss_filter/bin:$PATH"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

echo "GPU is available/Torch version:"
python3 -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'
#python3 -c 'import faiss; print(faiss.__version__)'

python descriptor_filter_batch_faiss.py \
  --input /scratch/project_2012660/MACE_training/descriptor_filtering/full_set/full_train_test_std_config_types.extxyz \
  --model /scratch/project_2012660/MACE_training/full_dataset_train/MACE_model.model \
  --output filtered_structures.extxyz \
  --threshold 0.01 \
  --device cuda \
  --nn_index 320 \
  --nn_k 5 \
  --plot_distances \
  --hist_output my_histogram.png

python save_extxyz_frames_train_test_split.py filtered_structures.extxyz

