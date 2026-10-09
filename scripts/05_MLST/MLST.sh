#!/usr/bin/bash
#SBATCH --job-name=mlst #作业名
#SBATCH --mem-per-cpu=2gb   #内存
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --output=%j.log
#SBATCH --partition=vip
mlst ./*.fna --scheme abaumannii_2 --csv >mlst_Ab.csv