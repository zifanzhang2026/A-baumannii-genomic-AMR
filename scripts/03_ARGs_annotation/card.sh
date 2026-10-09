#!/bin/bash
#SBATCH --job-name=card #作业名
#SBATCH --mem-per-cpu=2gb   #内存
#SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --output=%j.log
#SBATCH --partition=vip
abricate --db card --mincov 90 --minid 80  ./*.fna > card.result.tab
abricate --summary card.result.tab > card.summary.tab