#!/bin/bash
#SBATCH --job-name=xaa_qc #作业名 
#SBATCH --mem-per-cpu=6gb   #内存
##SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
##SBATCH --mem=200G
#SBATCH --output=xaa_qc_%j.log
#SBATCH --partition=vip
ASM_DIR=/public/home/wangjm01/zzf/Ab_project/xaa_dir OUTDIR=qc_xaa DO_ANI=1 REF_FOR_ANI=/public/home/wangjm01/zzf/Reference_genome/GCF_009035845.1_ASM903584v1_genomic.fna ANI_MIN=95.0 bash qc_assemblies.sh