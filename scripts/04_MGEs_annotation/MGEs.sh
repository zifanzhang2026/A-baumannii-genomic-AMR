#!/bin/bash
#SBATCH --job-name=xaa #作业名 
#SBATCH --mem-per-cpu=4gb   #内存
##SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
##SBATCH --mem=8G
#SBATCH --output=xaa_%j.log
#SBATCH --partition=vip
#SBATCH --thread-spec=112
#SBATCH --bb=1000
for i in xaa/*.fna; do
    db_name=$(basename "$i" | cut -d '.' -f 1)
    mkdir -p "./$db_name"
    
    # 创建BLAST数据库
    makeblastdb -in "$i" -dbtype nucl -out "./$db_name/${db_name}_db"

    # 运行BLASTn查询
    blastn -query MGEs_FINAL_99perc_trim.fasta -db "./$db_name/${db_name}_db" -outfmt 6 -perc_identity 90 -qcov_hsp_perc 85 -out "./$db_name/${db_name}_mge.txt"
done