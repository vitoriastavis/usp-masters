instalar sra conforme
https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

prefetch

fastq-dump --split-files --gzip SRR19374522

Se não precisar de tudo, use fasterq-dump com limite:
fasterq-dump --split-files --progress SRR19374522 -O output_folder
Caso tenha problemas de espaço, pode usar a opção --concatenate-reads para combinar as leituras.

