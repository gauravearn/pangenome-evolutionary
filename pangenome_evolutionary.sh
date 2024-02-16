#! /usr/bin/bash 
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-1-25
# an end to end workflow for the complete analysis 
# of the pangenomes from the sequenced genomes given an 
# proteome fasta files and the nucleotide fasta files
# it runs on the slurms and can be easily configured to the snakemake or the nextflow
read -r -p "please provide the batch:" batch
read -r -p "please provide the queue:" queue
read -r -p "please provide the number of the cores:" core
read -r -p "please provide the channel cores:" channel
read -r -p "please provide the memory allocation:" allocation
read -r -p "please provide the workdir:" workdir
read -r -p "please provide the user mail:" mail
    echo "catching the variables for the slurm"
echo "#!/bin/bash -l"
echo "#SBATCH -J "${batch}""
echo "#SBATCH --constraint="snb"|"hsw""
echo "#SBATCH -p "${queue}""
echo "#SBATCH -n "${core}""
echo "#SBATCH -c "${channel}""
echo "#SBATCH --mem="${allocation}""
echo "#SBATCH --workdir="${workdir}""
echo "#SBATCH --mail-type=END"
echo "#SBATCH --mail-user="${mail}""
    echo "the slurm configuration is complete"
    echo "this is an automated analysis of the pangenomes sequenced from the arabidopsis genomes"
    echo "this analysis uses either the complete genomes annotations and 
                provides a pangenome analysis for the core genes"
read -r -p "please kindly select the option:" option
read -r -p "please print the alignment and the tree calibration approaches:" approaches
if [[ "${option}" == "" ]] && 
     [[ "${approaches}" == "" ]]
then 
    echo "there are two options available"
    echo "1. give the fasta file"
    echo "2. give the ortho dir"
        echo "for the alignment of the fasta files for the model calibration and the alignment"
        echo "macse & mafft & prank is available and you need to provide the tool name as the input"
        echo "for the phylogenentic reconstruction there are two tools available use iqtree or use raxml"
        echo "thank you for the selection of the approaches"
fi
exit 1
if [[ ${option} == "fasta" ]]
then 
    echo "please provide the directory path:" 
    read -r -p "please provide the directory path:" dirpath
    read -r -p "please provide the path for the macse as a part of the alignment:" macse
    read -r -p "please provide the number of the threads:" threads
    read -r -p "please provide the alignment tools:" alignment
    read -r -p "please provide the phylogeny:" phylogeny
    read -r -p "please provide the path to the nucleotide fasta files:" nucleotide
    if [[ -d ${dirpath} ]] && 
            [[ ${macse} ]] &&
                [[ ${threads} ]] &&
                    [[ ${alignment} == "macse" ]] &&
                            [[ ${phylogeny} == "iqtree" ]] && 
                                                [[ ${nucleotide} ]]
        then
            echo "directory path verified"
            cd ${dirpath}
            echo "setting up the environment variable for the analysis"
            conda create -n pangenome && conda install -n pangenome -c bioconda orthofinder
            sudo apt-get install blast2
            conda activate pangenome
            conda install -n pangenome -c bioconda trimal
            conda install -n pangenome -c bioconda muscle
            conda install -n pangenome -c bioconda prank
            conda install -n pangenome -c bioconda mafft
            conda install -n pangenome -c bioconda iqtree
            conda install -n pangenome -c bioconda orthofinder
            conda install -n pangenome -c bioconda raxml
            conda deactivate
            cd "${dirpath}"
            export PATH="${macse}":$PATH
            echo $PATH
            echo "all the required configurations have been done"
            for i in ${dirpath}/*.faa
            do 
                grep ">" -c "{i}" >> number_of_proteins.txt
            done 
            echo "running the orthofinder for the orthology assignments"
            orthofinder -t "${threads}" -a 60 -M dendroblast -S diamond -M msa \
                -A mafft -T fasttree -n "${dirpath}"_analysis -p "${dirpath} \
                      >> ortho_analysis_log.txt 
            echo "orthology analysis finished for the pangenome"
            echo "making the alignments and the single core pangenome analysis"
            cd ..
            mkdir single_core_genes
            cp -r "${dirpath}"_analysis/single_genes /single_core_genes
            for i in single_core_genes/*.faa
                do 
                    echo gep ">" "${i}" > "${i}"_protein_id.txt
                done
            for i in "${i}"_protein_id.txt
            do 
                sed -i -e "s/>/d/g" "${i}"
            done
            ls -l *.faa > protein_faa.txt
            ls -l *.txt > proteinid.txt
            paste protein_faa.txt proteinid.txt | \
               while read col1 col2 col3; do 
                pullseq -i ${col1} -n ${col2} $$ "${col3}".extraction.fasta; done
            cd ..
            mkdir nucleotide_alignments
            cp -r /single_core_genes/*.extractions.fasta ./
            for i in $(pwd)/*.fasta;do 
                echo $PATH
                java -jar -Xmx100g macse -prog alignSequences \
                             -gc_def 12 -seq $i -out_AA ${i}%.AA -out_NT ${i}%.NT > ${i}.macse.run.log.txt
            done
            for i in $(pwd)/*.NT;do 
                mv "${i}" "${i}%.macse.fasta"
            done
            for i in *.macse.fasta;do 
                trimal -in "${i}" -out "${i}%.aligned.strict.fasta" --strict
            done
            for i in *.macse.fasta;do 
                trimal -in "${i}" -out "${i}%.aligned.strict_nogaps.fasta" --strict --nogaps
            done
            for i in *.aligned.strict.fasta;do 
                statal -in "${i}" -out "${i}%.strict.stats" 
            done
            for i in *.aligned.strict_nogaps.fasta;do 
                statal -in "${i}" -out "${i}%.strict_nogaps.stats" 
            done
            for i in $(ls -l *.aligned.strict.fasta); do 
                iqtree --seqtype DNA -s "${i}" --alrt 1000 -b 1000 -T "${threads}"
            done
            for i in $(ls -l *.strict_nogaps.fasta); do 
                iqtree --seqtype DNA -s "${i}" --alrt 1000 -b 1000 -T "${threads}"
            done
fi
if [[ -d ${dirpath} ]] && 
            [[ ${macse} ]] &&
                [[ ${threads} ]] &&
                    [[ ${alignment} == "macse" ]] &&
                            [[ ${phylogeny} == "raxml" ]] && 
                                                [[ ${nucleotide} ]]
        then
            echo "directory path verified"
            cd ${dirpath}
            conda create -n pangenome && conda install -n pangenome -c bioconda orthofinder
            sudo apt-get install blast2
            conda activate pangenome
            conda install -n pangenome -c bioconda trimal
            conda install -n pangenome -c bioconda muscle
            conda install -n pangenome -c bioconda prank
            conda install -n pangenome -c bioconda mafft
            conda install -n pangenome -c bioconda iqtree
            conda install -n pangenome -c bioconda orthofinder
            conda install -n pangenome -c bioconda raxml
            conda deactivate
            cd "${dirpath}"
            export PATH="${macse}":$PATH
            echo $PATH
            echo "all the required configurations have been done"
            for i in ${dirpath}/*.faa
            do 
                grep ">" -c "{i}" >> number_of_proteins.txt
            done 
            echo "running the orthofinder for the orthology assignments"
            orthofinder -t "${threads}" -a 60 -M dendroblast -S diamond -M msa \
                            -A mafft -T fasttree -n "${dirpath}"_analysis -p "${dirpath} >> ortho_analysis_log.txt 
            echo "orthology analysis finished for the pangenome"
            echo "making the alignments and the single core pangenome analysis"
            cd ..
            mkdir single_core_genes
            cp -r "${dirpath}"_analysis/single_genes /single_core_genes
            for i in single_core_genes/*.faa;do 
                    echo gep ">" "${i}" > "${i}"_protein_id.txt
                done
            for i in "${i}"_protein_id.txt; do 
                sed -i -e "s/>/d/g" "${i}"
            done
            ls -l *.faa > protein_faa.txt && ls -l *.txt > proteinid.txt
            paste protein_faa.txt proteinid.txt | \
               while read col1 col2 col3; do 
                            pullseq -i ${col1} -n ${col2} $$ "${col3}".extraction.fasta; done
            cd ..
            mkdir nucleotide_alignments
            cp -r /single_core_genes/*.extractions.fasta ./
            for i in $(pwd)/*.fasta; do 
                echo $PATH
                java -jar -Xmx100g macse -prog alignSequences \
                             -gc_def 12 -seq $i -out_AA ${i}%.AA -out_NT ${i}%.NT > ${i}.macse.run.log.txt
            done
            for i in $(pwd)/*.NT; do 
                mv "${i}" "${i}%.macse.fasta"
            done
            for i in *.macse.fasta; do 
                trimal -in "${i}" -out "${i}%.aligned.strict.fasta" --strict
            done
            for i in *.macse.fasta; do 
                trimal -in "${i}" -out "${i}%.aligned.strict_nogaps.fasta" --strict --nogaps
            done
            for i in *.aligned.strict.fasta; do 
                statal -in "${i}" -out "${i}%.strict.stats" 
            done
            for i in *.aligned.strict_nogaps.fasta; do 
                statal -in "${i}" -out "${i}%.strict_nogaps.stats" 
            done
            for i in $(ls -l *.aligned.strict.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done
            for i in $(ls -l *.strict_nogaps.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict_nogaps -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict_nogaps -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done
            done 
fi
if [[ -d ${dirpath} ]] && 
            [[ ${macse} ]] &&
                [[ ${threads} ]] &&
                    [[ ${alignment} == "muscle" ]] &&
                            [[ ${phylogeny} == "iqtree" ]] && 
                                                [[ ${nucleotide} ]]
        then
            echo "directory path verified"
            cd ${dirpath}
            conda create -n pangenome && conda install -n pangenome -c bioconda orthofinder
            sudo apt-get install blast2
            conda activate pangenome
            conda install -n pangenome -c bioconda trimal
            conda install -n pangenome -c bioconda muscle
            conda install -n pangenome -c bioconda prank
            conda install -n pangenome -c bioconda mafft
            conda install -n pangenome -c bioconda iqtree
            conda install -n pangenome -c bioconda orthofinder
            conda install -n pangenome -c bioconda raxml
            conda deactivate
            cd "${dirpath}"
            export PATH="${macse}":$PATH
            echo $PATH
            echo "all the required configurations have been done"
            for i in ${dirpath}/*.faa
            do 
                grep ">" -c "{i}" >> number_of_proteins.txt
            done 
            echo "running the orthofinder for the orthology assignments"
            orthofinder -t "${threads}" -a 60 -M dendroblast -S diamond -M msa \
                -A mafft -T fasttree -n "${dirpath}"_analysis -p "${dirpath} \
                      >> ortho_analysis_log.txt 
            echo "orthology analysis finished for the pangenome"
            echo "making the alignments and the single core pangenome analysis"
            cd ..
            mkdir single_core_genes
            cp -r "${dirpath}"_analysis/single_genes /single_core_genes
            for i in single_core_genes/*.faa; do 
                    echo gep ">" "${i}" > "${i}"_protein_id.txt
                done
            for i in "${i}"_protein_id.txt; do 
                sed -i -e "s/>/d/g" "${i}"
            done
          for i in $(ls -l *.aligned.strict.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done
            for i in $(ls -l *.strict_nogaps.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict_nogaps -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict_nogaps -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done  ls -l *.faa > protein_faa.txt && ls -l *.txt > proteinid.txt
            paste protein_faa.txt proteinid.txt | \
               while read col1 col2 col3; do 
                            pullseq -i ${col1} -n ${col2} $$ "${col3}".extraction.fasta 
                       done
            cd ..
            mkdir nucleotide_alignments
            cp -r /single_core_genes/*.extractions.fasta ./
            for i in $(pwd)/*.fasta; do 
                echo $PATH
                muscle -in "${i}" -out "${i}%".muscle.aligned.fasta > ${i}.macse.run.log.txt
            done
            for i in *.muscle.aligned.fasta; do 
                trimal -in "${i}" -out "${i}%.aligned.strict.fasta" --strict
            done
            for i in *.muscle.aligned.fasta; do 
                trimal -in "${i}" -out "${i}%.aligned.strict_nogaps.fasta" --strict --nogaps
            done
            for i in *.aligned.strict.fasta; do 
                statal -in "${i}" -out "${i}%.strict.stats" 
            done
            for i in *.aligned.strict_nogaps.fasta; do 
                statal -in "${i}" -out "${i}%.strict_nogaps.stats" 
            done
            for i in $(ls -l *.aligned.strict.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done
            for i in $(ls -l *.strict_nogaps.fasta); do 
                echo "building the phylogeny using the GTRGAMMA model"
                echo "building the phylogeny using the GTRCAT model"
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA 
                                    -p 12345 -n phylogeny_GTRGAMMA_strict_nogaps -T "${threads}" -N 50
                raxmlHPC-PTHREADS -s "${i}" --no-seq-check -O -m GTRGAMMA -p 12345 
                                             -n phylogeny_GTRGAMMA.bootstrap -T "${threads}" -N 50 -b 1000
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O \
                                    -m GTRCAT -p 12345 -n phylogeny_GTRCAT_strict_nogaps -T 10 -N 50
                    raxmlHPC-PTHREADS -s $FILE1 --no-seq-check -O -m GTRCAT \
                                 -p 12345 -n phylogeny_GTRCAT.bootstrap -T 10 -N 50 -b 1000
            done

else 
            echo "none of the options selected and you have to restart the analysis"
                printf "%s\n" "thank you for using the pangenomes"
fi
