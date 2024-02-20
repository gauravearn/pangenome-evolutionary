#! /usr/bin/bash
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-1-25
# updated 2024-2-15
# final update 2024-2-20
# an end to end workflow for the complete analysis
# of the pangenomes from the sequenced genomes given an
# proteome fasta files and the nucleotide fasta files
# it runs on the slurms and can be easily configured to the snakemake 
# or the nextflowecho "writing the slurm configurator"
echo "adding the support for the name process:"
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
    [[ "${approaches}" == "" ]]; then
    echo "there are two options available"
    echo " your protein and the nucletoide files should be formatted as the ones for the >gi|ID123456789"
    echo "if they are not formatted according to this criteria then please run the following"
    echo "for i in *.fa; do sed -i -e s"/ / | /"g ${i}; done"
    echo "for i in *.fa; do cut -f 1,2 -d " | " "${i}" > "${i%.*}".format.fasta; done"
    echo "1. give the fasta file"
    echo "2. give the ortho dir"
    echo "for the alignment of the fasta files for the model calibration and the alignment"
    echo "macse & mafft & prank is available and you need to provide the tool name as the input"
    echo "for the phylogenentic reconstruction there are two tools available use iqtree or use raxml"
    echo "thank you for the selection of the approaches"
fi
exit 1
########################################################
if [[ ${option} == "fasta" ]]
then
    echo "please provide the directory path"
    read -r -p "please provide the directory path:" dirpath
    read -r -p "please provide the path for the macse:" macse
    read -r -p "please provide the number of the threads:" threads
    read -r -p "please provide the path to the nucleotide fasta files:" nucleotide
    if [[ -d "${dirpath}" && "${macse}" && "${nucleotide}" ]]
    then
        echo "directory path verified"
        cd "${dirpath}"
        echo "setting up the environment variable for the analysis"
        conda create -n pangenome && conda install -n pangenome -c bioconda orthofinder
        conda install -n pangenome -c bioconda trimal blast2 diamond muscle prank mafft iqtree raxml -y
        conda activate pangenome
        cd "${dirpath}"
        export PATH="${macse}":$PATH
        echo $PATH
        echo "all the required configurations have been done"
        for i in ${dirpath}/*.faa; do
            grep ">" -c "{i}" >>number_of_proteins.txt
        done
        echo "formatting the headers for the alignments"
        for i in "${dirpath}"/*.faa; do
            awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  \
                                                   END {printf("\n");}' "${i}" >"${i%.*}".protein.fasta
            rm -rf *.faa
        done
        echo "formatting the headers for the super matrix construction"
        for i in "${nucleotide}"/*.fasta; do
            awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  \
                                                          END {printf("\n");}' "${i}" >"${i%.*}".nucl.fasta
            rm -rf *.fasta
        done
        echo "running the orthofinder for the orthology assignments"
        orthofinder -t "${threads}" -a 60 -M dendroblast -S diamond -M msa \
            -A mafft -T fasttree -n "${dirpath}"_analysis -f "${dirpath}" \
            >>ortho_analysis_log.txt
        echo "orthology analysis finished for the pangenome"
        echo "making the alignments and the single core pangenome analysis"
        cd ..
        mkdir single_core_genes
        cp -r "${dirpath}"_analysis/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/*.fa /single_core_genes
        for i in single_core_genes/*.fa; do
            grep ">" "${i}" >"${i}".ids.txt
        done
        cp -r *_protein_id.txt "${nucleotide}"
        cd "${nucleotide}"
        for i in *.ids.txt; do
            cut -f 2 -d "|" $i >"${i%.*}".short.txt
        done
        for i in *.nucl.fasta; do
            cat ${i%%.*}.format.ids.short.txt | while read line; \
                     do grep -A 2 $line ${i%%.*}.format.fasta >>${i%%.*}.select.fasta; done
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        for j in *.select.fasta; do
            java -jar -Xmx100g macse -prog alignSequences \
                -gc_def 12 -seq "${j}" -out_AA "${j}%".AA -out_NT \
                "${j}%".NT >"${j}".macse.run.log.txt
        done
        for i in *.NT; do
            mv "${i}" "${i%.*}".ntaligned.fasta
        done
        echo "renamed the aligned files as ntaligned"
        for i in *.ntaligned.fasta; do
            trimal -in "${i}" -out "${i%.*}".trimmed.fasta -nogaps
        done
        for i in *.ntaligned.fasta; do
            trimal -in "${i}" -out "${i%%.*}".trimmedstrict.fasta --strict
        done
            wget https://github.com/marekborowiec/AMAS/AMAS.py
            chmod 755 AMAS.py
            python3 AMAS.py -in *.trimmed.fasta -f fasta -d dna
            mv concatenated.out macsealignmentconcatenated.fasta
            mv partitions.txt macsealignmentpartitions.txt
            python3 AMAS.py -in *.trimmedstrict.fasta -f fasta -d dna
            mv concatenated.out macsealignmentconcatenated_strict.fasta
            mv partitions.txt macsealignmentpartitions_strict.txt
        iqtree --seqtype DNA -s macsealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
            -p 12345 -n macsephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s macsealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
            -p 12345 -n macsephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        iqtree --seqtype DNA -s macsealignmentconcatenated_strict.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s macsealignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
            -p 12345 -n macsephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s macsealignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
            -p 12345 -n macsephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        echo "finishing up the analysis"
    if [[ -d ${dirpath} ]] && ${nucleotide} && ${macse} == "" ]]
    then
        echo "directory path verified"
        cd "${dirpath}"
        echo "setting up the environment variable for the analysis"
        conda create -n pangenome && conda install -n pangenome -c bioconda orthofinder
        conda install -n pangenome -c bioconda trimal blast2 diamond muscle prank mafft iqtree raxml -y
        conda activate pangenome
        cd "${dirpath}"
        export PATH="${macse}":$PATH
        echo $PATH
        echo "all the required configurations have been done"
        for i in ${dirpath}/*.faa; do
            grep ">" -c "{i}" >>number_of_proteins.txt
        done
        echo "formatting the headers for the alignments"
        for i in "${dirpath}"/*.faa; do
            awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  \
                                                END {printf("\n");}' "${i}" >"${i%.*}".protein.fasta
            rm -rf *.faa
        done
        echo "formatting the headers for the super matrix construction"
        for i in "${nucleotide}"/*.fasta; do
            awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  \
                                                END {printf("\n");}' "${i}" >"${i%.*}".nucl.fasta
            rm -rf *.fasta
        done
        echo "running the orthofinder for the orthology assignments"
        orthofinder -t "${threads}" -a 60 -M dendroblast -S diamond -M msa \
                                -A mafft -T fasttree -n "${dirpath}"_analysis -f "${dirpath}" \
                                                                          >>ortho_analysis_log.txt
        echo "orthology analysis finished for the pangenome"
        echo "making the alignments and the single core pangenome analysis"
        cd ..
        mkdir single_core_genes
        cp -r "${dirpath}"_analysis/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/*.fa /single_core_genes
        for i in single_core_genes/*.fa; do
            grep ">" "${i}" >"${i}".ids.txt
        done
        cp -r *_protein_id.txt "${nucleotide}"
        cd "${nucleotide}"
        for i in *.ids.txt; do
             cut -f 2 -d "|" $i >"${i%.*}".short.txt
        done
        for i in *.nucl.fasta; do
            cat ${i%%.*}.format.ids.short.txt | while read line; \
                    do grep -A 2 $line ${i%%.*}.format.fasta >>${i%%.*}.select.fasta; done
        done
        echo "aligning the pangenome using the prank probabilistic alignment"
        for i in *.select.fasta; do
            p    rank -d "${i}" -o "${i%.*}".prankaligned.fasta
        done
        for i in *.prankaligned.fasta; do
             trimal -in "${i}" -out "${i%.*}".pranktrimmed.fasta -nogaps
             trimal -in "${i}" -out "${i%.*}".pranktrimmed_strict.fasta -strict
        done
            wget https://github.com/marekborowiec/AMAS/AMAS.py
            chmod 755 AMAS.py
            python3 AMAS.py -in *.pranktrimmed.fasta -f fasta -d dna
            mv concatenated.out prankalignmentconcatenated.fasta
            mv partitions.txt prankalignmentpartitions.txt
            python3 AMAS.py -in *.pranktrimmed_strict.fasta -f fasta -d dna
            mv concatenated.out prankalignmentconcatenated_strict.fasta
            mv partitions.txt prankalignmentpartitions_strict.txt
        iqtree --seqtype DNA -s prankalignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                                 -p 12345 -n prankphylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s prankalignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                            -p 12345 -n prankphylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
           iqtree --seqtype DNA -s prankalignmentconcatenated_strict.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s prankalignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
                                                    -p 12345 -n prankphylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s prankalignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
                                        -p 12345 -n prankphylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        echo "aligning the pangenome using the muscle alignment"
        for i in *.select.fasta; do
                muscle -in "${i}" -out "${i%.*}".musclealigned.fasta
        done
        for i in *.musclealigned.fasta; do
            trimal -in "${i}" -out "${i%.*}".muscletrimmed.fasta -nogaps
            trimal -in "${i}" -out "${i%.*}".muscletrimmed_strict.fasta -strict
        done
        wget https://github.com/marekborowiec/AMAS/AMAS.py
        chmod 755 AMAS.py
        python3 AMAS.py -in *.muscletrimmed.fasta -f fasta -d dna
        mv concatenated.out musclealignmentconcatenated.fasta
        mv partitions.txt musclealignmentpartitions.txt
        python3 AMAS.py -in *.muscletrimmed_strict.fasta -f fasta -d dna
        mv concatenated.out musclealignmentconcatenated_strict.fasta
        mv partitions.txt musclealignmentpartitions_strict.txt
        iqtree --seqtype DNA -s musclealignmentconcatenated.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s musclealignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                    -p 12345 -n musclephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s musclelignmentconcatenated.fasta --no-seq-check -O -m GTRGAMMA \
                                -p 12345 -n musclephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
        iqtree --seqtype DNA -s musclealignmentconcatenated_strict.fasta --alrt 1000 -b 1000 -T "${threads}"
        raxmlHPC-PTHREADS -s musclealignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
                                    -p 12345 -n musclephylogeny_GAMMA -T "${threads}" -N 50
        raxmlHPC-PTHREADS -s musclelignmentconcatenated_strict.fasta --no-seq-check -O -m GTRGAMMA \
                                            -p 12345 -n musclephylogeny_GTRCAT -T "${threads}" -N 50 -b 1000
    fi
fi 
