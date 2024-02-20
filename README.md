# pangenome_evolutionary_analysis
a complete workflow for analyzing the pangenomes from the core genesets. simply have to provide the fasta files and it will do everything and will make all the accessory information plots from the evolutionary analysis. It will also check for the breakage in the phylogeny and also will perform the repoint analysis. A ruby on rails application is also written which will single point take the slurm address and will give you the graphical analysis for all the evolutionary pangenome. There were some variable declaration to be fixed and if you use this then change the variables accrodingly and if not then check the updated version below with more parameters and support for the other analysis.  

2024-2-20 final release: Adding the supporting for the mixed linear modelling of the sequences and also for the supermatrix creation and following the phylogeny runs using the GTRCAT and GTRGAMMA phylogeny models.  An update fixing all the variable paths and adding support for the protein based as well as the nuceltodie based phylogenies and pangenomics. Made the code much shorter and within code, added support for the AWK filtering, so that external tools are not required.
```
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
```
it then loops over the multiple variables at once for the faster iterations. 
```
 for i in *.nucl.fasta; do
            cat ${i%%.*}.format.ids.short.txt | while read line; \
                    do grep -A 2 $line ${i%%.*}.format.fasta >>${i%%.*}.select.fasta; done
        done
```
Gaurav Sablok, \
Academic Staff Member, \
Bioinformatics, \
Institute for Biochemistry and Biology, \
University of Potsdam, \
Potsdam,Germany
