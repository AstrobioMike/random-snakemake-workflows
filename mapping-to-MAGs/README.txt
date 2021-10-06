Starts with a directory full of individual MAG files, `cat`s them together into one fasta, 
makes bowtie2 index and maps all specified samples to it, then gets contig-level coverage 
information (and produces a detection-filtered version too).

The config.yaml needs to have some things set, and there needs to be input files holding unique sample IDs and unique MAG IDs.


If MAGs come from different assemblies and therefore may share contig names, we can prepend the
MAG ID to the header names if wanted with `bit-rename-fasta-headers`, e.g.:

for MAG in $(cat workflow/unique-MAG-IDs.txt)
do 
    echo ${MAG}
    bit-rename-fasta-headers -w ${MAG} -i dRepd-MAG-fastas/${MAG}.fa --prefix -o dRepd-MAG-fastas/tmp.fa
    mv dRepd-MAG-fastas/tmp.fa dRepd-MAG-fastas/${MAG}.fa
done
 
bit: https://github.com/AstrobioMike/bit#bioinformatics-tools-bit
