# Supplementary Methods

Commands, bioinformatic programs and their parameters used for metagenomic aand metatranscriptomic data analysis

## Metagenomics Analysis

### Quality Control

#### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```bash
fastqc *_R1_001.fastq.gz *_R2_001.fastq.gz -o ./
```

#### [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) and [samtools](https://www.htslib.org/)

```bash
#Step 1: Trim off last bases
bbduk.sh in=$sample_R1_001.fastq.gz in2=$sample_R2_001.fastq.gz out=$sample_R1_001.lastbase_rm.fastq out2=$sample_R2_001.lastbase_rm.fastq ftm=5

#Step 2: trim off the partial adapter
bbduk.sh in=$sample_R1_001.lastbase_rm.fastq in2=$sample_R2_001.lastbase_rm.fastq out=$sample_R1_001.adapter_rm.fastq out2=$sample_R2_001.adapter_rm.fastq ref=~/data/Programs/Metagenomics/bbmap/resources/adapters.fa tbo tpe k=23 mink=11 hdist=1 ktrim=r

# Step 3: filter out contaminates
bbduk.sh in=$sample_R1_001.adapter_rm.fastq in2=$sample_R2_001.adapter_rm.fastq out=$sample_R1_001.dec.fastq out2=$sample_R2_001.dec.fastq outm=$sample_contaminates.fq ref=~/data/Programs/Metagenomics/bbmap/resources/phix_adapters.fa.gz k=31 hdist=1 stats=$sample_stats.txt

# Step 4: clipping off the low quality ends and low complexity regions (high entropy; AAAAAAA)
bbduk.sh in=$sample_R1_001.dec.fastq in2=$sample_R2_001.dec.fastq out=$sample_R1_001.qc.fastq out2=$sample_R2_001.qc.fastq qtrim=rl trimq=15 minlength=30 entropy=0.5
```

### Assembly

#### [MegaHit](https://github.com/voutcn/megahit)

```bash
#coassembly
megahit -1 CAT_ALL_R1_001.qc.fastq -2 CAT_ALL_R2_001.qc.fastq -o output/dir

#individual assemblies
megahit -1 sample_R1_001.qc.fastq -2 sample_R2_001.qc.fastq -o output/dir/S52

```

#### [MetaViralSpades](https://github.com/ablab/spades)

```bash
spades.py --metaviral -t 30 --pe1-1 CAT_S64_R1_001.qc.fastq --pe1-2 CAT_S64_R2_001.qc.fastq --pe2-1 CAT_S65_R1_001.qc.fastq --pe2-2 CAT_S65_R2_001.qc.fastq --pe3-1 CAT_S66_R1_001.qc.fastq --pe3-2 CAT_S66_R2_001.qc.fastq --pe4-1 CAT_S70_R1_001.qc.fastq --pe4-2 CAT_S70_R2_001.qc.fastq -o ./S64_65_66_70_viralSpades/
```

#### [MetaQUAST](https://quast.sourceforge.net/metaquast.html)

```bash
metaquast.py -m 500 -t 20 --fast -o metaquast -l viral contigs.fasta
```

### Mapping

#### [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/) & [samtools](https://www.htslib.org/) & [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

```bash
bbmap.sh ref= in=sample_R1_001.qc.fastq in2=sample_R2_001.qc.fastq out=./sample_bbmap.sam covstats=./sample_bbmap_covstats.txt scafstats=./sample_bbmap_scafstats.txt threads=20 minid=0.95 ambiguous=toss

samtools view -b ./sample_bbmap.sam | samtools sort -o ./sample_bbmap_sorted.bam

samtools index ./sample_bbmap_sorted.bam

jgi_summarize_bam_contig_depths --outputDepth ./depth.txt ./*.bam
```

### Binning

Binning from coassembly into Metagenome-Assembled-Genomes (MAGs) using [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

```bash
metabat2 -i final.contigs.fa.gz -o /Bin -a depth.txt --unbinned -t 30 -v
```

Individual binning of Ca. Sodalinema alkaliphilum MAGs

#### [MetaQUAST](https://quast.sourceforge.net/metaquast.html)

```bash
metaquast.py -m 500 -t 20 -o metaquast_full -l day0,day2,day4,day6,day8,day9 day0/final.contigs.fa day2/final.contigs.fa day4/final.contigs.fa day6/final.contigs.fa day8/final.contigs.fa day9/final.contigs.fa -r phormidium_ref.fasta
```

### MAG Quality

#### [CheckM2](https://github.com/chklovski/CheckM2)

```bash
checkm2 predict -t 30 -x fa --input ./allBins/bins/ --output-directory ./Checkm2
```

### Taxonomy

#### [Phyloflash](https://hrgv.github.io/phyloFlash/)

```bash
phyloFlash.pl -dbhome silva/138.1 -lib samples  -read1 sample_R1_001.qc.fastq -read2 sample_R2_001.qc.fastq -readlength 150
```

#### [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

Genome Taxonomy Database release 202 was used to classify MAG taxonomy

```bash
gtdbtk classify_wf -x fa --cpus 20 --genome_dir ./allBins/ --out_dir ./ --pplacer_cpus 20
```

### MAG Annotation

#### [MetaErg2](https://github.com/kinestetika/MetaErg)

```bash
apptainer exec -B /work/ebg_lab/ ~/metaerg_latest.sif metaerg --database_dir /referenceDatabases/metaerg_database/ --path_to_signalp /referenceDatabases/metaerg_database/ --path_to_tmhmm /referenceDatabases/metaerg_database/ --contig_file ./allBins/ --rename_genomes --rename_contigs --cpus 40 --file_extension .fna
```

### MAG and Contig Coverage and Relative Abundance

#### [CoverM](https://github.com/wwood/CoverM)

```bash
coverm genome -v -x fa -t 20 --methods trimmed_mean -1 *_R1_001.qc.fastq -2 *_R2_001.qc.fastq --genome-fasta-directory bins --bam-file-cache-directory ./coverm_bam -o out_genome_coverage.tsv

#for relative abundance
coverm genome -v -x fa -t 20 --methods relative_abundance --bam-files *.bam -s '~' -o output_relative_abundances.tsv
```

### CRISPR Identification

#### [MinCED](https://github.com/ctSkennerton/minced)

```bash
minced -spacers final.contigs.fa assembly.crisprs assembly.gff
```

### Viral Sequence Identification

#### [VirSorter2](https://github.com/jiarong/VirSorter2)

```bash
virsorter config --set HMMSEARCH_THREADS=30

virsorter run --keep-original-seq --min-length 1000 --min-score 0.5 -w ./pass1 -i ../final.contigs.fa.gz --include-groups dsDNAphage,ssDNA -j 30 all
```

### Viral Host Predictions

#### [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

```bash
blastn -db final.contigs.fa  -query viral_genomes.fna -out blastn.out
```

**CRISPR Spacers**

```bash
blastn -db viral_genomes.fna  -query ../assembly_spacers.fa -out blastn.out -task "blastn-short" -outfmt 6 -dust no
```
CRISPR-Blastn results were filtered by length >= 25 bp, percent identity >= 95% and e-value <= 1E-5 

**tetraucleotide frequency profiles**

Using an in-house script developed by [Marc Strous](https://github.com/kinestetika) 

#### [iPHoP](https://bitbucket.org/srouxjgi/iphop/src/main/iphop/)

```bash
iphop predict --fa_file ../viral_contigs_renamed.singleline.fna --db_dir /work/ebg_lab/referenceDatabases/iphop/db/ --out_dir iphop_output/
```

### Viral Taxonomy

An [in-house script](https://github.com/vmkhot/labjournal/blob/main/Scripts/Python/BSR_tree.py) was used to predict viral taxonomy against the viral reference sequence database - Refseq Release 203.

```bash
blastp -db cat_all_viruses2blast.faa -query cat_all_viruses2blast.faa -out VCproteins_refseq_selfblastp_hsp1.out -outfmt 6 -num_threads 20 -max_hsps 1
```
Blastp results were filtered by e-value <= 1E-5 and percent identity >= 40%

Script used to calculate Dice coefficients and create a pairwise alignment of viruses: [BSR_tree.py](https://github.com/vmkhot/labjournal/blob/main/Scripts/Python/BSR_tree.py)

The resulting network was visualized in Cytoscape (ref)

### Viral Sequence Annotation

#### [VirSorter2](https://github.com/jiarong/VirSorter2) & [DRAM-V](https://github.com/WrightonLabCSU/DRAM)

```bash
virsorter run --prep-for-dramv --keep-original-seq --min-length 1000 --min-score 0.5 -w ./pass1 -i viral_genomes.fasta --include-groups dsDNAphage,ssDNA -j 50 all

#Annotate
DRAM-v.py annotate -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab -o ./dramv/ --threads 30
```

## Metatranscriptome Analysis

### Quality Control

#### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```bash
fastqc -o ./ -t 30 ../raw_fastq/*.fastq.gz

multiqc -o ./ -n rawReads ./
```

Last base and quality trimming using [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

```bash
#trim off last base (round to 50bp)
bbduk.sh in=R1_001.fastq.gz in2=R2_001.fastq.gz out=R1_001.lastbase_rm.fastq.gz out2=R2_001.lastbase_rm.fastq.gz ftm=5

#quality filtering (for small RNA, minimum read length of 10 and min quality of 15)
bbduk.sh in=R1_001.lastbase_rm.fastq.gz in2=R2_001.lastbase_rm.fastq.gz  out=R1_001.qc.fastq.gz out2=R2_001.qc.fastq.gz qtrim=rl trimq=15 minlength=10
```

#### Sorting rRNA using [SortMeRNA](https://github.com/sortmerna/sortmerna)

```bash
sortmerna --ref referenceDatabases/sortmerna_db/smr_v4.3_default_db.fasta \
       --workdir ./sortmerna/ \
       --reads R1_001.qc.fastq.gz --reads R2_001.qc.fastq.gz \
       --aligned rRNA_reads/rRNA_reads.qc \
       --other non_rRNA_reads/non_rRNA_reads.qc \
       --sam --SQ --log --fastx --threads 40 --paired_in
```

### Mapping

Mapping transcriptome reads to nucleotide gene sequences from metagenome using ["seal.sh"](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/seal-guide/) from BBTools

```bash
seal.sh in=reads.fq nucl_seq.fa stats=sealstats.txt rpkm=sealrpkm.txt ambig=all
```

### Differential Expression Analysis

DE analysis was performed using the [*DESeq2*](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)  pipeline in R. 

The generic pipeline used to assess quality of raw counts from samples: [deseq2_script_sample_QC](https://github.com/vmkhot/Metatranscriptomics/blob/main/R-scripts/deseq2_script_sample_QC.R) 

R script for using *DESeq2* for time-course experiments for the *Ca. S. alkaliphilum*: [deseq2_time_course_script_cyano](https://github.com/vmkhot/Metatranscriptomics/blob/main/R-scripts/deseq2_time_course_script_cyano.R)

### Soft Clustering of Gene Expression Profiles

Genes were clustered by their expression profiles over time using the [*MFuzz*](http://mfuzz.sysbiolab.eu/) package in R

R script for using *MFuzz* for time-course experiments for the *Ca. S. alkaliphilum*: [mfuzz_clustering](https://github.com/vmkhot/Metatranscriptomics/blob/main/R-scripts/mfuzz_clustering.R)
