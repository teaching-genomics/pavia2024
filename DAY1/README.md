# BGPAV 2024 Day1

## Installation (use this if you want to replicate the practical but installing the tools yourself in your machine)

### Install needed software 

Let's install miniconda: 
1. go to https://docs.conda.io/projects/conda/en/latest/user-guide/install/
2. under regular installation, choose your operating system (Mac, Linux)
3. download the appropriate installer (also directly from this page https://docs.conda.io/en/latest/miniconda.html)
4. select the appropriate CPU architecture (32/64 bit), if you are unsure you can check out this page https://www.computerhope.com/issues/ch001121.htm
5. once the download is completed, move the installer in your workspace (the location of downloaded files may vary!):
```
mv ../../Downloads/Miniconda3-latest-MacOSX-x86_64.sh .
```
6. install conda
```
sh Miniconda3-latest-MacOSX-x86_64.sh
```
7. accept license with 'yes'. You can use the default install location in the prompt. At the end of the installation process, close the terminal and open a new window, check that conda is working
```
conda
``` 

#### Note 
most new computers will be 64! For Mac make sure you pick the right installer, ending with 'bash' (e.g. Miniconda3 MacOSX 64-bit bash)


### Install bowtie2
Let's install a NGS read aligner (bowtie2)
```
conda install -c bioconda bowtie2
```

### Install sra-tools
Let's install the dedicated tool from NCBI-SRA for downloading data
```
conda install -c bioconda sra-tools=3.0
```

### Install fastqc
Install and check it works by exploring the help page:
```
conda install -c bioconda fastqc
fastqc -h | less
```

### Install samtools
Let's install the software samtools
```
conda install -c bioconda "samtools>=1.10"
```

### Install freebayes
Let's install the freebayes variant caller
```
conda install bioconda::freebayes
```

### Let's install vcftools and bcftools, both very useful and fast programs for handling vcf files
```
conda install -c bioconda vcftools bcftools
```

### Let's install rasusa tool
```
conda install -c bioconda rasusa
```

### Let's install minimap2 
```
conda install biobuilds::minimap2
```

### Let's install blast
```
conda install -c bioconda blast
```

### Let's install hifiasm
```
conda install -c bioconda hifiasm
```

## Practical in gitpod (use this during the course practical; no need to install anything)


## Reference genome
1. go to https://www.ncbi.nlm.nih.gov/
2. search for 'e coli' (without quotes) 
3. under the 'Genomes' tab, click 'Genome'
4. select the first result
5. under 'NCBI RefSeq assembly' > RefSeq click on Actions--> see more files on ftp
6. Right click on the file 'GCF_000005845.2_ASM584v2_genomic.fna.gz' and select "copia indirizzo link" / "copy link address"


Retrieve E. coli reference genome (wget command and then paste the adress you just copied) in the gitpod workspace

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
```

Gunzip the genome and rename it for brevity
```
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna e_coli.fasta
```

Let's use gfastast to evaluate the reference 

```
gfastats e_coli.fasta 4641652
```
What information can you get about the assembly?

## Alignment
Now let's download some short reads to align to the reference genome we just downloaded.  
The idea is that everyone should pick a different sample
```
fasterq-dump SRR10428014  #fasterq-dump is a tool to download seequences from NCBI-SRA
```
Let's explore the fastq format. How can we visualize the first read?
```
cat SRR10428014_1.fastq | head -4
```
How many reads do we have in total? Let's count the number of lines 
```
cat SRR10428014_1.fastq | wc -l
```
Output example: 5270084
And divide by four:
```
echo $(( 5270084 / 4 ))
```
We will assess raw read quality using rdeval, a very useful software to analyse raw data.

```
rdeval SRR10428014_1.fastq 4641652   #the command needs the genome size in bp
rdeval SRR10428014_2.fastq 4641652   #the command needs the genome size in bp
```
Anayze the output trying to interpret the results according to what you learned during the theoretical lessons about NGS. 
And now let's align the reads to our reference. First we need to index the reference to make it faster to access:
```
bowtie2-build e_coli.fasta ecoli  #The command should print many lines of output then quit. When the command completes, the current directory will contain four new files that all start with 'ecoli' and end with .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. These files constitute the index - you're done!

```
Then we can align the reads:
```
bowtie2 -x ecoli -1 SRR10428014_1.fastq -2 SRR10428014_2.fastq -p 4 > alignment.sam
```
Now let's convert our alignment into a format that IGV can read (bam):
```
samtools view -b alignment.sam > alignment.bam
```
Let's sort by coordinate
```
samtools sort alignment.bam -o sorted_alignment.bam -@ 4
```
Let's generate the index for this file
```
samtools index sorted_alignment.bam
rm alignment.sam # we don't need the sam anymore
```

Let's also index the reference fasta file using samtools, so that the fasta can be loaded in IGV, with the command samtools faidx
```
samtools faidx e_coli.fasta
```

## Visualization in IGV
Let's download IGV viewer to look inside the alignment https://software.broadinstitute.org/software/igv/download (download the version with java included; choose a version appropriate for your OS)
There is also an online version in case you don't want to download it locally (https://igv.org/app/)
Run IGV and load the files to inspect them. You can also load the genome annotation as an additional track in the GFF format (available here https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz)


## Variant calling
Let's check that freebayes, the variant caller we will use, works by opening the help page
```
freebayes -h
```
Let's perform the variant calling using the alignment generated before
```
freebayes -f e_coli.fasta -p 1 sorted_alignment.bam > var.vcf
```
Now that you have the variant calling file from your alignment, you can try to visualize it in IGV by loading it as you did before for the alignment file.

## How to change sample name in the vcf file
Normally, sample name is assigned during alignment, where you can assign a sample ID for the read groups deriving from the same sample.  
Since we skipped this step of sample assignment in our alignment, for the tutorial purposes you can just assign a sample name later in the vcf using bcftools reheader like this:  
```
bcftools reheader --samples new_sample_name -o name_surname.vcf.gz var.vcf.gz  ### new_sample_name is a simple text file containing the new sample name to give (you can use your name_surname; you can create it using vim/nano as showed in classroom or also in a notebook and then transferring it to gitpod workspace); -o defines the name of the output file that will be generated (with the new sample name); the last part of the command is the input file (the vcf file obtained before using freebayes)
```


## Now let's learn how to explore a vcf file using the command line
Let's visualize the header of the file and learn some useful commands  
This command shows the file header; the header contains information on what has been done to the vcf
```
bcftools view -h var.vcf  
bcftools view -h var.vcf | tail     #to visualize the last ten lines 
bcftools view -h var.vcf | tail -1  #last line of the header. This last line is particularly important as it shows what each field in the main body of the vcf is and it also gives the individual names. Do you remember what information is given in the different fields?
bcftools view -H var.vcf   #to output every raw call
bcftools view -H var.vcf | head     #to visualize the first ten lines  
bcftools view -H var.vcf | head -1   # just the first site present in the file
```
We can index our vcf to make it more accessible. But we need to compress it first
```
bgzip var.vcf
bcftools index var.vcf.gz
```

Let's dive more into the vcf fields 
```
bcftools view var.vcf.gz | grep -m 1 -A 1 "#CHROM" | cut -f 1-7     #to extract the first seven columns: what information do we have here?
bcftools view -H var.vcf.gz | head -1 | cut -f 8 #now let's extract the info field for the first site in our vcf
bcftools view -h var.vcf.gz | grep "INFO" # What do these parameters mean? The header of the file defines every field present in the vcf. So let's grep the information from the header.
bcftools view -H var.vcf.gz | head -1 | cut -f 9 #Now let's explore the genotype information (present in the format field): what are the fields present in this section?
bcftools view -h var.vcf.gz | grep "##FORMAT" # What do they refer to? Let's grep the information from the header like before.
bcftools view -H var.vcf.gz | head -1 | cut -f 10 #Let's take a look at the first genotype call as an example:
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' var.vcf.gz | head #For another more synthetic way to look at the genotypes in our vcf, we can use the bcftools query utility, a very powerful and a useful tool:
bcftools view -H var.vcf.gz | wc -l   #How many variants are there? This command counts the lines referring to the genetic variants, excluding the header (remember the vcf has one variant per line, excluding the header)
```

Now let's try to calculate some vcf file statistics with vcftools
```
vcftools --gzvcf var.vcf.gz --depth #Calculate mean depth of coverage per individual (in this case we only have 1)
vcftools --gzvcf var.vcf.gz --site-quality #Calculate quality for each site
vcftools --gzvcf var.vcf.gz --minQ 30 --minDP 10 --maxDP 50 --recode --out var_filtered #Let's try to filter our variants for quality and min and max individual depth 
bcftools view -H var_filtered.recode.vcf | wc -l #How many variants do we have left after filtering? Have we somehow filtered our dataset?
```

# Optional: PCA using a vcf from different E. coli samples 
The idea is that everyone sends us their vcf files, so we can merge them in a population vcf file and then we can show you how a PCA can be performed
```
plink --vcf joint.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out ecoli_pop
```
For plotting, we'll have to switch to R studio... (we'll show you the results after the practical)

# HiFi dataset (alignment)
Now let's download the e coli Hifi data (from the google drive link, in the appropriate subfolder) in our gitpod workspace.  
Align Hifi reads to our reference (e_coli.fasta)
```
minimap2 -ax map-hifi e_coli.fasta SRR11434954.sample.fastq.gz > hifi_alignment.sam 
```
Convert to bam
```
samtools view -b hifi_alignment.sam > hifi_alignment.bam 
```
Sort and index the alignment
```
samtools sort hifi_alignment.bam -o sorted_hifi_alignment.bam -@ 4
samtools index sorted_hifi_alignment.bam
rm hifi_alignment.sam # we don't need the sam anymore
```
Now try to load this HiFi alignment on IGV as you did before and compare it with the alignment previously made with the Illumina dataset. What are the main differences? 


# HiFi dataset (genome assembly)
Now we will try to make a de novo assembly using our HiFi dataset. First, we will have to downsample our dataset to make the operation computationally feasible for use. We will use rasusa, a tool specifically designed to randomly subsample long-reads data.  
Subsample to 20x coverage
```
rasusa --coverage 20 --genome-size 4641652b --input SRR11434954.sample.fastq.gz -o 20x_SRR11434954.sample.fastq.gz 
```
Perform de novo assembly using our downsampled dataset (it will take few minutes)
```
hifiasm -o SRR11434954.asm -t 4 20x_SRR11434954.sample.fastq.gz
```

Extract FASTA sequence from GFA files (covert from gfa to fasta)
```
awk '/^S/{print ">"$2"\n"$3}' SRR11434954.asm.bp.hap1.p_ctg.gfa | fold > SRR11434954.asm.bp.hap1.p_ctg.fasta
```

Perform assembly evaluation using gfastats; you can also run it of the gfa, however the results will be different (why?)
```
gfastats SRR11434954.asm.bp.hap1.p_ctg.fasta
```

The file SRR11434954.asm.bp.hap1.p_ctg.gfa (not the fasta) can be visualized using the tool Bandage. You can download it from here https://rrwick.github.io/Bandage/ (choose a version appropriate for your OS).
Download the tool in your computer, download locally your Hifiasm assembly and open it in Bandage (remember to press "draw graph" after loading the gfa file).


# Optional: Blast alignment
We will create our own blast database using yeast reference sequence; download reference sequence from here (as you did for E coli) https://www.ncbi.nlm.nih.gov/genome/?term=saccaromyces+cerevisiae  
Rename and extract reference sequence:
```
mv GCF_000146045.2_R64_genomic.fna.gz yeast_genomic.fna.gz 
gunzip yeast_genomic.fna.gz
```
Make blast database
```
makeblastdb -in yeast_genomic.fna -parse_seqids -dbtype nucl 
```
Try to download Alk1 yeast gene sequence from NCBI (by looking within the gene database) and send it to a fasta file "alk1.fasta"  
Perform blast alignment:
```
blastn -outfmt 6 -query alk1.fasta -db yeast_genomic.fna
```
Take a look at the different columns of the blast output:

1.  qseqid      query or source (e.g., gene) sequence id

2.  sseqid      subject  or target (e.g., reference genome) sequence id

3.  pident      percentage of identical matches

4.  length      alignment length (sequence overlap)

5.  mismatch    number of mismatches

6.  gapopen     number of gap openings

7.  qstart      start of alignment in query

8.  qend        end of alignment in query

9.  sstart      start of alignment in subject

10.  send        end of alignment in subject

11.  evalue      expect value (the number of expected hits of similar quality (score) that could be found just by chance)

12.  bitscore    bit score (the required size of a sequence database in which the current match could be found just by chance. The bit-score is a log2 scaled and normalized raw-score)

Here we used as query a gene taken from the yeast annotated genome (same species) as an example. You can try building databases of another species and then "blasting" gene sequences from those species to try and detect homologous genes if you have some time left.

