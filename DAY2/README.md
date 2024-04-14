# BGPAV 2024 Day2

## Tools used 
For this section we will use the tools vg (https://github.com/vgteam/vg), ODGI (https://github.com/pangenome/odgi), Bandage, PGGB (https://github.com/pangenome/pggb#installation)  
As for DAY1, all tools will be already installed in your workspace.

## Pangenome graph construction (vg and odgi) and visualization
We will build our first graphs using vg construct.  
Run it without parameters to get information on its usage  
```
vg construct
```
Let’s start by constructing a graph from one single sequences. This will create a graph that just consists of a linear chain of nodes, each with 32 characters.
```
vg construct -r data/tiny.fa -m 32 > tiny.ref.vg
```
The -m flag tells vg to put at most 32 characters into each graph node (you might want to run it with different values and observe the different results).  
To visualize a graph, you can use vg view. As default, vg view will output a graph in GFA format.  
```
vg view tiny.ref.vg
```
Now try to vary the parameter passed to -m of vg construct and visualize the result.  
Now let’s build a new graph with some variants built into it. First, take a look at at data/tiny.vcf, which contains 2 SNV variants in VCF format.
```
vg construct -r data/tiny.fa -v data/tiny.vcf.gz -m 1024 > tiny.vg  ###we set -m at his maximum value (1024) to avoid graph splitting
```
Visualize the resulting graph in Bandage.  

## Saccharomyces cerevisiae pangenome graphs construction
In your workspace you will find  7 Saccharomyces cerevisiae assemblies (Yue et al., 2016) from the Yeast Population Reference Panel (YPRP) plus the SGD reference.  
Before constructing the pangenome graph, we need to rename the sequences according to the PanSN-spec specification, using ```fastix```. The sequence names have to follow such a scheme:  
```
DBVPG6044#1#chrI
DBVPG6044#1#chrII
...
S288C#1#chrI
S288C#1#chrII
...
SGDref#1#chrI
SGDref#1#chrII
...
```
To add the right prefix to all assemblies and concatenate all sequences in a single FASTA, you can use this command:
```
for file in *gz; do fastix -p "${file%.fa.gz}#1#" <(zcat $file) >> all.fa; done
```
Index the fasta file:  
```
samtools faidx all.fa
```
Let's create a folder to store our graphs:
```
mkdir graphs
```




pggb exposes parameters that allow users to influence the structure of the graph that will represents the input sequences. In particular, reducing the mapping identity (-p parameter) increases the sensitivity of the alignment, leading to more compressed graphs. It is recommended to change this parameter depending on how divergent are the input sequences.  
Assuming we will work with one chromosome at a time, we estimate the sequence divergence for each set of chromosomes. To partition the sequences by chromosome, execute:  
```
cut -f 1 scerevisiae8.fasta.gz.fai | cut -f 3 -d '#' | sort | uniq | while read CHROM; do
    CHR_FASTA=scerevisiae8.$CHROM.fasta.gz
    samtools faidx scerevisiae8.fasta.gz $(grep -P $"$CHROM\t" scerevisiae8.fasta.gz.fai | cut -f 1) | bgzip -@ 8 > $CHR_FASTA

    echo "Generated $CHR_FASTA"
done
```
For each set of chromosomes, to estimate the distance of each input sequence to every other sequence in the set, we use mash. In particular, the mash triangle command outputs a lower-triangular distance matrix. For example, to compute the distances between all mitochondrial sequences, execute:  
```
mash triangle scerevisiae8.chrMT.fasta.gz -p 8 > scerevisiae8.chrMT.mash_triangle.txt
cat scerevisiae8.chrMT.mash_triangle.txt | column -t
```
The distance is between 0 (identical sequences) and 1. Which is the maximum divergence? Find the maximum divergence for all chromosomes.  
From this analysis, chrVII, chrVIII, and chrXIII sets should show the higher sequence divergence, with maximum value of ~ 0.0639874. In general, we should set a mapping identity value lower than or equal to 100 - max_divergence * 100. That is, to analyze such a dataset, we have to specify -p lower than or equal to 93.60126. However, in order to account for possible underestimates of sequence divergence, and medium/large structural variants leading locally to greater divergence, we can set an even smaller mapping identity, like -p 90.    
Run pggb on each chromosome separately (including the mitochondrial one, chrMT), calling variants with respect to the SGDref assembly. For example, for chromosome 1, run:  
```
cd ~/day3_yeast
mkdir -p graphs
pggb -i assemblies/scerevisiae8.chrI.fasta.gz -s 20000 -p 90 -n 8 -t 8 -G 7919,8069 -o graphs/scerevisiae8.chrI -V SGDref:#
```
The -V SGDref:# parameter specify to call variants with respect to the path having as sample name SGDref; the sample name is identified by using # as separator (take a look at the PanSN-spec specification).  
Count the number of variants in the VCF file for each chromosomes. Which chromosome has the least named variants? Are there chromosomes on which 0 variants are called? If so, why? Hint: take a look at the alignments in the PAF file produced by pggb (take a look at the PAF format specification too). Try to understand what is happing and how the problem can be solved.  
Use odgi squeeze to put all the graphs together in the same file. Try to visualize the graph layout of the squeezed graph (the graph just generated with all chromosomes together) with odgi layout and odgi draw.  
Run pggb on all chromosomes jointly, giving in input the scerevisiae8.fasta.gz file. Call variants too with respect to the SGDref sample. Take a look at the graph layout. Is the layout of the graph obtained by running pggb separately on each chromosome (and then combining its results) similar to the graph obtained by running pggb on all chromosomes jointly?  
In the newly obtained VCF, count how many variants are called for each chromosome. Are these counts similar to those obtained by calling variants on each chromosome separately? Are there variants esclusive (that do not appear in the chromosome-based VCF files) to this new VCF file?  

## Sequence partitioning
Move to the directory created in the previous activity on Saccharomyces cerevisiae:  
```
cd ~/day3_yeast
```
We can’t really expect to pairwise map all sequences together and obtain well separated connected components. It is likely to get a giant connected component, and probably a few smaller ones, due to incorrect mappings or false homologies. This might unnecessarily increase the computational burden, as well as complicate the downstream analyzes. Therefore, it is recommended to split up the input sequences into communities in order to find the latent structure of their mutual relationship.  
We need to obtain the mutual relationship between the input assemblies in order to detect the underlying communities. To compute the pairwise mappings with wfmash, execute:  
```
cd assemblies
wfmash scerevisiae8.fasta.gz -p 90 -n 7 -t 8 -m > scerevisiae8.mapping.paf
```
Why did we specify -n 7?  
To project the PAF mappings into a network format (an edge list), execute:  
```
python3 ~/pggb/scripts/paf2net.py -p scerevisiae8.mapping.paf
```
The paf2net.py script creates 3 files:  
```
scerevisiae8.mapping.paf.edges.list.txt is the edge list representing the pairs of sequences mapped in the PAF;
scerevisiae8.mapping.paf.edges.weights.txt is a list of edge weights (long and high estimated identity mappings have greater weight);
scerevisiae8.mapping.paf.vertices.id2name.txt is the ‘id to sequence name’ map.
```
To identify the communities, execute:
```
python3 ~/git/pggb/scripts/net2communities.py \
    -e scerevisiae8.mapping.paf.edges.list.txt \
    -w scerevisiae8.mapping.paf.edges.weights.txt \
    -n scerevisiae8.mapping.paf.vertices.id2name.txt
```
How many communities were detected?  
The paf2net.py script creates a set of *.community.*.txt files one for each the communities detected. Each txt file lists the sequences that belong to the same community.  
Are there communities that contain multiple chromosomes? Which ones?  
Identity the communities again, but this time add the --plot option to visualize them too:  
```
python3 ~/git/pggb/scripts/net2communities.py \
    -e scerevisiae8.mapping.paf.edges.list.txt \
    -w scerevisiae8.mapping.paf.edges.weights.txt \
    -n scerevisiae8.mapping.paf.vertices.id2name.txt \
    --plot
```
Take a look at the scerevisiae8.mapping.paf.edges.list.txt.communities.pdf file.  
Write a little script that take the *.community.*.txt files in input and create the corresponding FASTA files, ready to be input to pggb. Run pggb on the communities with multiple chromosomes and compare the results (layout and variants) from the previous activities.  









A set of pairwise alignments implies a variation graph, so pangenome graphs can be obtained from alignments too. Use minimap2 and seqwish to build graphs from the HLA gene haplotypes  
```
ln -s ~/vg/test/GRCh38_alts/FASTA/HLA/
minimap2 -c -x asm20 -X -t 8 HLA/DRB1-3123.fa HLA/DRB1-3123.fa > DRB1-3123.paf
seqwish -s HLA/DRB1-3123.fa -p DRB1-3123.paf -g DRB1-3123.gfa
odgi build -g DRB1-3123.gfa -o - | odgi sort -i - -o DRB1-3123.og
odgi viz -i DRB1-3123.og -o DRB1-3123.png -x 2000
```
Use Bandage to visualize these graphs you just built.  

## Pangenome graph construction using pggb
In this section we will build HLA pangenome graphs using pggb  
```
mkdir day1_pggb
cd day1_pggb
ln -s ~/pggb/data
```
The human leukocyte antigen (HLA) system is a complex of genes on chromosome 6 in humans which encode cell-surface proteins responsible for the regulation of the immune system.  
Let’s build a pangenome graph from a collection of sequences of the DRB1-3123 gene:  
```
pggb -i data/HLA/DRB1-3123.fa.gz -n 12 -t 8 -o out_DRB1_3123
```
Run pggb without parameters to get information on the meaning of each parameter:  
```
pggb
```
Take a look at the files in the out_DRB1_3123 folder. Visualize the graph with Bandage.  
Why did we specify -n 12?  
How many alignments were executed during the pairwise alignment (take a look at the PAF output)? Visualize the alignments:  
```
cd out_DRB1_3123
~/wfmash/scripts/paf2dotplot png large *paf
cd ..
```
Use odgi stats to obtain the graph length, and the number of nodes, edges, and paths. Do you think the resulting pangenome graph represents the input sequences well? Check the length and the number of the input sequences to answer this question.  
How many blocks were selected and ‘smoothed’ during the two rounds of graph normalization (take a look at the *.log file to answer this question)?  
Try building the same pangenome graph by specifying a lower percent identity (-p 95 by default):  
```
pggb -i data/HLA/DRB1-3123.fa.gz -p 90 -n 12 -t 8 -o out2_DRB1_3123
```
Check graph statistics. Does this pangenome graph represent better or worse the input sequences than the previously produced graph?  
Try to decrease the number of mappings to reteain for each segment:  
```
pggb -i data/HLA/DRB1-3123.fa.gz -p 90 -n 6 -t 8 -o out3_DRB1_3123
```
How does it affect the graph?  
Try to increase the target sequence length for the partial order alignment (POA) problem (-G 4001,4507 by default):  
```
pggb -i data/HLA/DRB1-3123.fa.gz -p 90 -n 12 -t 8 -G 12000,13000 -o out4_DRB1_3123
```
How is this changing the runtime and the memory usage? How is this affecting graph statistics? How many blocks were selected and ‘smoothed’ during the two rounds of graph normalization?  
Try 1, 3 or 4 rounds of normalization (for example,by specifying -G 4001, -G 4001,4507,4547, or -G 4001,4507,4547, 4999). How does this affect graph statistics?  
Take the second pggb run and try to increase the segment length (-s 10000 by default):  
```
pggb -i data/HLA/DRB1-3123.fa.gz -s 20000 -p 90 -n 12 -t 8 -o out5_DRB1_3123
```
How is this affecting graph statistics? Why?  
pggb produces intermediate graphs during the process. Let’s keep all of them:  
```
pggb -i data/HLA/DRB1-3123.fa.gz -p 90 -n 12 -t 8 --keep-temp-files -o out2_DRB1_3123_keep_intermediate_graphs
```
What does the file with name ending with .seqwish.gfa contain? and what about the file with name ending with .smooth.1.gfa? Take a look at the graph statistics of all the GFA files in the out2_DRB1_3123_keep_intermediate_graphs folder.  
Choose another HLA gene from the data folder and explore how the statistics of the resulting graph change as s, p, n change. Produce scatter plots where on the x-axis there are the tested values of one of the pggb parameters (s, p, or n) and on the y-axis one of the graph statistics (length, number of nodes, or number of edges). You can do that using the final graph and/or the intermediate ones.  
Now we will try to build LPA pangenome graphs.  
Lipoprotein(a) (LPA) is a low-density lipoprotein variant containing a protein called apolipoprotein(a). Genetic and epidemiological studies have identified lipoprotein(a) as a risk factor for atherosclerosis and related diseases, such as coronary heart disease and stroke.  
Try to make LPA pangenome graphs. The input sequences are in data/LPA/LPA.fa.gz. Sequences in this locus have a peculiarity: which one? Hint: visualize the alignments and take a look at the graph layout (with Bandage and/or in the .draw_multiqc.png files).  


## Read mapping and variant calling
```
ls ~/vg/test/1mb1kgp
```
This directory contains 1Mbp of 1000 Genomes data for chr20:1000000-2000000. As for the tiny example, let’s’ build a graph that encodes the known sequence variation.
```
vg construct -r 1mb1kgp/z.fa -v 1mb1kgp/z.vcf.gz -m 1024 > z.vg
```
Let’s build the indexes needed for the mapping:
```
vg index -L -x z.xg z.vg
vg index -g z.gcsa -k 16 z.vg
```
Let's download the reads to map:
```
wget -c http://hypervolu.me/~guarracino/CPANG22/NA12878.20_1M-2M.30x.bam
samtools fastq -1 NA24385-20_1M-2M-30x_1.fq.gz -2 NA24385-20_1M-2M-30x_2.fq.gz NA12878.20_1M-2M.30x.bam 
```
The sample was used in the preparation of the 1000 Genomes results.
Map reads against the built graph:
```
vg map -d z -f NA24385-20_1M-2M-30x_1.fq.gz -f NA24385-20_1M-2M-30x_2.fq.gz > NA24385-20_1M-2M-30x.gam
```
For variant calling, compute the read support first:
```
vg pack -x z.xg -g NA24385-20_1M-2M-30x.gam -o NA24385-20_1M-2M-30x.pack -Q 5
```
Compute the snarls (using fewer threads with -t can reduce memory at the cost of increased runtime):
```
vg snarls z.xg > z.snarls
```
Call variants:
```
vg call z.xg -r z.snarls -k NA24385-20_1M-2M-30x.pack > NA24385-20_1M-2M-30x.vcf 
```

