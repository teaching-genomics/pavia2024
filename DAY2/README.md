# BGPAV 2024 Day2

## Tools used 
For this section we will use the tools vg (https://github.com/vgteam/vg), Bandage, PGGB (https://github.com/pangenome/pggb#installation)  
As for DAY1, all tools will be already installed in your workspace (except bandage, you'll have to install it locally).
For DAY2, you will have to delete the old workspace and create a new one.

## Simple examples of pangenome graph construction (vg and odgi) and visualization
We will build our first graphs using vg construct.  
First run it without parameters to get information on its usage  
```
vg construct
```
Move to /pavia2024/DAY2 folder (using cd command). Let’s start by constructing a graph from one single sequence. This will create a graph that just consists of a linear chain of nodes, each with 32 characters.
```
vg construct -r data/tiny/tiny.fa -m 32 > tiny.ref.vg
```
The -m flag tells vg to put at most 32 characters into each graph node (you might want to run it with different values and observe the different results).  
To visualize a graph, you can use ```vg view```. As default, vg view will output a graph in GFA format, that you can then visualize in Bandage.
```
vg view tiny.ref.vg  ### inspect the main fields of the gfa format
vg view tiny.ref.vg > tiny.ref.gfa   ###write it to an output file
```
Now let’s build a new graph with some variants built into it. First, take a look at at data/tiny.vcf, which contains 2 SNV variants in VCF format.
```
vg construct -r data/tiny/tiny.fa -v data/tiny/tiny.vcf -m 1024 > tiny_vcf.vg   ###we set -m at his maximum value (1024) to avoid graph splitting
vg view tiny_vcf.vg  > tiny_vcf.gfa   ###write it to an output file
```
Download the gfa locally and visualize the resulting graph in Bandage.  


## Saccharomyces cerevisiae pangenome graphs construction
In your workspace you will find 7 Saccharomyces cerevisiae assemblies (Yue et al., 2016) from the Yeast Population Reference Panel (YPRP) plus the SGD reference.  
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
cd data/yeast/
for file in *gz; do fastix -p "${file%.fa.gz}#1#" <(zcat $file) >> all.fa; done
```
Index the fasta file:  
```
samtools faidx all.fa
```
Let's extract mitogenomes. We will then generate some statistics and generate the fasta index:  
```
grep MT all.fa.fai | cut -f1 > mt.ls
gfastats -i mt.ls all.fa -o all.mt.fa
samtools faidx all.mt.fa
```
We will generate pangenome graphs from the mitoassemblies we just extracted. Let's create a folder to store our graphs:
```
cd ../..
mkdir -p graphs
```
Let's generate pangenome graphs from mitoassemblies using pggb
```
pggb -i data/yeast/all.mt.fa -p 90 -n 8 -t 4 -o graphs/scerevisiae8.mt 
```

## Read mapping and variant calling in Saccharomyces cerevisiae pangenome graphs
Now let's use the pangenome graph we just built to map short reads from yeast and perform variant calling.  
We will simulate a bunch of reads from one of the mitoassemblies used to build the graph and use these to perform the alignment.
```
cd data/yeast
gfastats Y12.fa.gz chrMT -o Y12-mt.fa      ### extract mitogenome from Y12 strain
wgsim -N 10000 Y12-mt.fa Y12-mt.read1.fq Y12-mt.read2.fq 
```
Now let's align these reads to the pangenome graph we just built.  

```
cd ../../graphs/scerevisiae8.mt
```
We will use the "all.mt.fa.bf3285f.11fba48.9977ae5.smooth.fix.gfa" graph. We will rename it for brevity and generate the proper indexes for vg map (the general-use mapper in the vg toolkit)
```
mv all.mt.fa.bf3285f.11fba48.9977ae5.smooth.fix.gfa all.mt.smooth.fix.gfa
vg mod -X 256 all.mt.smooth.fix.gfa > mod_all.mt.smooth.fix.gfa  ####chop nodes in the graph so they are not greater than 1024 bases long
vg convert -g mod_all.mt.smooth.fix.gfa --ref-sample SK1 -x > mod_all.mt.smooth.fix.xg  ####convert to xg format and assign one haplotype as reference (it will change haplotypes for this sample to reference paths)
vg index -g mod_all.mt.smooth.fix.gcsa -k 16 mod_all.mt.smooth.fix.xg ###generate index needed by vg map
```

Now let's map the simulated reads we generated
```
vg map -d mod_all.mt.smooth.fix -f ../../data/yeast/Y12-mt.read1.fq -f ../../data/yeast/Y12-mt.read2.fq -t 4 > alignment.gam
```

For variant calling, let's compute the read support first:
```
vg pack -x mod_all.mt.smooth.fix.xg -g alignment.gam -o all.mito.pack -Q 5
```
Compute the snarls:
```
vg snarls mod_all.mt.smooth.fix.xg > all.mito.snarls
```
Call the variants:
```
vg call mod_all.mt.smooth.fix.xg -r all.mito.snarls -k all.mito.pack -d 1 > mito.vcf 
```
Let's inspect the vcf to see which variants have been identified.