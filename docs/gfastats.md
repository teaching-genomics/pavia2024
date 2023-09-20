Let's get some test files first:
```
mv testFiles-gfastats/* .
```
Help:
```
gfastats -h
```
File:
```
cat random1.fasta
```
Summary statistics:
```
gfastats random1.fasta
```
Tabular output:
```
gfastats random1.fasta -t
```
Change locale:
```
gfastats large_input.fasta.gz --locale en_US.UTF-8
```
Full output:
```
gfastats random1.fasta --nstar-report
```
Report by sequence:
```
gfastats random1.fasta --seq-report
```
Original file:
```
gfastats random1.fasta -ofa
```
Line length:
```
gfastats random1.fasta -ofa --line-length 2
```
Subset:
```
gfastats random1.fasta Header2 -ofa
```
Subset with bed:
```
gfastats random1.fasta -e <(echo Header2) -ofa
```
cat random1.fasta.bed`
```
gfastats random1.fasta -ofa -e random1.fasta.bed
```
```
gfastats random1.fasta -ofa -i random1.fasta.bed
```
Size of components:
```
gfastats random1.fasta -s s
```
```
gfastats random1.fasta -s c
```
```
gfastats random1.fasta -s g
```
AGP:
```
gfastats random1.fasta -b a
```
BED coordinates:
```
gfastats random1.fasta -b s
```
```
gfastats random1.fasta -b c
```
```
gfastats random1.fasta -b g
```
Sorting:
```
gfastats random1.fasta -ofa --sort largest
```
```
gfastats random1.fasta -ofa --sort descending
```
```
gfastats random1.fasta -ofa --sort test.sort
```
GFA2:
```
gfastats random1.gfa2 -o gfa2
```
GFA2 to FASTA conversion:
```
gfastats random1.gfa2 -o fasta
```
GFA2 to GFA1 conversion:
```
gfastats random1.gfa2 -o gfa
```
GFA1:
```
gfastats random2.gfa -o gfa
```
GFA1 to FASTA:
```
gfastats random2.gfa -o fasta
```
GFA1 to GFA2:
```
gfastats random2.gfa -o gfa2
```
GFA1 no sequence:
```
gfastats random2.noseq.gfa -o gfa
```
GFA1 no sequence:
```
gfastats random2.noseq.gfa -o fa
```
Homopolymer compression:
```
gfastats random1.fasta --homopolymer-compress 1 -ofa
```
Find terminal overlaps:
```
gfastats random5.findovl.gfa -ogfa
```
```
gfastats random5.findovl.gfa --discover-terminal-overlaps 3 -ogfa
```
Discover paths:
```
gfastats random1.fasta -ogfa | grep -v "^P" > test.gfa
```
```
gfastats test.gfa -ogfa
```
```
gfastats test.gfa -ogfa2 --discover-paths
```
Superimpose AGP:
```
gfastats random1.fasta -a random1.agp -ofa
```
SAK reverse complement:
```
cat random1.rvcp.sak
```
```
gfastats random1.fasta -ofa
```
```
gfastats random1.fasta -k random1.rvcp.sak -ofa
```
Other SAK instructions:
```
cat random1.instructions.sak
gfastats random1.fasta -ofa
gfastats random1.fasta -ofa -k <(head -1 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -2 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -3 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -4 random1.instructions.sak)
gfastats random1.fasta -ogfa2 -k <(head -4 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -5 random1.instructions.sak)
gfastats random1.fasta -ogfa2 -k <(head -5 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -6 random1.instructions.sak)
gfastats random1.fasta -ogfa2 -k <(head -6 random1.instructions.sak)
gfastats random1.fasta -ogfa2 -k <(head -6 random1.instructions.sak)
gfastats random1.fasta -ogfa2 -k <(head -7 random1.instructions.sak)
gfastats random1.fasta -ofa -k <(head -8 random1.instructions.sak)
```
