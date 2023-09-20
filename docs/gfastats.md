Help:
```
gfastats -h
```
File:
```
cat testFiles/random1.fasta
```
Summary statistics:
```
gfastats testFiles/random1.fasta
```
Tabular output:
```
gfastats testFiles/random1.fasta -t
```
Change locale:
```
gfastats large_input.fasta.gz --locale en_US.UTF-8
```
Full output:
```
gfastats testFiles/random1.fasta --nstar-report
```
Report by sequence:
```
gfastats testFiles/random1.fasta --seq-report
```
Original file:
```
gfastats testFiles/random1.fasta -ofa
```
Line length:
```
gfastats testFiles/random1.fasta -ofa --line-length 2
```
Subset:
```
gfastats testFiles/random1.fasta Header2 -ofa
```
Subset with bed:
```
gfastats testFiles/random1.fasta -e <(echo Header2) -ofa
```
cat testFiles/random1.fasta.bed`
```
gfastats testFiles/random1.fasta -ofa -e testFiles/random1.fasta.bed
```
```
gfastats testFiles/random1.fasta -ofa -i testFiles/random1.fasta.bed
```
Size of components:
```
gfastats testFiles/random1.fasta -s s
```
```
gfastats testFiles/random1.fasta -s c
```
```
gfastats testFiles/random1.fasta -s g
```
AGP:
```
gfastats testFiles/random1.fasta -b a
```
BED coordinates:
```
gfastats testFiles/random1.fasta -b s
```
```
gfastats testFiles/random1.fasta -b c
```
```
gfastats testFiles/random1.fasta -b g
```
Sorting:
```
gfastats testFiles/random1.fasta -ofa --sort largest
```
```
gfastats testFiles/random1.fasta -ofa --sort descending
```
```
gfastats testFiles/random1.fasta -ofa --sort test.sort
```
GFA2:
```
gfastats testFiles/random1.gfa2 -o gfa2
```
GFA2 to FASTA conversion:
```
gfastats testFiles/random1.gfa2 -o fasta
```
GFA2 to GFA1 conversion:
```
gfastats testFiles/random1.gfa2 -o gfa
```
GFA1:
```
gfastats testFiles/random2.gfa -o gfa
```
GFA1 to FASTA:
```
gfastats testFiles/random2.gfa -o fasta
```
GFA1 to GFA2:
```
gfastats testFiles/random2.gfa -o gfa2
```
GFA1 no sequence:
```
gfastats testFiles/random2.noseq.gfa -o gfa
```
GFA1 no sequence:
```
gfastats testFiles/random2.noseq.gfa -o fa
```
Homopolymer compression:
```
gfastats testFiles/random1.fasta --homopolymer-compress 1 -ofa
```
Find terminal overlaps:
```
gfastats testFiles/random5.findovl.gfa -ogfa
```
```
gfastats testFiles/random5.findovl.gfa --discover-terminal-overlaps 3 -ogfa
```
Discover paths:
```
gfastats testFiles/random1.fasta -ogfa | grep -v "^P" > test.gfa
```
```
gfastats test.gfa -ogfa
```
```
gfastats test.gfa -ogfa2 --discover-paths
```
Superimpose AGP:
```
gfastats testFiles/random1.fasta -a testFiles/random1.agp -ofa
```
SAK reverse complement:
```
cat testFiles/random1.rvcp.sak
```
```
gfastats testFiles/random1.fasta -ofa
```
```
gfastats testFiles/random1.fasta -k testFiles/random1.rvcp.sak -ofa
```
Other SAK instructions:
```
cat testFiles/random1.instructions.sak
gfastats testFiles/random1.fasta -ofa
gfastats testFiles/random1.fasta -ofa -k <(head -1 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -2 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -3 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -4 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ogfa2 -k <(head -4 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -5 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ogfa2 -k <(head -5 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -6 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ogfa2 -k <(head -6 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ogfa2 -k <(head -6 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ogfa2 -k <(head -7 testFiles/random1.instructions.sak)
gfastats testFiles/random1.fasta -ofa -k <(head -8 testFiles/random1.instructions.sak)
```
