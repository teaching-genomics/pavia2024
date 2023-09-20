```
kreeq validate -f input.[fasta|fastq][.gz] -r reads1.fastq[.gz] reads2.fastq[.gz] [...] [-k 21]
```

It accepts multiple read files as input, separated by space. To check out all options and flags use `kreeq -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
kreeq validate -f testFiles/random1.fasta -r testFiles/random1.fastq
```

Importantly, the kreeq database can only be computed once on the read set, and reused for multiple analyses to save runtime:

```
kreeq validate -r testFiles/random1.fastq -o db.kreeq
kreeq validate -f testFiles/random1.fasta -d db.kreeq
```

Similarly, kreeq databases can be generated separately for multiple inputs and combined, with increased performance in HPC environments:

```
kreeq validate -r testFiles/random1.fastq -o random1.kreeq
kreeq validate -r testFiles/random2.fastq -o random2.kreeq

kreeq union -d random1.kreeq random2.kreeq -o union.kreeq
kreeq validate -f testFiles/random1.fasta -d union.kreeq
```