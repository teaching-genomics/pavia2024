The standard notation for using kreeq is as follows:
```
kreeq validate -f input.[fasta|fastq][.gz] -r reads1.fastq[.gz] reads2.fastq[.gz] [...] [-k 21]
```
It accepts multiple read files as input, separated by space. The two modes we will be using today are `validate` and `union`.
To check out all options and flags use:
```
kreeq -h
kreeq validate -h
kreeq union -h
```  

Let's get some test files first:
```
mv gfastar/docs/testFiles-kreeq/* .
```

We will test some typical usage with the files moved from the `testFiles` folder, e.g.:
```
kreeq validate -f random1.fasta -r random1.fastq
kreeq validate -f random2.fasta -r random1.fastq random2.fastq
```

Importantly, the kreeq database can only be computed once on the read set, and reused for multiple analyses to save runtime:

```
kreeq validate -r random1.fastq -o random1.kreeq
kreeq validate -f random1.fasta -d random_fa.kreeq
```

Similarly, kreeq databases can be generated separately for multiple inputs and combined, with increased performance in HPC environments:

```
kreeq validate -r random1.fastq -o random1.kreeq
kreeq validate -r random2.fastq -o random2.kreeq

kreeq union -d random1.kreeq random2.kreeq -o union.kreeq
kreeq validate -f random1.fasta -d union.kreeq
```

Now working with real sequencing data:

Let's start by running `gfastats` to get a sense of what we are evaluating:
```
gfastats input.fa
```

Now we are ready to run kreeq:
```
kreeq validate -r filtered.fastq -o filtered.kreeq
kreeq validate -r filtered2.fastq -o filtered2.kreeq

kreeq union -d filtered.kreeq filtered2.kreeq -o filtered_union.kreeq

kreeq validate -f input.fa -d filtered_union.kreeq
```
