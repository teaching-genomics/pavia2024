The standard notation of rdeval is:
``` 
rdeval input.[fasta|fastq|gfa][.gz] [expected genome size]
```

Let's begin by looking at the options of redeval:
```
rdeval -h
```

Now let's evaluate the contents of our fasta file:
```
rdeval random1.fasta
```

And filter sequences with a length greater than 10:
```
rdeval random1.fasta -f ">10"
```
