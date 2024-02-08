## Assembly of repeats from genome-skimming data

> genome-skimming (here): short or long read sequencing with an approximate 1× genome coverage



### Aim

Construct repeat assemblies from genome-skimming data to answer repeat-centered questions such as linkage or separation of 45S/35S-rDNA and 5S-rDNA.



### Preliminary work

No preliminary work necessary other than sequencing.



### Dependencies

- MEGAHIT (tested with 1.2.9)
- SPAdes (tested with 3.15.5)
- seqtk (tested with 1.3-r117 and 1.3-r106)
- seqkit (tested with 2.4.0 and 2.3.0)
- SRAtools (optional, recommended for bulk downloading of data)
- fastQC or similar (optional, recommended for quality control)
- trimmomatic, cutadapt or similar (optional, recommended for quality trimming)



### Data collection for repeat assembly

Depending on your specific question a genome coverage from 0.05× up to 1× might be sufficient. From our experience the detection of linkage or separation of rDNAs worked well in all tested ranges from 0.3× - 1× (lower not tested).

It is recommended to inspect your data using quality control and trimming tools (for suggested tools see Dependencies).



#### Genome-skimming assembly from short reads

Assemble repeats from short read data with a low coverage (~ 1×).

**Note**: The repeat assembly does not differ from a genome assembly using the stated tools. However, by using a low coverage it can be expected that only repeats (and organelle DNA) will be assembled in a meaningful manner. Therefore, it is recommended to not use more read data than a 1× coverage. Use tools like seqkit or seqtk to sample your read data to an appropriate size.

For the first round of the assembly MEGAHIT is used (or further details see ``megahit -h``, depending on your resources adjust the cpu and memory limits with ``-t`` and ``-m`` options):

```bash
mkdir megahit_out
megahit -1 <reads.for.fastq> -2 <reads.rev.fastq> -o meaghit_out --presets meta-large -t 6 --min-contig-len 5000
```

For the second round of the assembly SPAdes is used (or further details see ``spades -h``, depending on your resources adjust the cpu and memory limits with ``-t`` and ``-m`` options):

```bash
spades -1 <reads.for.fastq> -2 <reads.rev.fastq> -t 6 --trusted-contigs megahit_out/final.contigs.fa -o spades_out --isolate --cov-cutoff 20
```



#### Genome-skimming assembly from long reads

Assemble repeats from short read data with a low coverage (~ 1×).

**Note**: The repeat assembly does not differ from a genome assembly using the stated tools. However, by using a low coverage it can be expected that only repeats (and organelle DNA) will be assembled in a meaningful manner. Therefore, it is recommended to not use more read data than a 1× coverage. Use tools like [``long_read_sampling``](https://toolshed.g2.bx.psu.edu/view/petr-novak/long_reads_sampling/5596bafd2119) (RepeatExplorer2 Galaxy tools) sample your read data to an appropriate size.

For the first round of the assembly MEGAHIT is used (or further details see ``megahit -h``, depending on your resources adjust the cpu and memory limits with ``-t`` and ``-m`` options):

```bash
mkdir megahit_out
megahit -r <long.reads.fastq> -o megahit_out --presets meta-large -t 6 --min-contig-len 5000
```

For the second round of the assembly SPAdes is used (or further details see ``spades -h``, depending on your resources adjust the cpu and memory limits with ``-t`` and ``-m`` options):

```bash
spades -s <long.reads.fastq> -t 6 --trusted-contigs megahit_out/final.contigs.fa -o spades_out --isolate --cov-cutoff 20
```



### Visualization of bandage graphs

Use the bandage graph viewer according to its manual with ``.gfa`` files in the ``spades_out`` folder. You might want to use the built-in blast to search for desired repeats. Use the draw graph around blast hits function to get a less cluttered view of your "annotated" repeats.

