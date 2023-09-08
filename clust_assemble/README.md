## Repeat assembly from RE2 super-clusters

### Aim

Construct comprehensive consensus sequences from RepeatExplorer2 (RE2) repeat superclusters (SCL) (or clusters (CL)).



### Preliminary work

In advance to the following work-flow the RepeatExplorer2 pipeline needs to be run including data preparation and quality control. For details on how to run see the [RepeatExplorer2 protocol](https://rdcu.be/b80Gr). You need the output from either a local or galaxy RE2 run.



### Dependencies

- python3 (> 3.6, tested with 3.9.16 and 3.8.8)
- seqtk (tested with 1.3-r117 and 1.3-r106)
- seqkit (tested with 2.4.0 and 2.3.0)
- bowtie2 (tested with 2.4.1 and 2.2.5)
- samtools (tested with 1.17 and 1.12)
- MEGAHIT (tested with 1.2.9)
- SPAdes (tested with 3.15.5)
- DANTE (tested with Galaxy 1.1.0 and Anaconda 0.1.8)
- GNU parallel (optional but recommended, tested with 20230522)



### Data collection for repeat assembly

The assembly uses clustered read data as input. There are three ways to obtain appropriate read data (choose one):

1. Use the reads from superclusters in fasta format (quick)
2. Use the reads from superclusters in fastq format (recommended, faster with RE2 --keep-names)
3. Use the super-cluster contigs to retrieve reads in fastq format (extensive)

For each option the read collection could also be done using individual clusters. However, it is recommended to use superclusters. If you want to use individual clusters instead customize the documentation and scripts accordingly. 

**Important**: make sure that the manual cluster annotation was done properly and super-cluster information was customized if necessary (for details on the curation process see the [RepeatExplorer2 protocol](https://rdcu.be/b80Gr) Procedure 1&2 "Manual correction of automated repeat annotation", and the [video tutorial](https://www.youtube.com/watch?v=gvpJC6gC9Ck)).



#### Recommendations: What to choose as input?

Disclaimer: This is only personal experience and may vary for each individual data! The results may also differ for individual SCL, depending on the underlying repeat family and repeat diversity in the respective genome.

The overall statistically best results may be achieved using super-cluster reads in fastq format (2) if the --keep-names option was enabled during the RE2 run. If the --keep-names option was disabled using super-cluster reads in fasta format (1) will deliver faster and robust results. When interested in smaller super-clusters or in the repeat diversity using reads based on super-cluster contigs (3) is recommended.



#### 1. SCL reads (w/o quality information)

Use the reads from superclusters in fasta format (quick) for the repeat assembly.

Copy the `collect_RE2-SCL_reads.py` in to the RE2 output folder (same folder that also contains `seqclust` directory and `CLUSTER_TABLE.csv`).

```bash
cp /.../collect_RE2-SCL_reads.py /location/of/RE2/output/
```

Run the `collect_RE2-SCL_reads.py` script. By default, it uses the `CLUSTER_TABLE_with_final_annotation.csv` file. If your file is named differently change the input file name in the py-script. If storage space is limited on your machine use the `collect_RE2-SCL_reads-lessstorage.py` script instead (it will only create the read lists).

``` bash
python3 collect_RE2-SCL_reads.py
```

Use the read lists to extract the reads from `seqclust/reads/reads.fasta` 

```bash
for dir in clust_assembly/SCL*/; do
  seqtk subseq seqclust/reads/reads.fasta "${dir}reads.for.lst" > "${dir}reads.for.fasta"
  seqtk subseq seqclust/reads/reads.fasta "${dir}reads.rev.lst" > "${dir}reads.rev.fasta"
done
```

or use parallel to run multithreaded (Note: GNU parallel needs to be installed, adjust the number of cpu cores by modifying -j)

```bash
parallel -j 6 '
  dir={};
  seqtk subseq seqclust/reads/reads.fasta "${dir}reads.for.lst" > "${dir}reads.for.fasta"
  seqtk subseq seqclust/reads/reads.fasta "${dir}reads.rev.lst" > "${dir}reads.rev.fasta"
' ::: clust_assembly/SCL*/;
```



#### 2. SCL reads (w/ quality information)

Use the reads from superclusters in fastq format (recommended) for repeat assembly. 



##### 2.1 w/ --keep-names enabled

**Important**: This is only possible if the --keep-names option was activated while running RE2. If you still want to use fastq reads see section 2.2  or 3. below.

Provide the location of your original read files in fastq format from which the SCL reads can be collected.

```bash
# Directory
export SEQDIR="<path/to/read/files/>"
# File names
export READ_1="<file1.fastq>"
export READ_2="<file2.fastq>"
```

Copy the `collect_RE2-SCL_reads.py` in to the RE2 output folder (same folder that also contains seqclust directory and CLUSTER_TABLE.csv).

```bash
cp /.../collect_RE2-SCL_reads.py /location/of/RE2/output/
```

Run the `collect_RE2-SCL_reads.py` script. By default, it uses the `CLUSTER_TABLE_with_final_annotation.csv` file. If your file is named differently change the input file name in the py-script. If storage space is limited on your machine use the `collect_RE2-SCL_reads-lessstorage.py` script instead (it will only create the read lists).

``` bash
python3 collect_RE2-SCL_reads.py
```

Use the created read lists to extract the reads from your original reads in fastq format.

```bash
for dir in clust_assembly/SCL*/; do
seqtk subseq "${SEQDIR}${READ_1}" "${dir}reads.for.lst" > "${dir}reads.for.fastq"
seqtk subseq "${SEQDIR}${READ_2}" "${dir}reads.rev.lst" > "${dir}reads.rev.fastq"
done
```

or use parallel to run multithreaded (Note: GNU parallel needs to be installed, adjust the number of cpu cores by modifying -j INT, recommended)

```bash
parallel -j 6 '
  dir={};
  seqtk subseq "${SEQDIR}${READ_1}" "${dir}reads.for.lst" > "${dir}reads.for.fastq"
  seqtk subseq "${SEQDIR}${READ_2}" "${dir}reads.rev.lst" > "${dir}reads.rev.fastq"
' ::: clust_assembly/SCL*/;
```



##### 2.2 w/ --keep-names disabled

Following a work-around is shown if the option --keep-names was disabled during the RE2 run. Caution: This will not exactly get the same reads as with the reads list but will result in very similar functionality. Note that this will take substantially longer than extracting reads with a read name list. To get fastq reads the fasta reads are need first. Perform the fasta read collection as described in 1. SCL reads (w/o quality information).

Provide the location of your original read files in fastq format from which the SCL reads can be collected.

```bash
export SEQDIR="<path/to/read/files/>"
export READ_1="<file1.fastq>"
export READ_2="<file2.fastq>"
```

To extract the reads in fastq format use the provided shell script. Copy the script to the RE2 output folder for easier use.

```bash
cp /.../fastq_reads_without_keepnames.sh /location/of/RE2/output/
```

```bash
for dir in clust_assembly/SCL*/; do
  sh fastq_reads_without_keepnames.sh "${dir}" "${SEQDIR}${READ_1}" "${SEQDIR}${READ_2}"
done
```

or use parallel to run multithreaded (Note: GNU parallel needs to be installed, adjust the number of cpu cores by modifying -j INT, recommended)

```bash
parallel -j 6 '
  dir={};
  sh fastq_reads_without_keepnames.sh "${dir}" "${SEQDIR}${READ_1}" "${SEQDIR}${READ_2}"
' ::: clust_assembly/SCL*/;
```



#### 3. SCL reads based on contigs (w/ quality information)

Use the super-cluster contigs to retrieve reads in fastq format (extensive) for repeat assembly.

Cluster reads will be collected by mapping reads in fastq format against the cluster contigs (bowtie2). The reads to be used in the further analysis will be extracted from the mapping (samtools).

Define the location (complete path, and file names) of the read files with quality information to be used. This might be sub-sampled fastq reads or the original (quality trimmed) fastq reads (Note that using the original reads takes more time overall but might harvest better results, by catching a higher diversity of the repeats).

```bash
export SEQDIR="<path/to/read/files/>"
export READ_1="<file1.fastq>"
export READ_2="<file2.fastq>"
```

Copy the `collect_RE2-SCL_contigs.py` in to the RE2 output folder (same folder that also contains seqclust directory and CLUSTER_TABLE.csv).

```bash
cp /.../collect_RE2-SCL_contigs.py /location/of/RE2/output/
```

Run the `collect_RE2-SCL_contigs.py` script. By default, it uses the `CLUSTER_TABLE_with_final_annotation.csv` file. If your file is named differently change the input file name in the py-script. 

``` bash
python3 collect_RE2-SCL_contigs.py
```

Create bowtie index files from the collected SCL contigs.

```bash
for dir in clust_assembly/SCL*/; do
  bowtie2-build "${dir}collected.contigs.fasta" "${dir}collected.contigs.fasta"
done
```

Map the original reads in fastq format against all SCL contigs (adjust the number of cpu cores by modifying -p).

```bash
for dir in clust_assembly/SCL*/; do
  bowtie2 -x "${dir}collected.contigs.fasta" -1 $"${SEQDIR}${READ_1}" -2 $"${SEQDIR}${READ_2}" -S "${dir}scl_contigs_mapped.sam" --local --sensitive-local --no-unal -p 6
done
```

Extract the mapped reads from the alignment files. Run one of the following code blocks depending on if you want to use multiple cpu cores (again, using parallel).

```bash
for dir in clust_assembly/SCL*/; do
  samtools view -b -F 4 "${dir}scl_contigs_mapped.sam" |
  samtools fastq -1 "${dir}contigreads.for.fastq" -2 "${dir}contigreads.rev.fastq" -s "${dir}contigreads_single.fastq" -
 done
```

```bash
parallel -j 6 '
  dir={};
  samtools view -b -F 4 "${dir}scl_contigs_mapped.sam" |
  samtools fastq -1 "${dir}contigreads.for.fastq" -2 "${dir}contigreads.rev.fastq" -s "${dir}contigreads.single.fastq" -
' ::: clust_assembly/SCL*/;
```



### Repeat assembly using MEGAHIT

Assemble the reads for each supercluster using the [MEGAHIT](https://github.com/voutcn/megahit) assembler (for additional options see the MEGAHIT manual, adjust cpu cores used with -t).



#### 1. & 2. SCL reads

Adjust the file extension according to the reads you want to use.

```bash
for dir in clust_assembly/SCL*/; do
  megahit -1 "${dir}reads.for.fasta/q"  -2 "${dir}reads.rev.fasta/q" -o "${dir}mega_assembly" --presets meta-large -t 6
done

```



#### 3. SCL reads based on contigs

```bash
for dir in clust_assembly/SCL*/; do
  megahit -1 "${dir}contigreads.for.fastq"  -2 "${dir}contigreads.rev.fastq" -o "${dir}mega_assembly" --presets meta-large -t 6
done
```



### Repeat assembly using SPAdes

Assemble the reads for each super cluster using the [SPAdes](https://github.com/ablab/spades) assembler (for additional options see the SPAdes manual, adjust cpu cores used with --threads).

To assist the SPAdes assembly collect the SCL contigs.

```bash
python3 collect_RE2-SCL_contigs.py
```

Combine the RE2 SCL contigs with the final MEGAHIT contigs to get the most complete results.

```bash
for dir in clust_assembly/SCL*/; do
  cat "${dir}collected.contigs.fasta" "${dir}mega_assembly/final.contigs.fa" > "${dir}combined.contigs.fa"
done
```

**Note**: When running on a high core count computer the computation time can be reduced using a combination of GNU parallel and the --threads option of SPAdes:

```bash
# example for a 96-core server
time parallel -j 6 '
  dir={};
  spades --threads 16 [...INPUT and OPTIONS, see below]
  ' ::: clust_assembly/SCL*/;
```



#### 1. & 2. SCL reads

Run the assembly for each SCL.

```bash
for dir in clust_assembly/SCL*/; do
  spades --isolate --threads 6 --cov-cutoff auto \
  -1 "${dir}reads.for.fasta/q" \
  -2 "${dir}reads.rev.fasta/q" \
  --trusted-contigs "${dir}combined.contigs.fasta" \
  -o "${dir}assembly"
 done
```



#### 3. SCL reads based on contigs



```bash
for dir in clust_assembly/SCL*/; do
  spades --isolate --threads 6 --cov-cutoff auto \
  -1 "${dir}contigreads.for.fastq" \
  -2 "${dir}contigreads.rev.fastq" \
  -s "${dir}contigreads.single.fastq" \
  --trusted-contigs "${dir}combined.contigs.fasta" \
  -o "${dir}assembly"
 done
```



### Classification of assembled contigs

Classification of SPAdes assembled contigs (NODEs).

#### Classification using DANTE (RexDB)

Make sure to use the appropriate database for your organism (-D option, use `dante -h` for more information). Note that for DANTE in version 0.1.0 the command is slightly different, it is recommended to use the latest version. 

```bash
for dir in clust_assembly/SCL*/assembly/; do
  dante -q "${dir}contigs.fasta" -D Viridiplantae_v3.0 -o "contigs.domains.gff" -dir ${dir} -M BL80 -thsc 80 -c 6
done
```

Recommended to use GNU parallel over DANTE's -c option (only splits large files in chunks, rather than processing each folder in parallel).

```bash
parallel -j 6 '
  dir={};
  dante -q "${dir}contigs.fasta" -D Viridiplantae_v3.0 -o "contigs.domains.gff" -dir ${dir} -M BL80 -thsc 80 -c 1
' ::: clust_assembly/SCL*/assembly/;
```

Filter the DANTE output.

```bash
for dir in clust_assembly/SCL*/assembly/; do
  dante_gff_output_filtering.py -dg "${dir}contigs.domains.gff" -ouf "${dir}contigs.domains.filter.gff" -dps "${dir}contigs.domains.filter.faa" -sd All -thl 0.4 && sed -i "/^>/ s/ /#/g" "${dir}contigs.domains.filter.faa"
done
```



#### Classification using BLAST+

The retrieved NODEs can be characterized by similarity search against any desired database. Use blast with a tabular output option to use the provided script to transform to gff3 format to be used with FlexiDot later on. Output format option example:

```bash
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

You can transform (and filter) the tabular blast output to a GFF file using the provided `blast6_2_gff.py` script (see `shared_scripts` directory). Use the help function to get more information and use the provided table header from above as `--infmt` option.

```bash
python blast6_2_gff.py -h
```



### Visualization of SCL graphs and single nodes

#### SCL graphs using bandage assembly viewer

Use it according to manual with gfa files, might use domains.filter.faa to blast and get positions

#### Dotplots of individual nodes using FlexiDot

Use according to manual with domains.gff (provide, scripts for highest scoring fastas and to revise gffs for easier viz).

Here is an example of how to run flexidot using the highest scoring contigs (NODEs):

Collect the highest scoring NODEs (score = length * coverage):

```bash
time parallel -j 6 '
  dir={};
  python top_scoring_nodes.py "${dir}contigs.fasta" "${dir}contigs.top10.fasta"
' ::: clust_assembly/SCL*/assembly/;
```

Modify the DANTE annotation files (.gff) to enable the possibility to show distinct colors with flexidot (find the `modify_dante_gff_for_flexi.py` script in the `shared_scripts` directory):

```bash
time parallel -j 6 '
  dir={};
  python modify_dante_gff_for_flexi.py "${dir}contigs.domains.gff" "${dir}contigs.domains.flexi.gff"
' ::: clust_assembly/SCL*/assembly/;
```

Run the flexidot command. Adjust the options to your needs (use `flexidot.py -h` for more information).

```bash
for dir in clust_assembly/SCL*/assembly/; do
  flexidot.py  -i "${dir}contigs.top10.fasta" -p 2 -f 0 -x y -o "${dir}flexi.top10" -z 2 -E 7 -g "${dir}contigs.domains.flexi.gff" -G <.../flexidot/color_flexidot.config>
done
```

