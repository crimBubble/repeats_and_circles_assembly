## Assembly of eccDNA candidates from clustering

### Aim

Construct comprehensive circular sequences from ECCsplorer (clustering module only; RepeatExplorer2 (RE2)) candidate super-clusters (SCL).



### Preliminary work

In advance to the following work-flow the ECCsplorer pipeline (clustering modul with RepeatExplorer2) needs to be run including data preparation and quality control. For details on how to run see the [ECCsplorer tutorial](https://github.com/crimBubble/ECCsplorer/blob/master/tutorials/Mini-workshop.md).

Navigate to the output of the clustering module (`.../eccpipe_results/clustering_results`), run `ls -ltr` and your terminal should show something like this:

```bash
.../eccpipe_results/clustering_results$ ls -ltr
total 2792
-rwxrwxrwx+ 1 user user    1553 Jul  3  2020 HOW_TO_CITE.html
-rwxrwxrwx+ 1 user user 1478491 Jul  3  2020 contigs.fasta
-rwxrwxrwx+ 1 user user  118906 Jul  3  2020 SUPERCLUSTER_TABLE.csv
-rwxrwxrwx+ 1 user user  163985 Jul  3  2020 supercluster_report.html
-rwxrwxrwx+ 1 user user  291043 Jul  3  2020 tarean_report.html
-rwxrwxrwx+ 1 user user  212460 Jul  3  2020 COMPARATIVE_ANALYSIS_COUNTS.csv
...
-rwxrwxrwx+ 1 user user     959 Jul  3  2020 COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv
...

```



### Dependencies

- python3 (> 3.6, tested with 3.9.16 and 3.8.8)
- seqtk (tested with 1.3-r117 and 1.3-r106)
- seqkit (tested with 2.4.0 and 2.3.0)
- bowtie2 (tested with 2.4.1 and 2.2.5)
- Unicycler (tested with 0.5.0)
- DANTE (optional, tested with Galaxy/tool-shed 1.1.0 and Anaconda 0.1.8)
- GNU parallel (optional but recommended, tested with 20230522)



### Data collection for repeat assembly

The circle assembly uses clustered read data as input. First reads from super-clusters in fastq format will be collected. Secondly the reads are extended by original sequencing data (usually only a sub set of the original data is used for clustering).

#### Filter the eccDNA candidates (optional)

Sometimes it is helpful to remove unwanted eccDNA candidates before the assembly to save some work and time. In most cases the enriched organelle DNA is not the subject of interest and eccDNA candidates predicted as organelle DNA can be filtered out. Filter by automated annotation (super_CL_best_hit) 

**Note**: Manual refinement of the annotation before filtering is recommended:

```bash
awk -F"\t" '$5 !~ /organelle/' < COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv > COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES_filter.csv
```



#### 1. SCL reads (w/ quality information)

Use the reads from super-clusters in fastq format for circle assembly. 

Provide the location of your original read files in fastq format from which the SCL reads can be collected.

```bash
# Directory
export SEQDIR="<path/to/read/files/>"
# File names
export READ_1="<file1.fastq>"
export READ_2="<file2.fastq>"
```

Copy the `collect_ecc-SCL_reads.py` in to the clustering output folder (same folder that also contains `seqclust` directory and `COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv`).

```bash
cp /.../collect_ecc-SCL_reads.py /.../eccpipe_results/clustering_results/
```

Run the `collect_ecc-SCL_reads.py` script. By default, it uses the `COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES_filter.csv` file. If your file is named differently e.g., you did not filter, change the input file name in the py-script (`line 9: in_file="..."`).

``` bash
python3 collect_ecc-SCL_reads.py
```

Optional: Remove the prefix added within the ECCsplorer pipeline to match the original readnames.

```bash
for dir in ecc_assembly/SCL*/; do
sed -i 's/^[^_]*_//' "${dir}reads.lst"
done
```

Use the created read lists to extract the reads from your original reads in fastq format.

```bash
for dir in ecc_assembly/SCL*/; do
seqtk subseq "${SEQDIR}${READ_1}" "${dir}reads.lst" > "${dir}reads.for.fastq"
seqtk subseq "${SEQDIR}${READ_2}" "${dir}reads.lst" > "${dir}reads.rev.fastq"
done
```

or use parallel to run multithreaded (Note: GNU parallel needs to be installed, adjust the number of cpu cores by modifying -j INT, recommended)

```bash
parallel -j 6 '
  dir={};
  seqtk subseq "${SEQDIR}${READ_1}" "${dir}reads.lst" > "${dir}reads.for.fastq"
  seqtk subseq "${SEQDIR}${READ_2}" "${dir}read.lst" > "${dir}reads.rev.fastq"
' ::: ecc_assembly/SCL*/;
```



#### 2. SCL reads based on contigs (w/ quality information)

Use the super-cluster contigs to retrieve reads in fastq format (extensive) for eccDNA assembly. This can be used in addition to the super cluster reads from above.

Cluster reads will be collected by mapping reads in fastq format against the cluster contigs (bowtie2). The reads to be used in the further analysis will be extracted from the mapping (samtools).

Define the location (complete path, and file names) of the read files with quality information to be used. This might be sub-sampled fastq reads or the original (quality trimmed) fastq reads (Note that using the original reads takes more time overall but might harvest better results).

```bash
export SEQDIR="<path/to/read/files/>"
export READ_1="<file1.fastq>"
export READ_2="<file2.fastq>"
```

Copy the `collect_ecc-SCL_contigs.py` in to the clustering output folder (same folder that also contains `seqclust` directory and `COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES.csv`).

```bash
cp /.../collect_ecc-SCL_contigs.py /.../eccpipe_results/clustering_results/
```

Run the `collect_ecc-SCL_contigs.py` script. By default, it uses the `COMPARATIVE_CLUSTER_TABLE_eccCANDIDATES_filter.csv` file. If your file is named differently e.g., you did not filter, change the input file name in the py-script (`line 7: in_file="..."`).

``` bash
python3 collect_ecc-SCL_contigs.py
```

Create bowtie index files from the collected SCL contigs.

```bash
for dir in ecc_assembly/SCL*/; do
  bowtie2-build "${dir}collected.contigs.fasta" "${dir}collected.contigs.fasta"
done
```

Map the original reads in fastq format against all SCL contigs (adjust the number of cpu cores by modifying -p).

```bash
for dir in ecc_assembly/SCL*/; do
  bowtie2 -x "${dir}collected.contigs.fasta" -1 $"${SEQDIR}${READ_1}" -2 $"${SEQDIR}${READ_2}" -S "${dir}scl_contigs_mapped.sam" --local --sensitive-local --no-unal -p 6
done
```

Extract the mapped reads from the alignment files. Run one of the following code blocks depending on if you want to use multiple cpu cores (again, using parallel).

```bash
for dir in ecc_assembly/SCL*/; do
  samtools view -b -F 4 "${dir}scl_contigs_mapped.sam" |
  samtools fastq -1 "${dir}contigreads.for.fastq" -2 "${dir}contigreads.rev.fastq" -s "${dir}contigreads_single.fastq" -
 done
```

```bash
parallel -j 6 '
  dir={};
  samtools view -b -F 4 "${dir}scl_contigs_mapped.sam" |
  samtools fastq -1 "${dir}contigreads.for.fastq" -2 "${dir}contigreads.rev.fastq" -s "${dir}contigreads.single.fastq" -
' ::: ecc_assembly/SCL*/;
```



#### Combine reads

Combine reads from fastq and contig-fastq and remove duplicates (by read name).

```bash
for dir in ecc_assembly/SCL*/; do
  seqkit seq -i "${dir}reads.for.fastq" | cat - "${dir}contigreads.for.fastq" | seqkit rmdup -n > "${dir}merged.reads.for.fastq" && seqkit seq -i "${dir}reads.rev.fastq" | cat - "${dir}contigreads.rev.fastq" | seqkit rmdup -n > "${dir}merged.reads.rev.fastq"
done
```

Optional: Clean up your workspace. Some read files might be large depending on your raw data. 

```bash
for dir in ecc_assembly/SCL*/; do
  rm "${dir}reads.for.fastq" "${dir}contigreads.for.fastq" "${dir}reads.rev.fastq" "${dir}contigreads.rev.fastq" "${dir}scl_contigs_mapped.sam"
```



### EccDNA candidate assembly using Unicycler

Assemble the reads for each eccDNA candidate (super cluster (parts)) using the [Unicycler](https://github.com/rrwick/Unicycler) assembler (for additional options see the Unicycler manual, adjust cpu cores used with -t).

Make a new directory to keep things organized.

```bash
for dir in ecc_assembly/SCL*/; do
mkdir "${dir}cand_assembly"
done
```

Adjust the file name according to the reads you want to use (`reads, contigreads, merged`). The following code assume merged reads to be used. To use reads from super clusters or super cluster contigs change the commands accordingly.

```bash
for dir in ecc_assembly/SCL*/; do
  unicycler -1 "${dir}merged.reads.for.fastq" -2 "${dir}merged.reads.rev.fastq" -o "${dir}cand_assembly" --spades_path spades --mode normal -t 6 --min_fasta_length 1 --keep 2
done
```

#### Custom search for eccDNA candidate assemblies

Unicycler gives an indication whether circular sequences were assembled by the FASTA description "circular=true". However, when the sequence is very fragmented the circular sequences are visible in the bandage graphs but not present as complete sequence. In the next steps we try to find circles from these bandage graphs. You can also have a look at the graphs (see description below "Visualization").

Extract all sequences with the "circular=true"-tag as found by Unicycler from the assembly output.

```bash 
for dir in ecc_assembly/SCL*/cand_assembly/; do
seqkit grep -n -r -i --pattern "circular=true" "${dir}assembly.fasta" > "${dir}assembly.circ-true.fasta"
done
```

Filter bandage graphs (remove dead-ends) and save newly identified circular sequences. The new files can be recognized by there endings `.filter.gfa` and `.circles.fasta`. This script might be stored anywhere on you machine and uses some parameters (for more details see `--help` option).

```bash
for dir in ecc_assembly/SCL*/cand_assembly/; do 
python ../../scripts/find_circles_from_gfa.py -i "${dir}assembly.gfa"
done
```

**Note**: The last script as well as the "circular=true"-tag only provide an attempt to automatize the assembly process and the detection of resulting circular sequences. It is very much recommended to also take a look at the other bandage graph files (.gfa) to not miss any potential circles (use the build-in blast and depth range functionality ind bandage-viewer). Further it is recommended to test different inputs for the Unicycler (fastq-reads, contig-based-reads, merged reads). From our experience different inputs may result in different results especially if candidates share sequences partially (repeats).



__Happy eccsploring of your candidates__



### Classification of assembled circles

Suggestions for the classification of circular sequences. Note: These are just suggestions and by no means a complete overview.


#### Classification using DANTE (RexDB)

Use DANTE according to it`s manual. It is recommended to run DANTE only on circles of interest. If you want to run it automated do so as described in "Create comprehensive consensus sequences from repeat clusters" [here](..%2Fclust_assemble%2FREADME.md).


#### Classification using BLAST+

The retrieved `.circles.fasta` can be characterized by similarity search against any desired database. Use blast with a tabular output option to use the provided script to transform to gff3 format to be used with FlexiDot later on. Output format option example:

```bash
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

You can transform (and filter) the tabular blast output to a GFF file using the provided `blast6_2_gff.py` script (see `shared_scripts` directory). Use the help function to get more information and use the provided table header from above as `--infmt` option.

```bash
python blast6_2_gff.py -h
```


### Visualization of SCL graphs and single nodes


#### SCL graphs using bandage assembly viewer

Use it according to manual with ``.gfa`` files, you might use the ``<name>.domains.filter.faa`` file to blast and get the protein domain positions using the built-in blast option. The main results are presented in the ``<name>.filter.gfa`` file however also consider looking at the other ``.gfa`` files to gain a complete overview (sometimes additional circles can be found)


#### Dotplots of individual nodes using FlexiDot

Here is an example of how to run FlexiDot:

Modify the DANTE annotation files (.gff) to enable the possibility to show distinct colors with flexidot (find the `modify_dante_gff_for_flexi.py`-script in the `shared_scripts` directory):

```bash
time parallel -j 6 '
  dir={};
  python modify_dante_gff_for_flexi.py "${dir}contigs.domains.gff" "${dir}contigs.domains.flexi.gff"
' ::: clust_assembly/SCL*/assembly/;
```

Run FlexiDot:

```bash
  flexidot.py  -i <.../ecc_assembly/.../circles_of_interest.fasta> -p 2 -f 0 -x y -o circ -z 2 -E 7 -g <contigs.domains.flexi.gff> -G <.../flexidot/color_flexidot.config>
```

