# To convert bcl into fastq.

snRNA-seq data is saved in bcl format, to convert it into fastq file format, two packages are required: cellranger and bcl2fastq2. The conversion command is `cellranger mkfastq`. The FASTQ output generated will be the same as when running bcl2fastq directly. [cellranger mkfastq and bcl2fastq2](https://janis.readthedocs.io/en/latest/tools/bioinformatics/cellranger/cellrangermkfastq.html)


### How does bcl2fastq2 work?
[bcl2fastq mannual](https://sapac.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf)
**BCL to FASTQ Conversion Process**
The software uses input files, which are the output of a sequencing run, to convert BCL files into FASTQ files. For each cluster that passes filter (PF), the software writes one entry to one FASTQ file for each sample in
each read.
- For a single-read run, the software creates one Read 1 FASTQ file per sample.
- For a paired-end run, the software creates one Read 1 and one Read 2 FASTQ file per sample.
The sample FASTQ files are compressed and appended with the *fastq.gz extension. Thus, per-cycle
BCL files are converted into per-read FASTQ files that can be used as input for data analysis.

**Demultiplexing Process**
Multiplexing adds a unique index adapter sequence to each sample during library prep, generating uniquely
tagged libraries that can be identified and sorted for analysis. Demultiplexing then assigns clusters to a
sample based on the index adapter sequence of the cluster.
- To optimize demultiplexing results, choose index adapters that optimize color balance when performing library prep. For more information, see the Index Adapters Pooling Guide.
- The bcl2fastq2 Conversion Software demultiplexes multiplexed samples as part of the conversion process. If
samples are not multiplexed, the software skips demultiplexing and assigns all clusters in a flow cell lane to
one sample.

**Adapter Trimming and UMI Removal**
Depending on settings, the bcl2fastq2 Conversion Software trims adapter sequences and removes Unique
Molecular Identifier (UMI) bases from reads:
- Adapter trimming—The software determines whether a read extends past the DNA insert and into the
sequencing adapter. An approximate string matching algorithm identifies all or part of the adapter
sequence and treats inserts and deletions (indels) as one mismatch. Base calls matching the adapter
sequence and beyond are masked or removed from the FASTQ file.
- UMI removal—UMIs are random k-mers attached to the genomic DNA (gDNA) before polymerase chain
reaction (PCR) amplification. After the UMI is amplified with amplicons, the software can retrieve the
bases and include them in the read name in the FASTQ files. When the TrimUMI sample sheet setting is
active, the software can also remove the bases from the reads.
My understanding here is: Adapter trimming is performed as default setting if run cellranger mkfastq or bcl2fastq. In other words, "Base calls matching the adapter
sequence and beyond are masked or removed from the FASTQ file". Wheras, UMI removal will be performed only if you give UMI trimming command.

** To think: **
- How does the fastq file look like after bcl conversion?


To convert:
```
module load cellranger/3.0.2
module load bcl2fastq2/2.19.1
module load bcl2fastq2/2.20.0

cellranger mkfastq --qc --id=Justin_TZ_mkfastq \
                   --run=/data/bioinformatics/projects/sahar2021/Tingting/2_data_Justin_cellranger_scriptTest \
                   --csv=/data/bioinformatics/projects/sahar2021/Tingting/2_data_Justin_cellranger_scriptTest/SampleSheet_modified.csv
```
