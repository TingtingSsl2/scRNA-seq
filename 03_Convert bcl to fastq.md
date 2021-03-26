# To convert bcl into fastq.

snRNA-seq data is saved in bcl format, to convert it into fastq file format, two packages are required: cellranger and bcl2fastq2.

To convert:
```
module load cellranger/3.0.2
module load bcl2fastq2/2.19.1
module load bcl2fastq2/2.20.0

cellranger mkfastq --qc --id=Justin_TZ_mkfastq \
                   --run=/data/bioinformatics/projects/sahar2021/Tingting/2_data_Justin_cellranger_scriptTest \
                   --csv=/data/bioinformatics/projects/sahar2021/Tingting/2_data_Justin_cellranger_scriptTest/SampleSheet_modified.csv
```
