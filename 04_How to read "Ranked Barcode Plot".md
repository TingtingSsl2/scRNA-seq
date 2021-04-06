# How to read 10X Cellranger "Ranked Barcode Plot"
Below is showing a **Ranked Barcode Plot**, obtained after running 10X Genomics Cellranger's `cellranger count` function. 
![example](RankedBarcodePlot.png)

- **y-axis**: number of UMI counts mapped to each barcode.

- **x-axis**: number of barcodes below "that" UMI counts. 

- The steep drop-off in plot seperates cell-associated barcodes (on the left of x-axis) and empty partitions (on the right of x-axis). 

- Barcodes can be determined to be cell-associated based on their UMI count or by their RNA-profile. In other words, if there is enough UMI counts or RNAs associated with the barcode, this barcode is cell-associated. 

# How to interpret the "Fraction Reads in Cells" metric?
[How to interpret the "Fraction Reads in Cells" metric?](https://kb.10xgenomics.com/hc/en-us/articles/360003919491-How-to-interpret-the-Fraction-Reads-in-Cells-metric-)
**Fraction Reads in Cells**
The fraction of valid-barcode, confidently-mapped-to-transcriptome reads with cell-associated barcodes.

Question: I see a low value for the "Fraction Reads in Cells". How can I interpret this metric?

Answer: A low "Fraction Reads in Cells" value is typically explained by the following:

1) High ambient RNA (background) in your sample. This ambient RNA comes from lysed/dead cells in your sample. Cell Ranger is able to confidently align the reads from ambient RNA to the transcriptome but the reads are not associated with a valid cell-containing GEM.

2) The cell-calling heuristic did not apply. For example, there may be higher variation in RNA content than expected (more cells with lower RNA content). The current cell-calling heuristic assumes a ten-fold variation in RNA content.

Cell Ranger's algorithm for partitioning barcodes as cells versus background is based on the idea that barcodes for cells should have distinctly more transcript counts associated with them than the background barcodes. This can be visualized by the ranked barcode plot in the web_summary.html file. More details on the cell filtering algorithm can be found here.

If you suspect that Cell Ranger's cell calling algorithm did not work well for your sample, please re-run cellranger count again or cellranger reanalyze with --force-cells option to call the expected number of cells.

** The flow of FASTQ files.
- bcl files are converted to FASTQ files, this step generates FASTQ files with valid Index + Adapter sequences, and Undetermined FASTQ (sequences without valid Index + Adaptoer sequences).
- the determined FASTQ files will be mapped the the genome to generate count matrix. Ambient RNA (with Index + Adapter sequence) will also be mapped to the genome. But, these ambient RNA lack of cell-associated barcodes. 







