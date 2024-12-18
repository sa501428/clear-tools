# CLEAR Tools

Chromatin Landscape Exploration, Analysis, and Research Tools

# Main Tools

## APA2

Faster version of APA specifically for larger loop lists (>10k loops)

## ATA:

Aggregate Track Analysis (ATA) is designed for analyzing genomic signal data around defined DNA loci/peaks. This tool
reads signal data from a BigWig file, narrow peak data from a BED file, and aggregates the signal values (e.g.,
ChIP-seq) around the midpoints of peaks to generate a normalized output. The output is saved in a .npy format for
downstream analysis.

```
ata [--res int] <signal.bw> <peaks.bed> <outfile> <genome>

java -jar clear-tools.jar ata --window 1000 --res 1 signal.bw peaks.bed output.npy hg38
```

This aggregates the signal data from `signal.bw` around peaks in `peaks.bed` and saves the results as `output.npy`,
using the `hg38` genome.

`--res` int: (Optional) Resolution of signal aggregation. Defaults to `1` (and primarily designed and tested for 1 bp
analysis)
`signal.bw`: input BigWig file containing signal data.
`peaks.bed`: input BED file containing peak regions.
`outfile`: Output file prefix for aggregated results (saved as .npy).
`genome`: Genome assembly (e.g., `hg19`, `hg38`) used to map chromosomes.

## FUSE

Combines multiple bedpe files into one.
Flags available for duplicate removal, NMS, combining, and more.
Attributes do not need to be in the same order between the files.
Simplest usage is just to combine multiple bedpe files together without worrying about order of attributes etc.

## SPLIT

Split up a bedpe into multiple lists.
E.g. will split a file into 10 by randomly assigning features to 1 of the 10 outputs.
Useful when you want to parallelize a task on a cluster but need to split up the bedpe to do it.
Can use FUSE to combine the results into a final file.
Specifying a number <= 0 will split the bedpe by chromosome
(i.e. will make a new bedpe for each chromosome, and all features of that chromosome will be together
in the split file)

## INTERSECT

Intersect two bedpe files.
Flag available to subtract the files instead of intersecting them.
Also has flags for using exact matches, overlapping boundaries, etc.

## EXPAND

Expands the size of anchors for each loop.
Flags available for whether to explicitly set a size or only expand if the width
is currently smaller than the requested width.

# Dev or Deprecated Tools

## PINPOINT

Peak Intensity (or Identification) Near Points Of INTeraction

## FILTER

## FLAGS

Focal Loop Aggregation via Grid Search

## ENHANCE

Examining Neighborhoods at High-resolution by Aggregating Nearby Contact Events

### old name: amplifi

Aggregating Mixed Peaks Leading to Improved Focal Intensity

## RECAP

Ratios and Enrichments Compared Across Phenotypes

## COMPILE

COMParing Interactions from Loops across Experiments

## HOTSPOT

Highly Observable Transitions in Selective Pixels Of Tissues

## SIFT

Search and Identify Focal Targets
Simple Identification of Focal Targets

## SIEVE

Selecting Interactions Enriched Versus Environment

## SIMPLE

Simple Identification of Maximum Point of Local Enrichment

## additional tools

probability

clean

CAMP
Comparisons Across Multiple Phenotypes

ADAM
analyze differences across maps

SLASH
Summarizing Loop Aggregate Statistics and Hierarchies

TODO

ACE
Annotating Cliques of Enrichment

CORE
Cliques Of Regional Enrichments