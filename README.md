# A Proof-of-Concept Analytical Framework for Concurrent Microbial Profiling and Host Somatic Readout from Shotgun Biopsy Metagenomics in Colorectal Cancer

### Abstract
Comprehensive molecular characterization of colorectal cancer requires integrated analysis of both host and microbial features, but these are typically profiled using separate experimental approaches. Using paired tumor and tumor-adjacent biopsies from 19 Indian patients, we show that a single sequencing dataset can be leveraged to explore both the microbial community landscape and the host somatic genomic profile, offering a framework for dual-compartment analysis in colorectal cancer. Although microbial reads represented a small fraction of the total biopsy-derived sequences, species-level profiling revealed patient identity as the dominant driver of microbial variation, accounting for 59% of the observed variance. Despite this strong inter-individual effect, Enterocloster citroniae and Staphylococcus hominis emerged as significant discriminators of tumor tissue, suggesting their potential value as tissue-associated microbial markers. In parallel, host-derived reads enabled exploratory analysis of somatic alterations from the same metagenomic libraries. These data revealed a mutation landscape enriched in intergenic regions, recurrent alterations in KMT2C, and copy-number gains on chromosomes 7, 19, and 1. While host coverage was limited, the findings indicate that biopsy shotgun metagenomes can provide complementary information beyond microbiome profiling alone. Together, this study presents a proof-of-concept for using colorectal biopsy metagenomic data to interrogate both microbial and host compartments, highlighting a cost-effective strategy for future integrative cancer studies and population-specific biomarker discovery.

Our Analysis is Distributed in three major parts

1. Raw read quality processing and reference alignment against the hg38 genome to reconstruct host sequences from unaligned microbial reads.
2. Perform taxonomic and functional microbial metagenomic profiling followed by diversity and pathway analysis on the unaligned reads.
3. Identify high confidence human somatic nucleotide variants and copy number variations from the aligned reads.

Here we've overview of design of study and workflow

![Graphical Abstract](GraphicalAbstract.jpeg)
![Workflow](Workflow.tiff)
