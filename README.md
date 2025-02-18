# HCMV Transcriptome Analysis Pipeline
## Author: Emil Cacayan

### Table of Contents

- **[Introduction](#introduction)**
- **[Data](#data)**
- **[RNA-Seq Transcriptomic Pipeline Workflow](#rna-seq-transcriptomic-pipeline-workflow)**
- **[Usage](#usage)**
    + **[Dependencies](#dependencies)**
    + **[Running the Pipeline](#running-the-pipeline)**
- **[Disclaimer](#disclaimer)**
- **[References](#references)**

### Introduction

Transcriptomics is the study of the sum set of all RNA transcripts in a cell, organism, or other entity that manipulates genetic information. It is distinct from genomics as it provides a real-time view of the nexus between the information stored in genetic material and how it may be expressed. 

Researchers have characterized transcriptomes across many different organisms, but as they can provide a synchronous view of an organism's response to stimuli, it is useful to do transcriptomic analysis preceding and subsequent some event or as a time-series analysis. 

One of the vectors of genetic information of great interest to researchers are viruses - both in their clinical repercussions and clinical/nonclinical applications. Transcriptomic analysis of viruses produces data that allows researchers to more comprehensively study gene expression patterns, which gives the following potentially actionable insights:

1. **Characterization of viral life cycles.**
    - Over the course of infection, viruses undergo a number of changes in its genetic expression that can be tracked using transcriptomic analysis (Sudhagar et al., 2018).
    - Tracking these changes can provide researchers valuable information about potential drug targets (Mo et al., 2016).
2. **Identification of host response pathways.**
    - While this is not one of the objectives of this particular analysis, establishing the genetic expression profile of a host cell over the course of infection can help outline important infection response genes (Zhou et al., 2021).
3. **Annotation/Characterization of viral genes.**
    - Transcriptomic analysis allows for the classification of viral genes that may cause/escalate pathogenesis (Ivanov et al., 2023).
4. **Classification/prediction of disease severity.**
    - One of the hallmark results of bioinformatic analysis is identifications of biomarkers.
    - Transcriptomic analysis of host/viral transcriptomes allows researchers to identify biomarkers to elucidate disease progression/severity (Arriaga-Canon et al., 2022).
5. **Drug development and discovery.**
    - In understanding the information above, researchers can take advantage of the molecular mechanisms of viral infection to develop drugs to disrupt virus processes or interactions between the virus and its host (Zhou et al., 2021).

### Data

In this analysis, we wish to compare transcriptomes 2-days and 6-days post-infection (dpi) of patients infected with human herpesvirus 5, also known as human cytomegalovirus (HCMV), data which has been collected and published in a study by Cheng et al., 2017.

The sequence read archive data used in this analysis is used here:

- **Donor 1 (2 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896360
- **Donor 1 (6 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896363
- **Donor 3 (2 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896374
- **Donor 3 (6 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896375

Because this analysis was performed in a UNIX environment, the `wget` command was used to download the data, e.g.:
```
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896360
```
More information about the correct subdirectory to download these files to can be found in the [Usage](#usage) section under [Dependencies](#dependencies).

These transcriptomes will be indexed against the HCMV transcriptome (NCBI accession [NC_006273.2](https://www.ncbi.nlm.nih.gov/sra/SRX2896375)) (Gatherer et al., 2011).

### RNA-Seq Transcriptomic Pipeline Workflow

This pipeline (source code contained in `wrapper.py`) takes RNA-seq viral transcriptomic information and performs the following analyses (many of the results will be printed to the `PipelineProject_Emil_Cacayan/PipelineProject.log` file - the directory and file will be created if not already present). For the sake of demonstration, this pipeline will be explained using the data listed in the [data](#data) section. 

1. **Quantification of transcripts per million (TPM) using [kallisto](https://pachterlab.github.io/kallisto/about) (Bray et al., 2016).**
    - First, a transcriptome index for HCMV (see [data](#data)) is built.
    - Alongside the creation of the transcriptome index, the coding sequence (CDS) features will be extracted using [Biopython tools](https://biopython.org) and a `cds.fasta` file will be created which contains the coding sequence (CDS) features from the genome using the RefSeq `protein_id` as the header .
    - The `cds.fasta` file will be used as the input for the kallisto index command. To the `PipelineProject.log` file will be written the number of coding sequences in the genome in the following fashion:
    ```
    The HCMV genome (NC_006273.2) has # CDS.
    ```
    where `#` is replaced with the number of CDS features. 
    - For each sample and condition, the minimum, median, mean, and maximum TPM will be calculated from the results in the `abundance.tsv` kallisto output file. The output will be written as a tab-delimited table in `PipelineProject.log` in the following fashion:
    ```
    sample  condition   min_tpm   med_tpm   mean_tpm   max_tpm
    ```
2. **Detection of differentially expressed genes using `R` package [sleuth](https://pachterlab.github.io/sleuth/about) (Pimentel et al., 2017).**
    - The output from kallisto will be used as input for sleuth to find differentially expressed genes between the two conditions. 
    - A transcript will be considered significant if the calculated false discovery rate (FDR) falls below 0.05 ($\text{FDR} < 0.05$).
    - The following details for each significant transcript will be written as a tab-delimited table to `PipelineProject.log` in the following fashion:
    ```
    target_id   test_stat   pval   qval
    ```
6. **Filtering of transcriptome reads using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (Langmead & Salzberg, 2012).**
    - To compare strain similarity between patient samples, the transcriptome reads will be assembled (we don't expect the entire genome to be produced - just enough useful information to be used in alignment tools).
    - To ensure that only viral reads are included (as sometimes viral transcriptomic data contains reads from the host), we will index the reads against the genome index. Using Bowtie2, we will create a genome index using the reference genome for the target virus (HCMV).
    - For the subsequent assembly, only reads that map to the index will be used.
    - The number of reads in each transcriptome before and after the Bowtie2 mapping will be written to `PipelineProject.log` in the following fashion:
    ```
    Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 2 (6dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 3 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 3 (6dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    ```
6. **Assembly using [SPAdes](https://github.com/ablab/spades) (Bankevich et al., 2012).**
    - Using the output from Bowtie2, SPAdes will be utilized to generate two assemblies - one for each patient/donor with a k-mer size of 77.
    - The commands sent to the SPAdes assembler can be found in `PipelineProject.log`. 
7. **Determination of the strain of each assembly using [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Camacho et al., 2009).** 
    - The pipeline will then retrieve the longest contig from each SPAdes assembly. 
    - This contig will be used as a blast+ input to query the nr nucleotide database limited to members of the *Betaherpesvirinae* subfamily. This is done by creating a local database of just sequences from the *Betaherpesvirinae* subfamily. 
    - The blast+ runs only keep the best alignment (HSP) for any single query-subject pair of sequences.
    - For the top hit runs for each assembly, the following are written to `PipelineProject.log`: subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit score, E-value, and subject title in a tab-delimited table in the following fashion:
    ```
    Donor1:
    sacc   pident    length   qstart   qend   sstart   send   bitscore   evalue   stitle

    Donor3:
    saccc   pident   length   qstart   qend   sstart   send   bitscore   evalue   stitle
    ```

## Usage

### Dependencies

To obtain the raw transcriptomic reads, the RNA-Seq data was obtained from the sequence read archive database using `wget` after setting the current working directory as the root directory (`/python_pipeline_project`):

```
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896360 -P rna_seq_raw/
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896363 -P rna_seq_raw/
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896374 -P rna_seq_raw/
wget https://www.ncbi.nlm.nih.gov/sra/SRX2896375 -P rna_seq_raw/
```

This saves the downloaded files into a directory (`/python_pipeline_project/rna_seq_raw/`) - a new directory will be created if one is not already present in the root directory. 

### Running the Pipeline



## Disclaimer

This is a project for the computational biology class (COMP 483) taught by Dr. Heather Wheeler at Loyola University Chicago. I credit much of the structure of the pipeline to the objectives outlined in her class. I strongly discourage viewing the source code if you are currently in this class or similar classes and suggest you find and implement solutions to bioinformatic analysis using your own resources. I take no responsibility for any plagiarism that occurs from this repository as a result of my posting of this repository. 

## References

Arriaga-Canon, C., Contreras-Espinosa, L., Rebollar-Vega, R., Montiel-Manríquez, R., Cedro-Tanda, A., García-Gordillo, J. A., Álvarez-Gómez, R. M., Jiménez-Trejo, F., Castro-Hernández, C., & Herrera, L. A. (2022). Transcriptomics and RNA-Based Therapeutics as Potential Approaches to Manage SARS-CoV-2 Infection. International journal of molecular sciences, 23(19), 11058. https://doi.org/10.3390/ijms231911058

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology : a journal of computational molecular cell biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. BMC bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

Cheng, S., Caviness, K., Buehler, J., Smithey, M., Nikolich-Žugich, J., & Goodrum, F. (2017). Transcriptome-wide characterization of human cytomegalovirus in natural infection and experimental latency. Proceedings of the National Academy of Sciences of the United States of America, 114(49), E10586–E10595. https://doi.org/10.1073/pnas.1710522114

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics (Oxford, England), 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

Gatherer, D., Seirafian, S., Cunningham, C., Holton, M., Dargan, D. J., Baluchova, K., Hector, R. D., Galbraith, J., Herzyk, P., Wilkinson, G. W., & Davison, A. J. (2011). High-resolution human cytomegalovirus transcriptome. Proceedings of the National Academy of Sciences of the United States of America, 108(49), 19755–19760. https://doi.org/10.1073/pnas.1115861108

Ivanov, S. M., Tarasova, O. A., & Poroikov, V. V. (2023). Transcriptome-based analysis of human peripheral blood reveals regulators of immune response in different viral infections. Frontiers in immunology, 14, 1199482. https://doi.org/10.3389/fimmu.2023.1199482

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

Mo, Z. Q., Li, Y. W., Wang, H. Q., Wang, J. L., Ni, L. Y., Yang, M., Lao, G. F., Luo, X. C., Li, A. X., & Dan, X. M. (2016). Comparative transcriptional profile of the fish parasite Cryptocaryon irritans. Parasites & vectors, 9(1), 630. https://doi.org/10.1186/s13071-016-1919-1

Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature methods, 14(7), 687–690. https://doi.org/10.1038/nmeth.4324

Sudhagar, A., Kumar, G., & El-Matbouli, M. (2018). Transcriptome Analysis Based on RNA-Seq in Understanding Pathogenic Mechanisms of Diseases and the Immune System of Fish: A Comprehensive Review. International journal of molecular sciences, 19(1), 245. https://doi.org/10.3390/ijms19010245

Zhou, A., Dong, X., Liu, M., & Tang, B. (2021). Comprehensive Transcriptomic Analysis Identifies Novel Antiviral Factors Against Influenza A Virus Infection. Frontiers in immunology, 12, 632798. https://doi.org/10.3389/fimmu.2021.632798