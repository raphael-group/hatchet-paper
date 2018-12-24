# HATCHet paper

This repository contains the simulated data, the results of all the methods considered in the comparison, and the results of HATCHet ([HATCHet repository](https://github.com/raphael-group/hatchet)) on the published whole-genome multi-sample tumor sequencing datasets, all these are described in the HATCHet paper at:

[Simone Zaccaria and Ben Raphael, 2018](https://www.biorxiv.org/content/early/2018/12/17/496174)

## Contents

1. [Simulated data](#simulateddata)
    - [Tumors](#tumors)
    - [Patients and samples](#patientsandsamples)
    - [Results](#results)
      - [Fixed](#fixed)
      - [Free](#free)
2. [Cancer data](#cancerdata)
    - [Prostate cancer](#prostatecancer)
    - [Pancreas cancer](#pancreascancer)

## Simulated data
<a name="simulateddata"></a>

All the simulated data are included in the folder `simulation`. The simulated data comprises 256 mixed samples with 2-3 tumor clones for 64 patients (3-5 samples per patient), half with a whole-genome duplication (WGD). These simulated data have been generated using MASCoTE which has also been described in the HATCHet paper and is available here:

[MASCoTE](https://github.com/raphael-group/mascote)

Due to space limitations, we are unable to publish in this repository the sequencing reads from all samples. As such, for every sample we provide a BB file which encodes the read-depth ratio (RDR) and B-allele frequency (BAF) for every genomic bin of the reference genome in every sample. The BB files are produced by the pre-prosessing steps of HATCHet and summarize whole the basic input data needed by CNA-inference methods.

### Tumors
<a name="tumors"></a>

All these data are reported in the folder `data`. The patients are divided between those with a tumor without a WGD in `noWGD` folder and those with a tumor with a WGD in `WGD` folder. In both cases, a dataset corresponds to a collection of clones and is reported in the format `dataset_nX_sY`where `X` is the number of tumor clones (in addition there is a normal diploid clone) and `Y` is the random seed used by MASCoTE for reproducibility. The copy-number profiles of the tumor clones and the corresponding phylogenetic tree with CNAs and WGds are correspondingly reported in the two following files contained in the related subfolder `tumor`:

| Filename | Description | Format |
|----------|-------------|--------|
| `copynumber.cvs` | The allele and clone-specific copy-number profiles resulting from the CNAs and WGDs simulated by MASCoTE | The file is a tab-separated file with the following fields:<ul><li>`CHR`: the name of a chromosome</li><li>`START`: the genomic position in `CHR` determining the start of a genomic segment</li><li>`END`: the genomic position in `CHR` determining the end of the corresponding genomic segment</li><li>`cloneX`: the copy-number state of `cloneX` (with `X` from 0 to N-1) in the corresponding genomic segment. The copy number state of `cloneX` is given in the format <code>A&#124;B</code> where `A` and `B` are the two allele-specific copy numbers</li></ul>
| `tumor.dot` | The phylogenetic tree describing the tumor evolution where there is a node for every clone and the edges are labeled by the corresponding CNAs and WGDs. | The phylogenetic tree is encoded in the [DOT](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) format. The mutations are given in the following formats: <ul><li>A CNA in a edge is reported in the format `(START,END) del/tdup in P/M-CHR` where `START`, `END` are the genomic coordinates of the corresponding genomic segment, `del` or `tdup` indicate whether the corresponding CNA is a deletion or duplication respectively, `M` or `P` indicates whether the maternal or paternal copy has been, and `CHR` is the corresponding chromosome.</li><li>A chromosomal arm aberration is reported in the format `(START,END) del/tdup of P/M-CHR arm` where `START`, `END` are the genomic coordinates of the corresponding chromosomal arm, `del` or `tdup` indicate whether the corresponding aberration is a deletion or duplication respectively, `M` or `P` indicates whether the maternal or paternal copy has been affected, and `CHR` is the corresponding chromosome.</li><li>A chromosomal loss is given in the format `M/P-CHR loss` where `M` or `P` indicates whether the maternal or paternal copy has been lost and `CHR` indicates the corresponding chromosome.</li><li>A WGD is reported in the format `WGD`</li></ul>|

### Patients and samples
<a name="patientsandsamples"></a>

Each dataset includes two patients and for each patient a BB file describes the RDR and BAF of every genomic bin in all samples of the corresponding patient. The name of each BB file specifies the number of samples for the related patient, the number of clones, and the corresponding clone proportions. More specifically, the BB filename is given by a `_`-separated list where the first element preceeded by the letter `k` specifies the number of corresponding samples and each other element specifies the clone proportions of a sample, listed such that the first is the proportion of normal diploid clone and the clone proportion of any other tumor clone is given in corresponding order. The name of a sample is a `_`-separated list which starts with the noun `bulk` and each element specifies the clone proportion (without the dot) of every clone.

For example, `k4_01090_02008_00506035_00504055.bb.gz` is a BB file for a patient with 4 samples which incude 2 tumor clones (`clone0` and `clone1`) and a normal diploid clone `normal`. In particular, the samples have the following clonal compositions

| Name of sample | `normal` proportion | `clone0` proportion | `clone1` proportion |
|----------------|---------------------|---------------------|---------------------|
| `bulk_01normal_09clone0_Noneclone1` | `0.1` | `0.9` | Not present |
| `bulk_02normal_Noneclone0_08clone1` | `0.2` | Not present | `0.8` |
| `bulk_005normal_06clone0_035clone1` | `0.05` | `0.6` | `0.35` |
| `bulk_005normal_04clone0_055clone1` | `0.05` | `0.4` | `0.45` |

Another example, `k7_040600_010090_020008_0103060_0205003_0100504_01030303.bb.gz` is a BB file for a patient with 7 samples which incude 3 tumor clones (`clone0`, `clone1`, and `clone2`) and a normal diploid clone `normal`. In particular, the samples have the following clonal compositions

| Name of sample | `normal` proportion | `clone0` proportion | `clone1` proportion | `clone2` proportion |
|----------------|---------------------|---------------------|---------------------|---------------------|
| `bulk_04normal_06clone0_Noneclone1_Noneclone2` | `0.4` | `0.6` | Not present | Not present |
| `bulk_01normal_Noneclone0_09clone1_Noneclone2` | `0.1` | Not present | `0.9` | Not present |
| `bulk_02normal_Noneclone0_Noneclone1_08clone2` | `0.2` | Not present | Not present | `0.8` |
| `bulk_01normal_03clone0_06clone1_Noneclone2` | `0.1` | `0.3` | `0.6` | Not present |
| `bulk_02normal_05clone0_Noneclone1_03clone2` | `0.2` | `0.5` | Not present | `0.3` |
| `bulk_01normal_03clone0_03clone1_03clone2` | `0.1` | `0.3` | `0.3` | `0.3` |

Each BB file corresponds to a patint and is a tab-separated file describing the RDR and BAF of every genomic bin in all samples in the following format:

| Field | Description |
|-------|-------------|
| `CHR` | Name of a chromosome |
| `START` | Starting genomic position of a genomic bin in `CHR` |
| `END` | Ending genomic position of a genomic bin in `CHR` |
| `SAMPLE` | Name of a tumor sample |
| `RD` | RDR of the bin in `SAMPLE` |
| `#SNPS` | Number of SNPs present in the bin in `SAMPLE` |
| `COV` | Average coverage in the bin in `SAMPLE` |
| `ALPHA` | Alpha parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from A allele |
| `BETA` | Beta parameter related to the binomial model of BAF for the bin in `SAMPLE`, typically total number of reads from B allele |
| `BAF` | BAF of the bin in `SAMPLE` |

Due to space limitations, each BB file has been compressed using `gzip` with level of compression 9. The file can be easily decompressed with the command `gzip -d BBFILE`.

### Results
<a name="results"></a>

HATCHet has been compared with 4 current state-of-the-art methods for CNA inference:

| Method | Reference | Repository |
|--------|-----------|------------|
| Battenberg | [(Nik-Zainal et al., *Cell*, 2012)](https://www.cell.com/abstract/S0092-8674%2812%2900527-2) | [cgpBattenberg](https://github.com/cancerit/cgpBattenberg) and [Wedge-Oxford Battenberg](https://github.com/Wedge-Oxford/battenberg) |
| TITAN | [(Ha et al., *Genome Research*, 2014)](https://genome.cshlp.org/content/early/2014/07/24/gr.180281.114.abstract) | [TitanCNA](https://github.com/gavinha/TitanCNA) |
| THetA | [(Oesper et al., *Genome Biology*, 2013)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-7-r80) | [THetA/THetA2](https://github.com/raphael-group/THetA) |
| cloneHD | [(Fischer et al., *Cell Reports*, 2014)](https://www.cell.com/cell-reports/fulltext/S2211-1247(14)00373-8) | [cloneHD](https://github.com/andrej-fischer/cloneHD) |

Each of these methods and HATCHet has been applied on the simulated samples. More specifically, Battenberg, TITAN, and THetA have been applied on each sample individually, cloneHD has been applied jointly on all samples from the same patient, and HATCHet has been applied both on each sample individually (single-sample HATCHet) and jointly on all samples from the same patient. We consider two different settings when executing the methods on simulated data.

#### Fixed
<a name="fixed"></a>

First, every method has been applied on all 128 samples of the 32 patients without a WGD by providing the true value of the main parameters, including tumor ploidy, number of clones, and maximum copy number. In this case, the results obtained by every method are reported in the folder `fixed` and in the subfolder of the corresponding dataset. The results of Battenberg, TITAN, and THetA are specifically reported for every sample, the results of cloneHD are reported for every patient, and the results of HATCHet are specifically reported for every sample (when obtained by executing HATCHet on each sample inidividually) and specifically for every patient (when obtained by executing HATCHet jointly on all samples from the same patient).

#### Free
<a name="free"></a>

Second, every method has been applied on all 256 samples of the 64 patients with and without a WGD, requiring that each method infers all the relevant parameters, including tumor ploidy and number of clones, and setting the maximum copy number to 8. THetA has been excluded from this analysis as it does not automatically infer the presence/absence of a WGD. In this case, the results obtained by every method are reported in the folder `free` and in the subfolder of the corresponding dataset, which are divided according to either the presence or absence of a WGD. The results of Battenberg and TITAN are specifically reported for every sample, the results of cloneHD are reported for every patient, and the results of HATCHet are specifically reported for every sample (when obtained by executing HATCHet on each sample inidividually) and specifically for every patient (when obtained by executing HATCHet jointly on all samples from the same patient).

For every method, all the most important and relevant output files are reported. The largest of these files have been compressed due to space limitations using the command `gzip -9` and they can be easily decompressed by using the corresponding command `gzip -d`.

## Cancer data
<a name="cancerdata"></a>

HATCHet has been applied on two whole-genome multi-sample tumor sequencing datasets; the first dataset comprises 10 prostate cancer patients analyzed in [(Gundem et al., *Nature*, 2015)](https://www.nature.com/articles/nature14347) and the second dataset comprises 4 pancreas cancer patients described in [(Makohon-Moore et al., *Nature genetics*, 2017)](https://www.nature.com/articles/ng.3764).

### Prostate cancer
<a name="prostatecancer"></a>

The data for all prostate cancer patients are contained in the subfolder `prostate`. For each of the 10 prostate cancer patients (A10, A12, A17, A21, A22, A24, A29, A31, A32, and A34) the results inferred by HATCHet in a subfolder with the corresponding name. More specifically, the following files encode the results inferred by HATCHet for each prostate cancer patient:

| Name | Description | Format |
|------|-------------|--------|
| `best.seg.ucn` | Clone and allele-specific copy number profiles and clone proportions for every genomic segment | the format is described in the HATCHet repository [here](https://github.com/raphael-group/hatchet/blob/master/doc/doc_hatchet.md) |
| `best.bbc.ucn.gz` | Clone and allele-specific copy number profiles and clone proportions for every clustered bin with the corresponding RDR and BAF | the format is described in the HATCHet repository [here](https://github.com/raphael-group/hatchet/blob/master/doc/doc_hatchet.md). Due to space limitations, this file is compressed |
| `chosen.diploid.seg.ucn` | The best result inferred by HATCHet assuming there is no WGD | The format is the same of `best.seg.ucn` |
| `chosen.tetraploid.seg.ucn` | The best result inferred by HATCHet assuming there is a WGD | The format is the same of `best.seg.ucn` |

The mutations inferred from all samples of every prostate cancer patient are reported in a subfolder `mutations`. The SNVs and small indels are reported in two comma-separated files `indel_hc.csv` and `snv_hc.csv` with the following fields

| Name | Description |
|------|-------------|
| `Patient` | The name of a patient|
| `Sample` | A sample from the patient `Patient` |
| `chrom` | The name of a chromosome |
| `position` | The genomic position of a somatic-point mutation in `chrom` |
| `ref` | Number of sequencing reads coverging `position` with the reference allele |
| `var` | Number of sequencing reads coverging `position` with the alternating allele, i.e. harboring the mutation |
| `normal_reads1` | Reads supporting the reference allele of `position` in the matched-normal sample (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `normal_reads2` | Reads supporting the variant allele of `position` in the matched-normal sample (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `normal_var_freq` | Variant-allele frequency of `position` in the matched-normal sample |
| `normal_gt` | Genotype call for `position` in matched-normal sample | 
| `tumor_reads1` | Reads supporting the reference allele of `position` in the tumor sample `Sample` (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `tumor_reads2` | Reads supporting the variant allele of `position` in the tumor sample `Sample` (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `tumor_var_freq` | Variant-allele frequency (VAF) of `position` in the tumor sample `Sample` |
| `tumor_gt` | Genotype call for `position` in the tumor sample `Sample` |
| `somatic_status` | Status of the variant (Germline, Somatic, or LOH). Here, all the mutations are Somatic |
| `variant_p_value` | Significance of variant read count compared to baseline error rate |
| `somatic_p_value` | Significance of tumor read count compared to normal read count |

### Pancreas cancer
<a name="pancreascancer"></a>

The data for all pancreas cancer patients are contained in the subfolder `pancreas`. For each of the 4 pancreas cancer patients (Pam01, Pam02, Pam03, and Pam04) the results inferred by HATCHet in a subfolder with the corresponding name. More specifically, the following files encode the results inferred by HATCHet for each pancreas cancer patient:

| Name | Description | Format |
|------|-------------|--------|
| `best.seg.ucn` | Clone and allele-specific copy number profiles and clone proportions for every genomic segment | the format is described in the HATCHet repository [here](https://github.com/raphael-group/hatchet/blob/master/doc/doc_hatchet.md) |
| `best.bbc.ucn.gz` | Clone and allele-specific copy number profiles and clone proportions for every clustered bin with the corresponding RDR and BAF | the format is described in the HATCHet repository [here](https://github.com/raphael-group/hatchet/blob/master/doc/doc_hatchet.md). Due to space limitations, this file is compressed |
| `chosen.diploid.seg.ucn` | The best result inferred by HATCHet assuming there is no WGD | The format is the same of `best.seg.ucn` |
| `chosen.tetraploid.seg.ucn` | The best result inferred by HATCHet assuming there is a WGD | The format is the same of `best.seg.ucn` |

The mutations inferred from all samples of every pancreas cancer patient are reported in a subfolder `mutations`. The SNVs and small indels are reported in two comma-separated files `indel_hc.csv` and `snv_hc.csv` with the following fields

| Name | Description |
|------|-------------|
| `Patient` | The name of a patient|
| `Sample` | A sample from the patient `Patient` |
| `chrom` | The name of a chromosome |
| `position` | The genomic position of a somatic-point mutation in `chrom` |
| `ref` | Number of sequencing reads coverging `position` with the reference allele |
| `var` | Number of sequencing reads coverging `position` with the alternating allele, i.e. harboring the mutation |
| `normal_reads1` | Reads supporting the reference allele of `position` in the matched-normal sample (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `normal_reads2` | Reads supporting the variant allele of `position` in the matched-normal sample (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `normal_var_freq` | Variant-allele frequency of `position` in the matched-normal sample |
| `normal_gt` | Genotype call for `position` in matched-normal sample | 
| `tumor_reads1` | Reads supporting the reference allele of `position` in the tumor sample `Sample` (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `tumor_reads2` | Reads supporting the variant allele of `position` in the tumor sample `Sample` (the corresponding fields with plus/minus are specific to reads belonging to +/- strand) |
| `tumor_var_freq` | Variant-allele frequency (VAF) of `position` in the tumor sample `Sample` |
| `tumor_gt` | Genotype call for `position` in the tumor sample `Sample` |
| `somatic_status` | Status of the variant (Germline, Somatic, or LOH). Here, all the mutations are Somatic |
| `variant_p_value` | Significance of variant read count compared to baseline error rate |
| `somatic_p_value` | Significance of tumor read count compared to normal read count |
