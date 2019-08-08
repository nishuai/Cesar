[![Build Status](https://travis-ci.org/lima1/PureCN.svg?branch=master)](https://travis-ci.org/lima1/PureCN)
[![Coverage](https://img.shields.io/codecov/c/github/lima1/PureCN.svg)](https://codecov.io/gh/lima1/PureCN)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0) 

# 💡 Cesar

A tool developed for tumor-only copy number estimation for targeted capture sequencing data with segmentation and anchor-based recalibration.

**一款在靶向panel中寻找特定基因拷贝数变异的软件。**

* [中文版说明文档](https://github.com/nishuai/Cesar/blob/master/docs/Casar%E6%A3%80%E6%B5%8BNSCLC%E6%82%A3%E8%80%85ctDNA%E4%B8%AD%E7%9A%84CNV.docx?raw=true)


## 🏘️ Installation

To install Cesar, simply copy the repository to your local destination

```
wget https://github.com/nishuai/Cesar/archive/master.zip
unzip master.zip && mv master Cesar
```
And make sure the `MASS` package is installed in your local environment, in R:
```
install.packages('MASS')
```

## ⚡ How to use

Cesar is designed for use in a specific target capture sequencing panel with all genes to test for CNV status. It detectets abnormal CNV status by learning from a bunch of normal samples with no CNV. 

A usual learning paradigm is to take 10-100 normal samples with site specific coverage information in a pileup format. The first 3 columns in a pileup file should be Chromosome (chr1), position (122975) and read depth (3006). However, smaller number of training samples (like 5) is also allowed. If you would like to get your hands on it before using your own data, you can try to run an example demo comes with Cesar:

### Example run for trainning Cesar:

#### Training Cesar
```
Rscript Training_anchors.R inputdata/example.bed mpileups/normal_pileups/ output_dir
```
Cesar takes a sorted and non-overlaping bed file, if you are not sure if your bed file meets the Cesar standard, you can provide `Training_anchors.R` with it and Cesar will tell you if the bed file is OK. Cesar also comes with a function to revise your bed file, simply use `R/revise_bed.R` to make sure your bed file is sorted and all checked.

To reduce perturbations in depth estimation, we suggest not to include bed regions with a span of less than 30 bp. `R/revise_bed.R` also comes with a function to remove bed regions with less than a user-specivied number of nucleotides.

#### Revise your bed file
```
Rscript R/revise_bed.R inputdata/example.bed 30 output_dir
```

`R/revise_bed` takes 1-3 arguments, if only the first one is give, it will automatically remove bed regions with span less than 30 base pairs. If output_dir is not provided, the revised bed file will be stored just beside the input bed file in the same directory.

#### Training result

Depending on the number of training samples given, the training process will typically take about some few minuts to complete (0.2 minute per sample for 2M sample). This step will generate 2 model files for Cesar.R in Rdata format. the `output_dir/model_anchors.rda` and `output_dir/model_parameters.rda`.

#### Detecting CNV using Cesar
After learning form normal samples, Cesar can be used to detect CNV in a test sample:
```
Rscript Cesar.R mpileups/test.mpileup output_dir/model_anchors.rda output_dir/model_parameters.rda outputdir
```
The bed file used here should be ensically the same as that used in training, However, in this step Cesar will automatically detect gene names at the 4th column of the bed file, trying to call gene-level CNVs. If the 4th column is not given, Cesar will simply call CNV in a global manner.

Cesar will generate a pdf file to visualize global CNV status in the sample. If the corresponding gene name is given for each bed region, it will also generate a table with genes having the most CNV. 



 

