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

The training usually requires learning from 10-100 normal samples with site specific coverage information in a pileup format. The first 3 columns in a pileup file should be Chromosome (chr1), position (122975) and read depth (3006). To get your hands on Cesar, you can try to run an example demo comes with Cesar:

### Example run for trainning Cesar:
```
Rscript Training_anchors.R inputdata/example.bed mpileups/normal_pileups/ output_dir
```
Depending on the number of training samples given, the training process will typically take about some few minuts to complete. This step will generate 2 model files for Cesar.R in Rdata format. the `output_dir/model_anchors.rda` and `output_dir/model_parameters.rda`.

After learning form normal samples, Cesar can be used to detect CNV in a test sample:
```
Rscript Cesar.R mpileups/test.mpileup output_dir/model_anchors.rda output_dir/model_parameters.rda outputdir
```
This will generate a file called `{input_sample_name}.cnv`. Where it contains the CNV states for genes we specified.




 

