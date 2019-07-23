[![Build Status](https://travis-ci.org/lima1/PureCN.svg?branch=master)](https://travis-ci.org/lima1/PureCN)
[![Coverage](https://img.shields.io/codecov/c/github/lima1/PureCN.svg)](https://codecov.io/gh/lima1/PureCN)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0) 

# Cesar

A tool developed for tumor-only copy number estimation for targeted capture sequencing data with segmentation and anchor-based recalibration.

## Installation

To install Cesar, simply copy the repository to your local destination

```
wget https://github.com/nishuai/Cesar/archive/master.zip
unzip master.zip && mv master Cesar

```


## How to use

Cesar is designed for use in a specific target capture sequencing panel with limited number of genes to test for CNV status, if one needs to train Cesar for their own panel, go to line 122 in Training_anchors_clinic.R and change the target regions for CNV detection, the numbers in line 122 corresponds to line number in the given bed file. 

The training takes a directory with at least 5 processed pileup files, the pileups shoud be ending with *.freq. Only the first 3 columns is used for Cesar, the first 3 columns should be Chromosome, position and read depth. If you are all set, you can try with an example comes with Cesar:

### Example run for trainning Cesar:
```
Rscript Training_anchors_clinic.R inputdata/segmented_bed_30k.bed mpileups/normal/ output_dir
```
Depending on the number of training samples give, the training process will typically take about some few minuts to complete. This step will generate 2 model files for Cesar.R as input. the output_dir/model_anchors.rda and output_dir/model_parameters.rda.

After learning form normal samples, Cesar can be used to detect CNV in abnormal samples:
```
Rscript Cesar.R mpileups/met2.1875ERBB2.625/703-13.sort.bam.mpileup.freq output_dir/model_anchors.rda output_dir/model_parameters.rda outputdir
```
This will generate a file called {input_sample_name}.cnv. Where it contains the CNV states for genes we specified.



[说明文档](https://github.com/nishuai/Cesar/blob/master/docs/Casar%E6%A3%80%E6%B5%8BNSCLC%E6%82%A3%E8%80%85ctDNA%E4%B8%AD%E7%9A%84CNV.docx?raw=true)).




