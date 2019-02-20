# DEG-RNA-seq

## Basic Usage
<p align="center">
<img src=" " alt="" width="600" height="650">
</p>

<div align="center">
  <sub>The pipline is build by Amina Bedrat at harvard T.H.chan School of public health
  (<a href="http://lemoslab.org">Lemos Lab</a>)
  </a>
</div>
 
![PERL](https://img.shields.io/badge/perl-v5.18.2-blue.svg)
[![License](https://img.shields.io/badge/license-GNU_v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.fr.html)


## Features

- Written in  perl and shell.
- No installation necessary.
- To use, you need to change the paths of the run file.
- Works on Mac and Linux (never tested on Windows).

## Needed files
- 
- 

## Used tools 
- bamUtil
- BWA 
- samtools
- GATK
- perl

## Usage
The pipeline runs on two steps:

### Alignenent and depth calculation
#### To run the first step:

```bash
$ cd PATH/TO/RUN1.sh
$ RUN1.sh
```
This step will: 
- 
- 

After this step, X files are generated. They are used through the second step.

### Copy number calculation:
#### To run the second step:

```bash
#YOU NEED TO CHANGE THE:
#UUID =>Line 4
#Specify the path of the results output folders (line 9)
#THE PATH OF THE 14 FILES (from line 20 to 35)
#Know you can run 

$ cd PATH/TO/RUN2.sh
$ RUN2.sh
```
This step will:
- 


## Results:
The rDNA copy number is summerized in the file rDNA_CN/CN-sampleID.txt

```
rDNA_subunit	Background_Depth	rDNA_Depth	rDNA_CN	Type

```

## Contributing

#### Bug Reports & Feature Requests

Please use the [issue tracker](https://github.com/karan/joe/issues) to report any bugs or file feature requests.


#### Contact

This is a handy script that automates a lot of developing steps.
