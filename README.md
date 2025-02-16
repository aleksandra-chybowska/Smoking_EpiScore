# Blood- and brain-based epigenome-wide association studies of smoking
Scripts accompanying "Blood- and brain-based epigenome-wide association studies of smoking" paper.

Five folders containing annotated scripts to perform the following five main analyses described in the manuscript:

1. Brain-based EWAS of smoking (Brain)

2. Predictor of pack years of smoking (EpiScore)

3. Bayesian EWAS of pack years of smoking (BayesR) and high resolution EWAS of smoking (TWIST, ONT)

4. GWAS of pack years of smoking (GWAS)

Additionally, the repository contains a demo package.

5.  Demo/replication package with a small simulated dataset to demonstrate the code (Replication)

## System Requirements

Demo package can be executed on a standard computer. It was tested on a following machine:

RAM: 16 GB
CPU: 4 cores @ 2.6GHz
OS: Windows

Other scripts were executed on a high performance computing cluster with no "non-standard" hardware. It is not advisable to run them on a local machine ("normal" computer). The following specs are advisable for optimal performance:

RAM: 256 GB
CPU: 32 cores @ 2.7GHz
OS: Linux

All code is written in R (version 4.3.1) or python (version 3.9.7). Packages used and versions detailed in manuscript.
In addition, BayesR+ version 0.1.0 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7280277/) and LDSC version 1.0.1 was used.

## Installation Guide

Demo package: Please install R. The demo package scripts contain installation commands, which should take care of any missing libraries. Typical install and run time on a "normal" computer is less than 15 minutes.

Other scripts: Please follow installation instructions of BayesR+ (https://github.com/ctggroup/BayesRRcmd/blob/master/README.md, version 0.1.0), 
Genome-wide Complex Trait Analysis (GCTA) version 1.93.2 Linux (https://yanglab.westlake.edu.cn/software/gcta/) 
and LDSC, version 1.0.1 (https://github.com/bulik/ldsc). 


## How to Run the Demo

To run the demo, please go to Replication/scripts/ and execute the scripts in order dictated by their names. Further instructions are available in Replication/How_to_run_the_demo.docx.

## MIT License

MIT License

Copyright (c) 2024 Aleksandra Daria Chybowska

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.