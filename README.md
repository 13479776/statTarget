# statTarget 

![GitHub release](https://img.shields.io/badge/statTarget-Good-blue.svg)
![GitHub release](https://img.shields.io/badge/releases-v%201.6.0-yellow.svg)
![GitHub release](https://img.shields.io/badge/downloads-top%2020%25-green.svg)

![GitHub release](https://img.shields.io/badge/Dependents-R%203.3.0%20-brightgreen.svg)

```
statTarget2 was released at https://stattarget.github.io/docs/

```
==============

**Citation

Please cite the following article when using statTarget or QC-RFSC algorithm:

Luan H., Ji F., Chen Y., Cai Z. (2018) statTarget: A streamlined tool for signal drift correction and interpretations of quantitative mass spectrometry-based omics data. Analytica Chimica Acta. dio: https://doi.org/10.1016/j.aca.2018.08.002

Luan H., Ji F., Chen Y., Cai Z. (2018) Quality control-based signal drift correction and interpretations of metabolomics/proteomics data using random forest regression. bioRxiv 253583; doi: https://doi.org/10.1101/253583

**Frequently Asked Questions**
```
Q3) ERROR: RGTK2 installation error for Mac OS paltform

    The R 3.3.0 and RGtk2 2.20.31 sailed through the test. We recommend the R 3.3.0 and RGtk2 2.20.31 for mac OS paltform. If you have multiple versions of R framework installed, RSwitch - a small GUI that allows you to switch between R versions quickly. Download it from https://r.research.att.com

Q2) Error: in imsamFP[, 2:ncol(imsamFP)]: subscript out of bounds

    In the statTarget version 1.4.X, the 1st column name of Profile Files is prefixed and should be labelled as "name".

Q1) Error: in knnimp.internal(...), NA/NaN/Inf in the foreign function call (arg 1)

    Only unique value or character could be typed into the 1st column of Profile Files and labelled as the metabolites ID.

```
**News Reports**
```
News  2017-4-1: [Development version 1.5.7] (http://bioconductor.org/packages/devel/bioc/html/statTarget.html "version 1.5.7")  have been released.

The result of Fold Change will be calculated accroding to the raw data in version 1.4.12 and dev_v1.5.7 (following the missing value estimation step).

2017-1-21: [Development version 1.5.6] (http://bioconductor.org/packages/devel/bioc/html/statTarget.html "version 1.5.6")  have been released.

transX() is to generate statTarget inputs from Mass Spectrometry Data softwares, like XCMS.
transX() directly read the .tsv file from diffreport function in XCMS software.

2017-1-15: [Release version 1.4.10] (http://bioconductor.org/packages/release/bioc/html/statTarget.html "version 1.4.8")  have been released.

2016-12-17: [Development version 1.5.5] (http://bioconductor.org/packages/devel/bioc/html/statTarget.html "version 1.5.2")  have been released. (statTarget Alert added)

2016-11-12: [Development version 1.5.1] (http://bioconductor.org/packages/devel/bioc/html/statTarget.html "version 1.5.1")  have been released. (Bugs fix)
```
Description
-----------------

An `easy to use tool` provides `graphical user interface` for quality control based `signal correction`, `integration of metabolomic data` from multi-batch experiments, and the `comprehensive statistic analysis` in non-targeted or targeted metabolomics.


Link to Bioconductor: http://bioconductor.org/packages/devel/bioc/html/statTarget.html

The Manual: http://www.bioconductor.org/packages/devel/bioc/vignettes/statTarget/inst/doc/statTarget.html

The main `GUI of statTarget` has two basic components. The first is signal correction. It includes `quality control-based signal correction` that is a widely accepted method for removal of inter-batch and intra-batch unwanted variations and integration of large-scale metabolomic data from multiple analytical batches (Dunn WB., et al. 2011; Luan H., et al. 2015).

`statTarget - Shift Correction` provide QC-RLSC algorithm that fit the QC data, and each metabolites in the true sample will be normalized to the QC sample. Additionally, LOESS based generalised cross-validation (GCV) would be automatically applied to avoid overfitting of the observed data, when the QCspan was set at 0 (Default value).

`statTarget - Statistical Analysis` provide features including Data preprocessing, Data descriptions, Multivariate statistics analysis and Univariate analysis.


Data preprocessing : 80-precent rule, glog transformation, KNN imputation, Median imputation and Minimum values imputation.


Data descriptions : Mean value, Median value, Sum, Quartile, Standard derivatives, etc.


Multivariate statistics analysis : PCA, PLS-DA, VIP, Permutation test, S-plot, Random forest.


Univariate analysis : Student T-test, Shapiro-Wilk normality test and Mann-Whitney tests. 


Biomarkers analysis : ROC, Odd ratio, P-value, and Volcano plot.
 

Requirements
-----------------

Dependent on R (>= 3.3.0)

Packages should be installed:

randomForest,plyr,pROC,rrcov,RGtk2,pls,gWidgets2,gWidgets2RGtk2,pdist,impute

Steps and Data Frame
-----------------
![github](https://github.com/13479776/Picture/blob/master/statTarget.png "13479776")

Usage
-----------------

1 If you did not install the R software yet,you can download R >= 3.3.0  from https://www.r-project.org

2 Install the package "statTarget" at the Bioconductor
 
  For Windows PC, copy this code into R 
  
  > source("https://bioconductor.org/biocLite.R") 
  
  > biocLite("statTarget")
  
  > library(statTarget) 
    
  > statTargetGUI() 
    
  
  Copy this code into R
  
  > source("https://bioconductor.org/biocLite.R")
  
  > biocLite("statTarget")
  
  > library(statTarget)  ## `Load statTarget`
  
  > statTargetGUI()  ## `Execute statTarget GUI` 
  

4 Input data and run. 



    *GTK+ installation
  
    For Win PC, If GTK+ is not available, a notice will be shown. please press 'OK'.
  
    For mac PC,  XQuartz instead of X11 should be installed. Download it from https://www.xquartz.org. if the you have problems of RGtk2 installation, *R 3.3.3* and  *RGtk2 2.20.31* will be fine. Download RGtk2 2.20.31 from https://cran.r-project.org/bin/macosx/mavericks/contrib/3.3/RGtk2_2.20.31.tgz.
  

***
See the work flow: Shift Correction, `example(shiftCor)` or Statistical Analysis, `example(statAnalysis)`


transX() is to generate statTarget inputs from Mass Spectrometry Data softwares, like XCMS.


transX() directly read the .tsv file from diffreport function in XCMS software.
***


Tutorial
-----------------

Download the [statTarget tutorial](https://raw.githubusercontent.com/13479776/Picture/master/%20statTarget_tutorial%20.pptx "statTarget tutorial .pptx") and [example data](https://github.com/13479776/Picture/blob/master/Data_example.zip "Data_example.zip") .


Author
-----------------

Hemi Luan, luanhm@sustc.edu.cn, hemi.luan@gmail.com


References
-----------------
Dunn, W.B., et al.,Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. Nature Protocols 2011, 6, 1060-83.

Luan H., LC–MS-Based Urinary Metabolite Signatures in Idiopathic Parkinson’s Disease. J Proteome Res., 2015, 14 (1),467–478.

Luan H., Non-targeted metabolomics and lipidomics LC–MS data from maternal plasma of 180 healthy pregnant women. GigaScience 2015 4:16
