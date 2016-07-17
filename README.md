# statTarget
==============
Description
-----------------
A graphical user interface, easy to use tool provide quality control based sig-
nal correction, integration of metabolomic data from multiple batches, and the
comprehensive statistic analysis for non-targeted and targeted approaches.

The main GUI of statTarget has two basic components.
The rst is shift correction. It includes quality control-based robust LOESS
signal correction (QC-RLSC) that is a widely accepted method for quality con-
trol based signal correction and integration of metabolomic data from multiple
analytical batches (Dunn WB., et al. 2011; Luan H., et al. 2015).

statTarget - Shift Correction provide QC-RLSC algorithm that t
the QC data, and each metabolites in the true sample will be normalized to the QC sample. Additonally, LOESS based generalised cross-validation (GCV) would be automatically applied to avoid overtting of the observed data, when
the QCspan was set at 0.

statTarget - Statistical Analysis provide features including Data pre-
processing, Data descriptions, Multivariate statistics analysis and Univariate
analysis.
Data preprocessing : 80-precent rule, log transformation, KNN imputation, Me-
dian imputation and Minimum values imputation.
Data descriptions : Mean value, Median value, Sum, Quartile, Standard deriva-
tives, etc.
Multivariate statistics analysis : PCA, PLSDA, OPLSDA, VIP, Random forest.
Univariate analysis : Student T-test, Shapiro-Wilk normality test and Mann-
Whitney tests.
Biomarkers analysis for Clinical research : ROC, Odd ratio.


Requirements
-----------------

Dependent on R (>= 3.3.0), GTK2

Packages should be installed:
randomForest,plyr,LMGene,pracma,pROC,robustbase,grDevices,graphics,stats,utils,rrcov,RGtk2,pls

Steps and Data Frame
-----------------
![github](https://github.com/13479776/Picture/blob/master/statTarget1.png "13479776")
![github](https://github.com/13479776/Picture/blob/master/statTarget2.png "13479776")

Usage
-----------------

1 Download gtk+-bundle 2.22.1-20101229 win64.zip from
http://ftp.gnome.org/pub/gnome/binaries/win64/gtk+/2.22/ .

2 Download the package statTarget_win_1.2.0.zip under the window PC.

3 install.package("statTarget_win_1.2.0.zip")

Author
-----------------

Hemi Luan, hemi.luan@gmail.com
