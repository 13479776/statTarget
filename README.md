# statTarget
==============
Description
-----------------
A graphical user interface, easy to use tool provide quality control based sig-
nal correction, integration of metabolomic data from multiple batches, and the
comprehensive statistic analysis for non-targeted and targeted approaches.

The main GUI of statTarget, shown in Figure 1, has two basic components.
The rst is shift correction. It includes quality control-based robust LOESS
signal correction (QC-RLSC) that is a widely accepted method for quality con-
trol based signal correction and integration of metabolomic data from multiple
analytical batches (Dunn WB., et al. 2011; Luan H., et al. 2015).

statTarget - Shift Correction (Figure 1) provide QC-RLSC algorithm that t
the QC data, and each metabolites in the true sample will be normalized to the QC sample. Additonally, LOESS based generalised cross-validation (GCV) would be automatically applied to avoid overtting of the observed data, when
the QCspan was set at 0.
Requirements
-----------------

Dependent on R (>= 3.3.0)
Packages should be installed:
randomForest,plyr,LMGene,pracma,pROC,robustbase,grDevices,graphics,stats,utils,rrcov,RGtk2,pls

Steps and Data Frame
-----------------
![github](https://github.com/13479776/Picture/blob/master/statTarget1.png "13479776")
![github](https://github.com/13479776/Picture/blob/master/statTarget2.png "13479776")

Usage
-----------------
Download the package statTarget_win_1.0.1.zip under the window PC.

install.package("statTarget_win_1.2.0.zip")

Author
-----------------

Hemi Luan, hemi.luan@gmail.com
