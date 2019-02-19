
# statTarget2 
<p align="">
 <img src="https://github.com/13479776/Picture/blob/master/statTarget_label_biocticker.png" height="200" title="statTarget 2.0">
 </p>


![GitHub release](https://img.shields.io/badge/statTarget-Good-blue.svg)
![GitHub release](https://img.shields.io/badge/releases-v%201.11.2-yellow.svg)
![GitHub release](https://img.shields.io/badge/downloads-top%2020%25-green.svg)
![GitHub release](https://img.shields.io/badge/Dependents-R%203.3.0%20-brightgreen.svg)
![GitHub release](https://img.shields.io/badge/downloads-9020/total-brightgreen.svg)



For details and latest version see https://stattarget.github.io/docs/


Package vignettes: [Vignettes](https://stattarget.github.io/docs/my-new-doc/) 


Reference manual: [Manual](https://github.com/13479776/Picture/blob/master/statTarget-manual.pdf)


Demo data: [Data](https://stattarget.github.io/docs/demo/)


Example reports: [Reports](https://stattarget.github.io/docs/demo/)


Binary Package: [statTarget2_WindowsOnly.zip](https://github.com/13479776/Picture/raw/master/statTarget_2.0.0.zip)

###  Citation

 Please cite the following article when using statTarget or QC-RFSC algorithm:
 
 Luan H., Ji F., Chen Y., Cai Z. (2018) statTarget: A streamlined tool for signal drift correction and interpretations of quantitative mass spectrometry-based omics data. Analytica Chimica Acta. dio: https://doi.org/10.1016/j.aca.2018.08.002
 
 Luan H., Ji F., Chen Y., Cai Z. (2018) Quality control-based signal drift correction and interpretations of metabolomics/proteomics data using random forest regression. bioRxiv 253583; doi: https://doi.org/10.1101/253583

### System requirements
--------------------------------------------------------------------

> Dependent on R (>= 3.3.0)

> If you did not install the R software yet,you can download R >= 3.3.0  from https://www.r-project.org


### Installation
--------------------------------------------------------------------
     
> Install “statTarget2” at the Bioconductor

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("statTarget", version = "3.9") 
    library(statTarget)
    statTargetGUI()
 
     
> Install the version of "statTarget2" at Github. copy this code into R

    library(devtools)
    
    devtools::install_github("statTarget/statTarget2")
    
    library(statTarget)
    
    statTargetGUI()


    
> IMPORTANCE: for mac PC,  XQuartz instead of X11 support should be installed for the Graphical User Interface (GUI). Download it from https://www.xquartz.org. 


> RGTK2 is a binding for R to the GTK2 library and dependent libraries, and a multi-platform package for creating graphical user interfaces. If you have any problems about RGTK2 installation, see the related installation information for R and GTK on Windows/Mac OS at https://gist.github.com/sebkopf/9405675. 


> We recommend the R 3.3.0 and RGtk2 2.20.31 for mac OS paltform. `The R 3.3.0 and RGtk2 2.20.31 sailed through the test.` If you have multiple versions of R framework installed, RSwitch  - a small GUI that allows you to switch between R versions quickly. Download it from https://r.research.att.com

