## Introduction

This is the code base for my Masters thesis, titled _Sectoral dependence and contagion in the BRICS grouping: an application of the R-vine copulas_. The study provides a contribution to the literature that distinguishes between contagion and interdependence between financial assets by proposing the application of the R-Vine methodology in testing the extent of correlation between the extreme joint distribution of the different markets or economies. 

### R and RStudio versions

Ensure that you have the latest R version installed. The version information for R for this project is:  
version.string R version 3.5.1 (2018-07-02)
nickname: Feather Spray
platform: x86 64-w64-mingw32
arch: x86 64 
os: mingw32 
system: x86_64, mingw32
major: 3
minor: 5.1
year: 2018
month: 07
day: 02
svn rev: 74947
language: R


The RStudio version is Version 1.1.463

### Project setup

A .Rproj file is used for the project setup. From RStudio, go to 'File', 'Open Project...' and open the .proj file accociated with this project. RStudio will open the necessary files.
If you are not using RStudio but base R instead, set your current working directory to where the main.R file is within the project.

### Library dependencies

The libraries that are imported into the project are:
```r
library(dplyr)
library(MTS)
library(psych)
library(tseries)
library(rugarch)
library(VineCopula)
library(dplyr)
library(kdecopula)
```
### Running the project

Simply run main.R. All the relevant results will be placed into the finalResults folder.
