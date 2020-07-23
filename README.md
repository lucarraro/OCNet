# OCNet
An R-package to generate and analyze Optimal Channel Networks

[![Linux Build Status](https://travis-ci.com/lucarraro/OCNet.svg?branch=master)](https://travis-ci.com/lucarraro/OCNet.svg?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/OCNet)](http://CRAN.R-project.org/package=OCNet) 
[![DOI](https://zenodo.org/badge/219014909.svg)](https://zenodo.org/badge/latestdoi/219014909)

## Overview

OCNet enables the creation and analysis of Optimal Channel Networks (OCNs). These are oriented spanning trees (built on rectangular lattices made up of square pixels) that reproduce all scaling features characteristic of real, natural river networks. As such, they can be used in a variety of numerical and laboratory experiments in the fields of hydrology, ecology and epidemiology. 

OCNs are obtained by minimization of a functional which represents total energy dissipated by water flowing through the network spanning the lattice. Such a formulation embeds the evidence that morphological and hydrological characteristics of rivers (in particular, water discharge and slope) follow a power-law scaling with drainage area. 

For further details, please see Carraro et al. (2020). Generation and application of river network analogues for use in ecology and evolution. *Ecology and Evolution*. doi:10.1002/ece3.6479.

## A minimal working example

Set the random seed to 1 and create an OCN in a 30x20 lattice with default options:
```
set.seed(1)
OCN <- create_OCN(30,20)
````
Draw the so-obtained OCN:
```
draw_simple_OCN(OCN)
````
![OCN 30x20](/inst/extdata/OCN_3020.png)

## Installation

```
# install devtools (if previously not installed)
if (!("devtools" %in% installed.packages())) {install.packages("devtools")}

# install OCNet from GitHub
devtools::install_github("lucarraro/OCNet", build_vignettes = TRUE)
```

## Installation issues and workarounds

### Windows

An error might occur when version `0.100.47` of package `rgl` is installed. Installation with `rgl_0.100.30` works fine. 

### Linux

Installing packages `rgdal` and `rgl` (imported by OCNet) gives rise to errors. This can be solved by running 

```
sudo apt install libftgl2 libcgal-dev libglu1-mesa-dev libglu1-mesa-dev apt-get libx11-dev libfreetype6-dev libgdal-dev 
```  

## Authors

Luca Carraro (maintainer), Florian Altermatt, Emanuel A. Fronhofer, Reinhard Furrer, Isabelle Gounand, Andrea Rinaldo, Enrico Bertuzzo
