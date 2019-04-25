# Ch4-Opakapaka-Growth
Analysis and Manuscript Associated with Dissertation Chapter 4

Estimating Opakapaka Growth Parameters using Bayesian and Maximum Likelihood frameworks
This repository contains code for the manuscript 'Revised Growth Estimates for Pristipomoides Filamentosus in the Hawaiian Islands using mark-recapture studies, bayesian analysis, and integrative data approaches'
Project goals were to fit vonBertalanffy growth parameters to mark-recapture data collected by HDAR's Opakapaka Tagging Program. 
There are two primary implementations of this:
	1) A JAGS based Bayesian approach
	2) A maximum likelihood based integrative data approach
The later approach also includes direct aging data from Ralston and Miyamoto (1983), DeMartini et al. (1994), Andrews et al. (2012) as well as length-frequency data collected by Moffitt & Parrish (1996). 

Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Prerequisites

Prior to running code in this repository, you will need to download and install the following software:
	- R (https://www.r-project.org)
	- JAGS (http://mcmc-jags.sourceforge.net)

It is recommended you also install RStudio (http://rstudio.com)

Once R is installed and running, there are several packages that must be installed
	- doParallel
	- mixtools
	- R2jags
	- lattice
	- coda
	- ggplot2
	- forcats
	
These can be installed using the R command install.packages('package_name') by replacing package_name with the name of the relevant package.

After installing the relevant packages, all analysis can be performed on the user's local machine, though use of a high performance cluster is strongly suggested as maximum likelihood models can take literally weeks to run.

Bayesian and maximum likelihood analyses are handled by two files. These can be found in the src folder within the project directory. 

Bayesian analysis handled by the notebook file: "Bayesian VBGF Fitting.Rmd"
Maximum likelihood handled by the script file "Integrative_Growth_Model_Analysis.R"

R notebooks handle relative paths weird when knitting, prior to running this file, change the variable "proj_dir" in the first chunk to reflect your path to the "Analysis" folder within the main project directory.

After this is completed, script/notebook files should run their respective analyses.


Author(s)
Stephen Scherrer

License
This project is licensed under the MIT License
All rights preserved, all wrongs traversed


Acknowledgments
Frank Parrish, Robert Moffitt, Stephen Ralston, Garret Miyamoto, Allen Andrews, Edward DeMartini, Jon Brodziak, Ryan Nichols, and Robert Humphreys were involved with prior growth studies used in integrative maximum likelihood modeling
Annette Tagawa for provided the OTP mark recapture data used in this analysis. 
We would also like to thank Zane Zhang and Paige Eveson for providing code used to fit Bayesian and maximum likelihood models.

