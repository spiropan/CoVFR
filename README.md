This repository contains code and processed data used for the paper 
"COVID-19 vaccination and age-stratified all-cause mortality" available at 
https://www.researchgate.net/publication/355581860_COVID_vaccination_and_age-stratified_all-cause_mortality_risk

The InputFiles subfolder contains spreadsheets downloaded from the CDC that were filtered and sorted by State abbreviation. These spreadsheets were imported into Matlab and used in downstream analyses. 
The proc_all.m file contains all MATLAB code used to import, process and analyze data from the spreadsheets. 

The Tables subfolder contains monthly vaccination, COVID cases, and Y20, Y21 age-stratified death counts. These Tables are written out by the Matlab code (first cell of the proc_all.m file).
Readers who would like to inspect and replicate the results or run different variations of the analyses may find it easier to double check these tables against the original CDC data and then work off of these tables rather than start from the files in the InputFiles subfolder. 

The Results subfolder contains subfolders named after each model type, which in turn contain .csv files with results for each outcome variable used.   

The Figures subfolder contains subfolders named after each model type, which in turn contain plots labeled according to month, outcome variable and age group.

For links to the original data sources, please see the Data Sources section in our manuscript (see link above).  

Please email me or open an issue with any questions or requests for clarifications.  
