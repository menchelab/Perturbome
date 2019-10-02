# The Perturbome
Disease associated processes, as well as therapeutic interventions, can be identified with perturbations within the interactome map of molecular interactions in the cell. Here, we introduce a mathematical framework for mapping out precisely how different perturbations influence each other in terms of 12 distinct interaction types. We performed a large-scale imaging screen of cells perturbed by diverse chemical compounds and identified a perturbome network of 242 perturbations and 1832 interactions. Our analysis of the chemical and biological features of the perturbing compounds reveals distinct molecular fingerprints for each interaction type. The interactome overlap between perturbations determines how they influence each other, either increasing or decreasing each others’ effects, or resulting in the emergence of entirely new phenotypes. Our framework can be applied to other key challenges, such as dissecting the combined impact of genetic variations or predicting the effect of a drug on a particular disease phenotype.  

The following code can both be used to apply some or all aspects on own data as well as for the reproduction of the individual plots.

## Installation and System Requirements
### Python packages used for analysis
Most of the data extraction and analysis was performed using python 2.7. Analysis of the core periphery structure was conducted using the cpalgorithm module in python 3.7.0 
#### Python version: 2.7.10  
Python packages included:  
- numpy [v 1.15.1]
- skimage [v 0.14.1]
- sklearn [v 0.20.1]
- networkx [v 2.1.0]
- MySQLdb [v 1.4.1]
- matplotlib [v 2.2.3]
- pandas [v 0.24.2]
- PIL [v 5.4.1]
- scipy [v 0.19.1]
- seaborn [v 0.9.0]
- mygene [v 2.3.0]
- urllib2 [v 2.7.0]
- json [v 2.0.9]
- mayavi [v 4.6.2]
- sympy [v 1.0.0]

#### Python version: 3.7.0
Python packages included: 
- cpalgorithm [v 0.0.1]

### System requirements
Code was written and execuded on a MacBookPro
- CPU: Intel Core i5
- Memory: 8GB
- Graphics: Intel Iris Graphics 6100
- MacOS: 10.14.3

### Installation instructions
All python packages should be easily being installed via pip [https://pypi.org/project/pip/]  
e.g. pip install <package_name>  
Typically time to install all packages should be between 30min and 1h depending on the amount of previously installed packages.


## Data
All data to run the individual scripts should be within the respective /data/<NameOfAnalysis> folders. Note that some of the data is redundant (i.e. same data files in various /data/<NameOfAnalysis> folders. The reason for that is to keep each analysis modular and executable without the others.
Additional information and data can be downloaded under: https://sites.google.com/view/menchelab/perturbome

Raw fluorescent images are  available under: https://idr.openmicroscopy.org/ [idr0069]

## Code execution and run time
The majority of code runs within several minutes. Only exception is the calculation of interactions [9_Calculate_Interactions.ipynb] that can take from several hours up to one day depending on the underlying system. 
The final output of the individual analysis parts matches the output of the paper: "Mapping the perturbome network of cellular perturbations", with the only exception of some heuristic code sections. 

## Code description and help
The majority of all code is written in jupyter notebook format, including detailed information in markdown about the purpose and result of the individual sections of code.  
In case anything is not clear please send an email to: mcaldera@cemm.at
