# GarnerKindLai_MeTu
GarnerKindLai_MeTu paper
https://www.biorxiv.org/content/10.1101/2023.11.29.569241v1
(Currently under revision.)

## OS
All EM analyses were performed in MacOS using R and Python. All calcium imaging analyses were performed in Windows using Matlab.

## Required packages

### R
Natverse	https://natverse.org	DOI: 10.7554/eLife.53350	collection of R packages for neuroanatomical analysis

Tidyverse	https://www.tidyverse.org/

alphashape3d	https://CRAN.R-project.org/package=alphashape3d	 

### Python 3

numpy https://numpy.org/install/

pandas https://pandas.pydata.org/docs/getting_started/index.html

matplotlib https://matplotlib.org/stable/

seaborn https://seaborn.pydata.org/installing.html

pydantic https://docs.pydantic.dev/latest/install/

scipy https://scipy.org/install/

navis https://navis.readthedocs.io/en/latest/source/install.html

caveclient https://caveclient.readthedocs.io/en/latest/guide/intro.html

cloudvolume https://github.com/seung-lab/cloud-volume

pymaid	https://github.com/schlegelp/PyMaid

FAFBseg	https://github.com/flyconnectome/fafbseg-py

bokeh https://docs.bokeh.org/en/latest/docs/first_steps.html#first-steps

holoviews https://www.holoviews.org/

hvplot https://hvplot.holoviz.org/getting_started/installation.html

neuprint https://connectome-neuprint.github.io/neuprint-python/docs/

### Matlab

Computer Vision Toolbox

Image Processing Toolbox

Optimization Toolbox

Statistics and Machine Learning Toolbox


## Installation

All necessary packages for R, Python, and Matlab should be installed. Installation instructions can be found at homepages of those packages. We used the most recent packages as of March, 2024. Installation time for these packages depend on the computer specification. Please see the homepage of each packages.

There is no required hardware other than a computer that can run R, Python, and Matlab (most recent versions as of March, 2024).

Git clone or download the current repository. Cloning or downloading time is typically less than 1 minute.


## Running the code.

EM analyses: Run the main script in the folder to replicate analyses and produce figures. Each section of code analyze the data and generate figures. The analyses time varies widely (some analyses, such as scanning the entire database to identify visual pathways to the central complex, may take a few hours).

Two-photon calcium imaging: Run plot_script files in each of the figure folder to generate plots in the figures. (Raw imaging files will be provided upon request. Intermediate data files are available in each folder.)


## Other notes

Users must have a FlyWire account and an API token to run EM analyses code. Visit https://flywire.ai for more details.





