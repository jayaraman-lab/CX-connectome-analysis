# Central complex connectome analysis
This repository contains the R notebooks used to generate the figures in the paper: **A connectome of the Drosophila central complex reveals network motifs suitable for flexible navigation and context-dependent action selection**, Brad K. Hulse, Hannah Haberkern, Romain Franconville, Daniel B. Turner-Evans, Shinya Takemura, Tanya Wolff, Marcella Noorman, Marisa Dreher, Chuntao Dan, Ruchi Parekh, Ann M. Hermundstad, Gerald M. Rubin, Vivek Jayaraman.

The code has been written by Hannah Haberkern, Brad Hulse, Dan Turner-Evans, Romain Franconville and Marcella Noorman.
This code has not been peer reviewed.

## Prerequesites
You will need a recent (>4.0) installation of R. Depending on the situation, you could chose between two installation routes:
### Minimal, manual installation

If you just plan to browse some of the notebooks and pick what you need, we recommend you simply install the [**neuprintrExtra**](https://github.com/jayaraman-lab/neuprintrExtra) package which was developed for this project:
 ```r
  # install
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jayaraman-lab/neuprintrExtra")

# use 
library(neuprintrExtra)
  ```
   as this package and its dependencies are used by most notebooks. You will likely have to install required additional packages as you make your way through different notebooks. 

### Full, reproducible installation
If you plan to run all notebooks, or prefer not to worry with manual installations, or you are going through this code years after its release, you can recreate the environment we used thanks to the [**renv**](https://rstudio.github.io/renv/articles/renv.html) package.

## Organization of the folders
The repository is organized by folders, each folder covering either a section of the paper (in most cases) or a type of analysis. Each folder contains an individual *Readme* file that describes the content of the notebooks.

### Folders related to paper sections:
**fiberTractsFigures** contains the code related to figure 1 -- figure supplement 2 and 3 that summarize fiber tract diameters of some of the main CX types

**validation** contains the code related to figures 3 and 4, about the evaluation of the effect of the level of proofreading on the results

**Inputs** contains code related to figures 5 to 9 and 25 to 27, that describe input pathways into the CX, generally (figure 5), through the antero-ventral pathway (AVP, figures 6 to 9), and through the NO (figures 25 to 27).

**EBFigures** contains code related to figures 10 to 15, 17 and 18, which describe ellipsoid body networks.

**PBFigures** contains code related to figures 20 to 24, which describe PB connectivity and neuron types.

**FBFigures** contains code related to figures 29 to 53, which describe the columnar organization and intra connectivity of the FB (figures 29 to 39), the FB tangential neurons (figures 40 to 45), the pathways linking the mushroom body to the FB (figures 46 and 47) and the circuits related to sleep (figures 48 to 53) 

**Outputs** contains code related to figures 54 to 65 and 74, which describe central complex output pathways

**DiscFigures** contains code related to figures 70, 73 and 75. It runs the bump propagation simulations of figures 70 and 73 and generates the graph representation of figure 75.

### Folders related to particular analysis:
**NeuropileInnervationPlots** contains the code generating all the neuropile innervation plots throughout the paper (starting with figure 7A)

### Other folders:
**colormaps** contains consistent colormaps reused in multiple figures

**R** contains .R files of functions reused throughout the analysis
