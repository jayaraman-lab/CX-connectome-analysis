# Code for generating figure panels in the ellipsoid body (EB) section
Below is a guide for identifying which notebooks were used for which figure panels.

### Main figures and supplemental information
**EB_synapseDistributions.Rmd**: Generate visualizations of synapse locations shown in Figure 10 and Figure 10 figure supplements.
2D histograms (see section 3a below)
* Figure 10 B (ER4m)
* Figure 10 figure supplement 1 (ring neurons)
* Figure 10 figure supplement 2 (columnar neurons)
* Figure 10 figure supplement 3 (EXR neurons)

Normalized synapse density plots in slices through the EB (see section 3b)
* Figure 10 C-E
* Figure 10 figure supplement 4

**EB_connectivityGraphsAndMatrix.Rmd**: Generate connectivity graph and matrix of connectivity
* Figure 10 F: connectivity graph
* Figure 11 A: connectivity matrix
* Figure 11 D: connectivity matrix
* Figure 13 A: connectivity matrix
* Figure 13 C: connectivity graph ordered by EPG input
* Figure 13 figure supplement 1 A,B: connectivity matrix

**EB_typeComparison.Rmd**: Visualize normalized contributions of different ring neuron types to EL and EPG neurons (Figure 13 B).

**EPGmorphology.Rmd**: ...  (maybe rename since the same notebook is used for EL?)

**ExR_connectivity.Rmd**: Connectivity and similarity matrices, bar graph of partners of ExR neurons
* Figure 14 B (similarity matrices)
* Figure 14 C (connectivity matrices) TODO
* Figure 14 figure supplement 2A (similarity matrices)
* Figure 14 figure supplement 2B (bar graph) TODO
* Figure 14 figure supplement 3A,B (connectivity matrices) TODO

**ExR-EB2EXMotifs-Preparation.Rmd** and **ExR-EB2EXMotifs-Plots.Rmd**: generate plots in figure 15. The "preparation" notebook does the 
computationally intensive work of gathering the output pathways of the ExR neurons and stores the results into a local "data" folder. The "plots" 
notebook just does the plotting. The "preparation" notebook only needs to be run once.
* Figure 15 B: motifs prevalence plot
* Figure 15 C: breakdown of motifs weights
* Figure 15 D-G: motifs circular plots for 4 ExR neurons

### Supplemental information only
**....Rmd**: ...
* ...
