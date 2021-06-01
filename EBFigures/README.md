# Code for generating figure panels in the ellipsoid body (EB) section
Below is a guide for identifying which notebooks were used for which figure panels.

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

**modularityAnalysis/generateModularityTable.m**: Compute wedge-specific modularity of connectivity btw ring neurons and EPGs. Running generateModularityTable() in MATLAB will generate the plots in Figure 11 figure supplement 1 C and will save the p-values to a csv file. This function relies on two other functions (included in the folder '/modularityAnalysis/'): (1) the function computeModularity.m computes the modularity of inputs from (or outputs to) a given ring neuron type; this is performed both for the measured connectivity matrix and for many shuffled versions of the connectivity matrix; (2) the function Shuffle.m (https://www.mathworks.com/matlabcentral/fileexchange/27076-shuffle) speeds up the process of shuffling connectivity matrices.
* Figure 11 figure supplement 1 C (histograms and similarity matrices)

**EB_typeComparison.Rmd**: Visualize normalized contributions of different ring neuron types to EL and EPG neurons (Figure 13 B).

**EB_electrotonicDistances.Rmd**: Analysis of distance between synapses and putative spike initiation zone for EPGs and ELs in EB. Generates all plots in Figure 12 and its supplements.
* Figure 12 A: example EPG and its rootpoint

CDFs of distance between synapse and rootpoint grouped by modality
* Figure 12 C (normalized electrotonic distance, EPG)
* Figure 12 supplement 1E (raw distance along arbors, EPG)

Medians of distance distributions
* Figure 12 E (normalized electrotonic distance, EPG)
* Figure 12 supplement 3B (normalized electrotonic distance, EL)
* Figure 12 supplement 1F (raw distance along arbors, EPG)
* Figure 12 supplement 2 (normalized electrotonic distance, by neuron type, EPG)

Rank ordering of medians
* Figure 12 supplement 1B (EPG)
* Figure 12 supplement 3C (EL)

Median distance between distributions
* Figure 12 supplement 1C (EPG)

Distribution of synapses with mean distance contours
* Figure 12 C left (normalized electrotonic distance, EPG)
* Figure 12 supplement 3A top (normalized electrotonic distance, EL)
* Figure 12 supplement 1D (raw distance along arbors, EPG)
* Figure 12 C right (synapse density by modality, EPG)
* Figure 12 supplement 3A bottom (synapse density by modality, EL)
* Figure 12 supplement 1A (synapse density by modality, plotted separately, EPG)

**ExR_connectivity.Rmd**: Connectivity and similarity matrices, bar graph of partners of ExR neurons
* Figure 14 B (similarity matrices)
* Figure 14 C (connectivity matrices)
* Figure 14 figure supplement 2A (similarity matrices)
* Figure 14 figure supplement 2B (bar graph)
* Figure 14 figure supplement 3A,B (connectivity matrices)

**ExR-EB2EXMotifs-Preparation.Rmd** and **ExR-EB2EXMotifs-Plots.Rmd**: generate plots in figure 15. The "preparation" notebook does the
computationally intensive work of gathering the output pathways of the ExR neurons and stores the results into a local "data" folder. The "plots" notebook just does the plotting. The "preparation" notebook only needs to be run once.
* Figure 15 B: motifs prevalence plot
* Figure 15 C: breakdown of motifs weights
* Figure 15 D-G: motifs circular plots for 4 ExR neurons

**EB Columnar Figure.RMD**: Generate connectivity matrices and graph
* Figure 17 C: connectivity matrix
* Figure 17 D: connectivity matrix
* Figure 17 E: connectivity graph laid out as a force-directed network

**EB-EBtFigure.Rmd**: Generate connectivity matrices and plot the total number of synapses per ROI per neuron
* Figure 18 B: connectivity matrix
* Figure 18 C: number of synapses per ROI per neuron
* Figure 18 D: connectivity matrix
* Figure 18 E: connectivity matrix