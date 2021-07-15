# Code for generating figure panels in the fan-shaped body (FB) section
Below is a guide for identifying which notebooks were used for which figure panels.

**FB_GetAndSave_Synapses.Rmd**:  Does not generate any figures but loads and preprocesses synapse location data, which gets stored to a configurable directory and used by several other notebooks

**FB_ColumnWidth.Rmd**: Plots PB-FB-X columnar neuron arbor width and inter-column distance 
* Figure 29 D: Arbor width 
* Figure 29 E: Inter-column distance

**FB_ColumnLocations.Rmd**: Plots the median location of all FB columnar neurons 
* Figure 29 C: PFNp_a, PFNa, PFGs
* Figure 29 figure supplement 1: All PB-FB-X types
* Figure 31: vDelta and hDelta overview (vDeltaA_a, vDeltaB, vDeltaH, hDeltaK, hDeltaA, hDeltaH).
* Figure 31 figure supplement 1: All vDelta types
* Figure 31 figure supplement 2: All hDelta types 
* Figure 32: FX overview (FR1, FS1A, FC1E)
* Figure 32 figure supplement 1: All FR/FS types
* Figure 32 figure supplement 2: All FC types  

**FB_PhaseShifts.Rmd**: Plots the discrete PB-FB phase shifts and PB-FB to intra-FB connectivity matrices
* Figure 30 Aiii: Discrete phase shifts for PFGs and PFR_a
* Figure 30 Biii: Discrete phase shifts for PFNp_a and PFNa
* Figure 34 B: Discrete phase shifts for PFNp_a and PFNa
* Figure 34 C: Connectivity matrices showing that phase shifts impact FB connectivity
* Figure 39 Aiii: Discrete phase shifts for PFL2
* Figure 39 Biii: Discrete phase shifts for PFL1
* Figure 39 Ciii: Discrete phase shifts for PFL3

**FBFigure-ColConnectivity.Rmd**: Connectivity graphs and matrices, function to display one, two, or three step connections, connectivity thresholds, dendrograms, and clustering
* Figure 33 A: Connectivity graph
* Figure 33 B: Step counting between types
* Figure 33 C: Connectivity matrix
* Figure 33 figure supplement 1 A: Connectivity matrix
* Figure 33 figure supplement 1 B: Number of types in a connectivity matrix as a function of the % of total neurons within that type that must be connected
* Figure 33 figure supplement 2 A: Dendrograms
* Figure 33 figure supplement 2 B: Clustering by cosine distance
* Figure 33 figure supplement 2 C: Connectivity graph grouped by clusters
* Figure 33 figure supplement 3 B: Connectivity matrix

**FB_ColumnAngles.Rmd**: Quantitative analysis of FB columnar phase and PB-FB phase shift magnitude
* Figure 34 D: Anatomical phase as a function of each neuron's medial-lateral location
* Figure 34 E: PB-FB phase shift magnitude
* Figure 34 F: PFN phase shift histogram according to number of glomeruli sampled
* Figure 34 figure supplement 1 C: PFNa and PFNp_c phase shift histograms according to number of glomeruli sampled

**FB_Lateralization.Rmd**: Analyzes whether left and right PB-FB types target distinct neurons/types in FB (i.e. lateralization)
* Figure 35 B: Scatter plot showing input from left and right PB-FB-X populations for postsynaptic FB types
* Figure 35 D: Scatter plot showing input from left and right PB-FB-X populations for postsynaptic FB neurons
* Figure 35 E: Scatter plot of connection consistent versus lateralization for PB-FB-X to intra-FB connection matrices

**FBFigure-AB.Rmd**: Connectivity matrices and mean number of neurons per FB column or per neuron for a given neuron type
* Figure 36 D: Connectivity matrix
* Figure 36 E: Mean number of synapses per FB column per type
* Figure 36 F: Connectivity matrix
* Figure 36 figure supplement 1 C: Mean number of synapses per neuron
* Figure 36 figure supplement 1 D: Connectivity matrix

**FB_Motifs.Rmd**: Analysis of intra-FB connectivity motifs
* Figure 37 B: Connectivity matrices showing example motifs
* Figure 37 C: Scatter plot of PCA space of connectivity matrices
* Figure 37 figure supplement 1 D: Scatter plot of PCA space of connectivity matrices along with additional example connectivity matrices
* Figure 37 figure supplement 2 A: Plot of variance explained by each principal component
* Figure 37 figure supplement 2 B: Matrices showing first 6 PCs along with proportion of variance explained by each

**FBFigure-InfoGoesUp.Rmd**: Plot the density of synapses for a given type (or types) in a chosen ROI
* Figure 38 A: 2D synapse density histogram
* Figure 38 C: 2D synapse density histogram

**FBFigure-TanAnalyses.Rmd**: Type counting analyses (per ROI and overall numbers), connectivity matrices, and synapse location plotting
* Figure 40 C: Number of types that received input or give output in a given ROI
* Figure 40 D: Number of neurons per types
* Figure 41 B: Connectivity matrix
* Figure 41 C: Connectivity matrix
* Figure 42: Connectivity matrix
* Figure 43: Connectivity matrix
* Figure 44 A: Connectivity matrix
* Figure 44 B: Connectivity matrix
* Figure 45 A: Connectivity matrix
* Figure 45 B: Synapse location

**FB_Mbon_to_Cx.Rmd**: Analysis of MBON to CX connectivity
* Figure 46 A: Network graph showing direct MBON to CX connections
* Figure 46 figure supplement 1 A: Scatter plot comparing relative weight to raw weight for MBON-to-CX connections
* Figure 46 figure supplement 1 B: Bar graph showing each MBON's strongest connection to the CX
* Figure 46 figure supplement 1 C: Histogram of MBON-to-CX connection strengths
* Figure 47 A: Matrix showing the number of FB neurons targeted (per layer) by each MBON (threshold = 0.01), through at most one interneuron
* Figure 47 B: Matrix showing the number of FB neurons targeted (per layer) by each MBON (threshold = 0.02), through at most one interneuron
* Figure 47 C: Network graph showing indirect connections (through one interneuron) from visual MBONs to FB tangential neurons
* Figure 47 D: Network graph showing indirect connections (through one interneuron) from thermo/hygro MBONs to FB tangential neurons
* Figure 47 E: Network graph showing indirect connections (through one interneuron) from MBONs to non-FB tangential CX neurons
* Figure 47 figure supplement 1 A: Connectivity matrix from MBON to interneurons targeting CX types ('first hop')
* Figure 47 figure supplement 1 B: Connectivity matrix from interneurons targeting CX types to those CX types ('second hop')

**FB_SleepSection.Rmd**: Analysis of dorsal FB sleep-wake network 
* Figure 48 C: Matrix showing cosine similarity of connectivity for neurons contained in R23E10 
* Figure 50 A: Neuron-to-neuron connectivity matrix of sleep-wake types
* Figure 50 B: Network graph showing type-to-type connectivity of sleep-wake types
* Figure 51 A: Connectivity matrices showing downstream targets of sleep-wake types
* Figure 51 B: Bar plot showing the number of synapses onto downstream types from sleep-wake types 
* Figure 51 C: Bar plot showing the number of types downstream of sleep-wake types
* Figure 52 A: Connectivity matrices showing upstream inputs to sleep-wake types
* Figure 52 B: Bar plot showing number of synapses onto sleep-wake types from upstream neurons 
* Figure 52 C: Bar plot showing the number of types upstream of sleep-wake types
* Figure 53: Network graph showing major connections bewtween the dFB and EB 
