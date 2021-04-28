# Outputs figure guide
The notebooks in this folder cover figures 54 to 65 and figure 74.

**Outputs-Core.Rmd** does not generate any figures, but builds the objects that are going to be used by all other notebooks. Those objects are stored in a local `data` folder which is accessed by the other noteboooks. It is fairly computationally intensive (execution of all the cells takes over 1h30mins) but only needs to be run once.

**Outputs-Figure-CXOuts.Rmd** generates the plots in figure 54:
- figure 54 A: Modified neuron innervation plot of the putative CX output neurons
- figure 54 B: 3d rendering of output neuron synapse locations
- figure 54 C: Synapse counts for the putative CX output neurons

**Outputs-Figure-DivConv.Rmd** generates the plots in figure 55:
- figure 55 A: diagram of pathway flow
- figure 55 B: quantification of divergence/convergence/recurrence
- figure 55 C,D: composition of output pathways per layer
- figure 55 E: CX neurons reach

**Outputs-Figure-GA-ROB-RUB.Rmd** generates the plots in figure 56 and supplements:
- figure 56 A: 3d rendering of the CX to CX fraction contributed by individual output synapses
- figure 56 B: CX target composition for neurons of the GA/BU/ROB/RUB
- figure 56 D: connectivity matrix of Gall types
- figure 56 F: connectivity graph of the same types in the GA and the EB
- figure 56 figure supplement 1 A: Full connectivity matrix in the GA
- figure 56 figure supplement 2 A: Connectivity matrix of PFR->PFR connections
- figure 56 figure supplement 2 B: Network diagram of strong PFR outputs
- figure 56 figure supplement 3 A: Connectivity matrix of FR->FR connections
- figure 56 figure supplement 3 B: Connectivity matrix of FR outputs
- figure 56 figure supplement 3 C: Network diagram of strong FR outputs

**Outputs-Figure-CX2CXMotifs.Rmd** generates the plots in figure 57 and supplements:
- figure 57 A: CX target composition in the LAL/WED/PS/SMP/CRE
- figure 57 B: Pathway connectivity of CX to CX connections (excluding FBt neurons)
- figure 57 C: Network diagrams of PFL to NO pathways
- figure 57 figure supplement 1 A: CX to CX total pathway weights by kind of connection
- figure 57 figure supplement 1 B: Pathway connectivity of all CX to CX connections

**Outputs-Figure-CX2CXMotifs.Rmd** generates the plots in figure 58:
- figure 58 B,C,D: example circular motif plots
- figure 58 E: Motifs prevalence

**Outputs-Figure-Modularity.Rmd** generates the plots in figure 59 and supplements:
- figure 59 A: Total pathway weights contributed by CX type
- figure 59 B: Distribution of pathway weights
- figure 59 C: Pathway weight connectivity table to strong targets
- figure 59 D: Downstream neuropil innervation of CX targets
- figure 59 E: 3d rendering of synapse locations of main CX targets
- figure 59 figure supplement 1: Clustering of CX types by their output connectivity at different path length
- figure 59 figure supplement 2, panel A,B: graph of the main CX targets
- figure 59 figure supplement 2, panel C: Connectivity matrix between modules of the output network

**Outputs-Figure-PathsToKnown.Rmd** generates the plots in figure 60:
- figure 60 A,B: Fraction of output pathways reaching known types, per CX output type
- figure 60, panel C: Pathway weights to the main known targets of the CX

**Outputs-Figure-MBONsPPL.Rmd** generates the plots in figure 61 and supplements:
- figure 61 A: Pathway weights to the main MBON, DAN and Antennal lobe targets
- figure 61 B: Network graph of pathways the strongest of those targets
- figure 61, figure supplement 1, A: Network graph of pathways to OviIN and MBON27

**Outputs-FigureVisual.Rmd** generates the plots in figure 62 and supplements:
- figure 62 A: Pathway weights to the main VPN targets
- figure 62 B,D,F,H: Network graphs of pathways to example targets
- Selection of neurons for figure 62 C,E,G,I
- figure 62, figure supplement 1: Network graph of pathways from PFR_b to VPNs
- figure 62, figure supplement 2 A: Connectivity matrix of PFL3 to LC33 connections
- figure 62, figure supplement 2 B: Connectivity matrix of PFL3 and LC33 targets

**Outputs-Figure-DNs.Rmd** generates the plots in figure 63 and supplements:
- figure 63 A: Pathway weights to the main DNs targets
- figure 63 B: Network graph of the pathways to the strongest DNs targets
- figure 63, figure supplement 1 A,C,E: Network graphs of extra DN pathway examples

**Outputs-Figure-PFL3.Rmd** generates the plots in figure 64:
Neuron to neuron connectivity matrix of PFL3 direct strong targets

**Outputs-Figure-FS4A.Rmd** generates the plots in figure 65:
Neuron to neuron connectivity matrix of FS4A direct strong targets

**Outputs-Figure-FirstSteps.Rmd** is used to help make figure 74:
- Panels inspiring figure 74
- Figure 74, figure supplement 1: network graph of PFL1 strong target