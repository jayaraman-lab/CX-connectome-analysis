# Outputs figure guide
The notebooks in this folder cover figures 54 to 65 and figure 74.

**Outputs-Core.Rmd** does not generate any figures, but builds the objects that are going to be used by
all other notebooks. Those objects are stored in a local `data` folder which is accessed by the other noteboooks.
It is fairly computationally intensive (execution of all the cells takes over 1h30mins) but only needs to be ran once.

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




