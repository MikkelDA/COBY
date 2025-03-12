# COBY Tutorial

For a detailed introduction to the software, please refer to the [COBY Documentation](../COBY_Documentation.pdf).

For a quick reminder of available commands, you can use the [Cheat Sheet](../COBY_Cheat_Sheet.pdf).

Three tutorials are available as Jupyter Notebooks: basic, advanced, and a tutorial creating the showcased systems from the manuscript.

Tutorial ToC:

## **1. Basic Tutorial**

* **General command explanation**: A simple POPC membrane with a protein in solvent (water + 0.15 M NaCl)

<img src="figures/GeneralCommandExplanation.png" width="500">

* **Membranes**
    * **Change APL**
    * **Complex symmetric membrane**
    * **Complex asymmetric membrane**: consisting of POPC, POPE, and CHOL with different APL values for each leaflet

<img src="figures/Membranes3_ComplexAsymmetricMembrane.png" width="500">

* **Proteins**
    * **Move and rotate**
    * **Centering methods**
    * **Protein with ligands**
* **Solvation**
    * **Changing the ion concentration**
* **Box Types**
    * **Hexagonal box**
    * **Skewed hexagonal box**
    * **Dodecahedron**

## **2. Advanced Tutorial**

* **Multiple arguments**
* **Membranes**
    * **Phase separation**: featuring a checkerboard pattern of membrane patches

<img src="figures/Membranes4_PhaseSeparation.png" width="500">

*    * **Monolayers**: consisting of two monolayers with solvent bridging the headgroups and the vacuum bridging the tails

<img src="figures/Membranes5_Monolayers.png" width="500">

*    * **Lipid optimisation**
* **Nanodisc**: a disc with DMPC lipids contains within the protein rings

<img src="figures/Nanodisc.png" width="500">

* **Pores (holes)**
    * **A membrane with a solvated hole**
    * **A membrane with an unsolvated hole**
    * **Multiple holes**
    * **Polygons**: A membrane with multiple holes defined as polygonal shapes

<img src="figures/Holes4_Polygons.png" width="500">

* **Membrane Patches**
    * **A circular patch with the solvent in the membrane plane**
    * **A circular patch without the solvent in the membrane plane**
    * **Multiple membrane patches**

<img src="figures/Patches3_MultiplePatches.png" width="500">

* **Membrane holes and patches**
    * **Phase separation**
    * **Modifying shapes**
    * **Matryoshka membrane**: membrane concentric circles with a protein in the centre

<img src="figures/HolesAndPatches3_MatryoshkaMembrane.png" width="500">

* **Stacked membranes**
    * **Three bilayers**

<img src="figures/StackedMembranes1_ThreeBilayers.png" width="500">

* **Solvation**
    * **Changing the solvent concentration**
    * **Charge neutralisation methods**
    * **Mixed solvent**: A membrane solvated with regular and small water beads in specified ratios
    
<img src="figures/Solvation4_MixedSolvent.png" width="500">

*    * **Phase-separated solvent**: A membrane solvated with two solvent volumes containing different salt concentrations

<img src="figures/Solvation5_PhaseSeparatedSolvent.png" width="500">

* **Flooding**
	* **From the library**
	* **Import a solute molecule**: A membrane protein flooded with imported sucrose molecules

<img src="figures/Flooding3_ImportMultipleSolutes.png" width="500">

*	* **Import multiple solutes**
*	* **Import a lipid**
* **Custom unit cell**
* **COBY logo**

<img src="figures/COBY_Logo.png" width="500">

## **3. Manuscript tutorial**

* **Neuronal membrane**: Complex asymmetric membrane containing 58 lipid types.
* **Monolayers**: Two monolayers with solvent bridging the headgroups and vacuum bridging the tails.
* **Multiple solvent spaces**: Two voltage-gated potassium channels with modelled ion gradients.
* **Phase-separated membrane**: A phase-separated giant unilamellar vesicle membrane featuring liquid-ordered and liquid-disordered patches.
* **SARS-CoV-2 spike protein extended intermediate**: A viral and a host membrane bridged by an extended intermediate of the SARS-CoV-2 spike protein.
* **Gram-negative bacteria outer membrane and periplasm**: An example of a crowded system, featuring an asymmetric bilayer, four protein types, and four different solutes.
* **Multilamellar system with lipid-DNA complexes**: A type of a system used as an artificial transfection agent, featuring four stacked bilayers and DNA molecules in the inter-membrane spaces. 
* **Nanodic**
* **SIRAH force field**: Implementation of SIRAH molecules in COBY. 
* **Protein with benzene**: A benzene flooding setup surrounding a human K-Ras protein. 
* **Ionic liquid**: A solvent consisting of four types of cations and one type of anions in a liquid phase. 

