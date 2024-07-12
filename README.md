# Coarse-grained System Builder (COBY)

The program as well as documentation and a tutorial is avialable on the [GitHub](https://github.com/MikkelDA/COBY).

COBY is a Python-based software for building flat membranes in coarse-grained resolution. It handles asymmetric membranes, phase-separated systems, multiple bilayers, protein insertion, solvation and flooding with one or multiple solute molecules of choice.

COBY is versatile, fast, easy-to-use and it also allows for high-level of customisation. It can be used either as a **Python package** or directly from the **terminal command line**. It neatly handles multiple parameter libraries (even within the same system), making it developer-friendly. 

COBY is continuously under development and we welcome suggestions for new features. 

![](https://github.com/MikkelDA/COBY/raw/master/figures/COBY_Logo.png)

## Citation

Andreasen _et al._ (2024) TBD.

## Installation 

**Using pip (requires python>=3.9)**

    conda create --name COBY python==3.9 ipykernel

    conda activate COBY

    pip install COBY

    python -m ipykernel install --user --name=COBY

Ipykernel allows the user to use COBY environment within the Jupyter notebook.

## Basic usage 

For a detailed introduction to the software, please consult [COBY Documentation](https://github.com/MikkelDA/COBY/blob/master/COBY_Documentation.pdf).

For a quick reminder of available commands, you can use the [Cheat Sheet](https://github.com/MikkelDA/COBY/blob/master/COBY_CHEAT_SHEET.pdf).

Three [tutorials](https://github.com/MikkelDA/COBY/tree/master/Tutorial) are available as Jupyter Notebooks, one covering the basics, another covering more advanced functionalities and the final covering the systems shown in our soon-to-be released paper.

The tutorials includes (amongst other systems):

* **Simple membrane with protein**: A simple POPC membrane with a protein in solvent (water + 0.15 M NaCl) which explains the various arguments in COBY
* **Asymmetric membrane**: An asymmetric complex membrane (POPC, POPE and CHOL) in solvent with different area per lipid values for each membrane
* **Phase separated Membrane**: A phase separated membrane
* **Monolayers**: Two monolayers with solvent between them and vacuum over the pbc
* **Nanodisc**: A nanodisc with DMPC lipid contained within and with solvent surrounding it
* **Holes**: A membrane with multiple manually defined holes
* **Patches**: Multiple manually shaped membrane patches
* **Matryoshka membrane**: A matryoshka membrane with a protein in the center
* **Stacked membranes**: Three vertically stacked membranes
* **Mixed solvent**: A symmetric membrane solvated with regular and small water beads in specified ratios
* **Phase separated solvent**: A symmetric membrane solvated with two solvent volumes containing different salt concentrations
* **Flooding of imported solutes**: A membrane system that has been flooded with imported solute molecules, followed by solvation
* **COBY logo**: Our logo

## Licence

COBY is preserved under the [Apache License 2.0](https://github.com/MikkelDA/COBY/blob/main/LICENSE).

![](https://github.com/MikkelDA/COBY/raw/master/figures/membrane_protein.png)

