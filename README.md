# Coarse-grained System Builder (COBY)

The publication is available open-source at [https://doi.org/10.1021/acs.jcim.5c00069](https://doi.org/10.1021/acs.jcim.5c00069).

The source code, documentation, and tutorials are avialable on [GitHub](https://github.com/MikkelDA/COBY).

COBY is a Python-based software tool for building flat membranes at a coarse-grained resolution. It supports asymmetric membranes, phase-separated systems, multiple bilayers, protein insertion, solvation, and flooding with one or multiple solute molecules of choice.

COBY is versatile, fast, easy to use and offers a high degree of customisation. It can be used either as a **Python package** or directly from the **terminal command line**. It efficiently manages multiple parameter libraries (even within the same system), making it developer-friendly. 

COBY is continuously under development, and we welcome suggestions for new features, which should be submitted under [Issues](https://github.com/MikkelDA/COBY/issues). 

![](https://github.com/MikkelDA/COBY/raw/master/figures/COBY_Logo.png)

## Citation (JCIM)
```
@article{Andreasen2025,
	author = {Andreasen, Mikkel D. and Souza, Paulo C. T. and Schiøtt, Birgit and Zuzic, Lorena},
	title = {Creating Coarse-Grained Systems with COBY: Toward Higher Accuracy of Complex Biological Systems},
	journal = {Journal of Chemical Information and Modeling},
	volume = {0},
	number = {0},
	pages = {null},
	year = {0},
	doi = {10.1021/acs.jcim.5c00069},
	URL = {https://doi.org/10.1021/acs.jcim.5c00069},
	eprint = {https://doi.org/10.1021/acs.jcim.5c00069}
}
```

## Installation 

**Using pip (requires python>=3.9)**

    conda create --name COBY python==3.9 ipykernel

    conda activate COBY

    pip install COBY

    python -m ipykernel install --user --name=COBY

Ipykernel allows the user to use COBY environment within the Jupyter notebook.

## OS compatability

COBY has been developed and thoroughly tested on Ubuntu Linux. It has also been tested to a lesser extent on Windows 11 and macOS, with all tutorial systems running without issue.

## Testing

Installation can be tested using the tutorial notebooks, which cover the majority of the functions available in COBY.

## Basic usage 

For a detailed introduction to the software, please refer to the [COBY Documentation](https://github.com/MikkelDA/COBY/blob/master/COBY_Documentation.pdf).

For a quick reminder of available commands, you can use the [Cheat Sheet](https://github.com/MikkelDA/COBY/blob/master/COBY_CHEAT_SHEET.pdf).

Three [tutorials](https://github.com/MikkelDA/COBY/tree/master/Tutorial) are available as Jupyter Notebooks: one covering the basics, another focusing on more advanced functionalities, and the final tutorial showcasing the systems from the manuscript.

## Known issues

* Recursion depth crash on Mac
	* Problem: It has been reported that running COBY with newer versions of Shapely on a Mac may lead to crashes caused by exceeding the recursion depth limit. The reason for why or how this happens is unclear but is likely not something that can be fixed by the COBY developer. The error may be relatod the the macOS version as the problem is known to occur on Sonoma but not on Sequoia.
	* Solution 1: Downgrade Shapely to version 2.0.2, which is the newest version known to be stable on Mac.
	* Solution 2: Update your Mac to Sequoia or an even newer version.

## Licence

COBY is preserved under the [Apache License 2.0](https://github.com/MikkelDA/COBY/blob/main/LICENSE).

![](https://github.com/MikkelDA/COBY/raw/master/figures/membrane_protein.png)
