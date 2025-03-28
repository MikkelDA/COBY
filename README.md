# Coarse-grained System Builder (COBY)

The paper is currently available as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.23.604601v1).

The source code, documentation, and tutorials are avialable on [GitHub](https://github.com/MikkelDA/COBY).

COBY is a Python-based software tool for building flat membranes at a coarse-grained resolution. It supports asymmetric membranes, phase-separated systems, multiple bilayers, protein insertion, solvation, and flooding with one or multiple solute molecules of choice.

COBY is versatile, fast, easy to use, offering a high degree of customisation. It can be used either as a **Python package** or directly from the **terminal command line**. It efficiently manages multiple parameter libraries (even within the same system), making it developer-friendly. 

COBY is continuously under development, and we welcome suggestions for new features, which should be submitted under [Issues](https://github.com/MikkelDA/COBY/issues). 

![](https://github.com/MikkelDA/COBY/raw/master/figures/COBY_Logo.png)

## Citation (bioRxiv)
```
@article {Andreasen2024COBY,
	author = {Andreasen, Mikkel D. and Souza, Paulo C. T. and Schi{\o}tt, Birgit and Zuzic, Lorena},
	title = {Creating Coarse-Grained Systems with COBY: Towards Higher Accuracy in Membrane Complexity},
	elocation-id = {2024.07.23.604601},
	year = {2024},
	doi = {10.1101/2024.07.23.604601},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/07/23/2024.07.23.604601},
	eprint = {https://www.biorxiv.org/content/early/2024/07/23/2024.07.23.604601.full.pdf},
	journal = {bioRxiv}
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

## Licence

COBY is preserved under the [Apache License 2.0](https://github.com/MikkelDA/COBY/blob/main/LICENSE).

![](https://github.com/MikkelDA/COBY/raw/master/figures/membrane_protein.png)
