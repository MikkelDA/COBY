
from COBY.version import __version__, __version_changes__
from COBY.citation import __doi__, __citation__

### Prints the current version
class version():
    def __init__(self):
        print(__version__)

### Prints the version changes for the current version
class version_changes():
    def __init__(self):
        print(__version_changes__)
changes   = version_changes
changelog = version_changes

### Prints the doi for the paper (https://doi.org/10.1021/acs.jcim.5c00069)
class doi():
    def __init__(self):
        print(__doi__)
DOI = doi

### Prints the doi for the paper (https://doi.org/10.1021/acs.jcim.5c00069)
class citation():
    def __init__(self):
        print(__citation__)

