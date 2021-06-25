
## LATTE-LAMMPS Functions

A collection of Python functions useful for working with the LATTE and LAMMPS research codes.

These functions are intended to:
- generate LAMMPS pairwise potential files
- generate tight binding parameterizations in the form of `.skf` files, which LATTE can convert to a usable parameterization with its built-in converter `DLtab.py`


## Working with the provided functions

The main object needed to work with these functions is the "parameterization dictionary", which contains information about species, and their tight binding and pairwise interactions.

**Note that all units in the dictionary are atomic** (Bohr radii, Hartrees, etc.)

This object is Python dictionary which has the following form:

```
homonuclearDictionary={ #defines interaction between atoms of same species
    "mass": mass of atomic species,
    "gridDist": spacing of grid on which TB interaction is computed, #[Bohr radii]
    "nGridPoints": number of grid point on which to compute interactions,
    "type":'homonuclear',
    "elementFunction": function pointer to function that takes distance and returns Hamiltonian and overlap matrix elements,
    "domainTB": two element list specifying range on which TB parameterization is viable, #domain of viability, [r_min, r_cut] [Bohr radii]
    "EVec":three element list of on-site energies, #[E_d, E_p, E_s] 
    "SPE": spin polarization error (unused),
    "UVec": three element list of Hubbard U values, #[U_d, U_p, U_s]
    "fVec": occupations of neutral atom in ground state #[f_d, f_p, f_s]
    "cVec": 8 element list of pairwise correction polynomial coefficient, #[c_2,c_3,...,c_9]
    "pairFunction": function pointer to function which takes distance and returns pairwise energy and force,
    "domainPair": two element list specifying range on which pairwise correction is viable, #domain of viability, [r_min, r_cut] [Bohr radii]
    "pairKeyword": pair keyword for LAMMPS potential,
    "pairDescription": string used at LAMMPS potential description,
    "contributor": appropriate author(s)
    }

heteronuclearDictionary={ #defines interaction between atoms of different species
    "gridDist": spacing of grid on which TB interaction is computed, #[Bohr radii]
    "nGridPoints": number of grid point on which to compute interactions,
    "type":'heteronuclear',
    "elementFunction": function pointer to function that takes distance and returns Hamiltonian and overlap matrix elements,
    "domainTB": two element list specifying range on which TB parameterization is viable, #domain of viability, [r_min, r_cut] [Bohr radii]
    "cVec": 8 element list of pairwise correction polynomial coefficient, #[c_2,c_3,...,c_9]
    }
```

The input-output format of the functions pointers assigned to `"elementFunction"` and `"pairFunction"` are:
```
def nameOfFunctionToComputeHandSElement(r):
    #do some r-dependent computations of matrix elemements
    elementDict={
        #0:sigma, 1:pi, 2:delta
        #Hamiltonian elements
        "Hss0": computed value,
        "Hsp0": computed value,
        "Hsd0": computed value,
        "Hpp0": computed value,
        "Hpp1": computed value,
        "Hpd0": computed value,
        "Hpd1": computed value,
        "Hdd0": computed value,
        "Hdd1": computed value,
        "Hdd2": computed value,
        #overlap matrix elements
        "Sss0": computed value,
        "Ssp0": computed value,
        "Ssd0": computed value,
        "Spp0": computed value,
        "Spp1": computed value,
        "Spd0": computed value,
        "Spd1": computed value,
        "Sdd0": computed value,
        "Sdd1": computed value,
        "Sdd2": computed value,
        }
    return elementDict

def nameOfFunctionToComputePairwiseEandF(r):
    #do some r-dependent computations of energy and force
    return energy, force
```


## Documentation

```
def getLATTEProperty(prprty,fileName):
    """
    Searches LATTE output text file and returns the specified
    property as a float.
    NOTE: does not yet work for nonscalar values like tensors,
        do to this one would need to have explicity checks for
        the names associated with tensor inputs
    ---Inputs---
    prprty: name of desired property (must be exact), string
    fileName: name of file in which to search
    ---Outputs---
    propertyValue: value(s) of property, data type depends on
        property (could be a scalar, tensor, etc.)
    """

def makeSKF(outputFileName,parameterizationDictionary):
    """
    Writes a .skf formatted file corresponding to the given parameterization
    dictionary
    ---Inputs---
    outputFileName: name of .skf file, string
    parameterizationDictionary: dictionary with all info about the parameterization
    ---Outputs---
    NONE: structured file of specified name is created
    """

def plotSKF(fileName,domain):
    """
    Plots the elements of Hamiltonian and overlap matrix against distance r
    ---Inputs---
    fileName: name of a .skf file, string
    domain: domain (r values) on which to plot elements, list [r_min, r_max]
    ---Outputs---
    NONE: makes and shows plots
    """

def makeLAMMPSPairwiseTable(outputFileName,parameterizationDictionary):
    """
    Writes a pairwise potential as a LAMMPS pairwise potential table.
    All of parameterization dictionary should be in atomic units (Bohr radii,
    Hartrees, etc.). This function generates a table in metal units (eV, Angstroms)
    and does the conversion internally
    ---Inputs---
    outputFileName: name of .table file
    parameterizationDictionary: dictionary with all info about the parameterization
    ---Outputs---
    NONE: structured file of specified name is created
    """
```

This documentation has been made as complete as possible, but please see these links for more information about the file formats and expected inputs:
- [`.skf` format](https://dftb.org/fileadmin/DFTB/public/misc/slakoformat.pdf)
- [LAMMPS pairwise table format](https://docs.lammps.org/pair_table.html)
