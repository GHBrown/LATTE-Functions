
import copy
from datetime import datetime
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

"""
A collection of functions useful when performing tight binding simulations
with the LATTE-LAMMPS research codes.
"""

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
    fileObject=open(fileName,'r')
    propertyNameLength=len(prprty)
    lines=fileObject.readlines()
    fileObject.close()
    #if statements for tensor/vector properties go here
    for line in lines:
        if (line[0:propertyNameLength]==prprty):
            propertyString=line.split()[-1] #get last element in split (the value)
            return float(propertyString)
    return 'Specified property not found.'
    
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
    pD=parameterizationDictionary #for brevity
    fileObject=open(outputFileName,'w')
    dVec=[0]*10 #just nonsense values since ds are placeholders
    gridPointsVec=pD["gridDist"]*np.array(range(pD["nGridPoints"]))
    lineListList=[]

    #grid and atomic information
    lineListList.append([pD["gridDist"],pD["nGridPoints"]]) #line 1
    if (pD["type"]=="homonuclear"):
        lineListList.append(pD["EVec"]+[pD["SPE"]]+pD["UVec"]+pD["fVec"]) #homonuclear line 2
        lineListList.append([pD["mass"]]+pD["cVec"]+[pD["domainTB"][1]]+dVec) #homonuclear line 3
    elif (pD["type"]=="heteronuclear"):
        lineListList.append([12345]+pD["cVec"]+[pD["domainTB"][1]]+dVec) #heteronuclear line 2, mass (first entry) is placeholder so use obvious dummy value
    else:
        print("ERROR: parameterization type must be \"heteronuclear\"")
        print("       or \"homonuclear\"")

    #integral table
    for r in gridPointsVec:
        if (r < pD["domainTB"][0]): #write 1.0s for r < r_min (what LATTE expects)
            tempLine=[1.0]*20
        else:
            eD=pD["elementFunction"](r) #elementDict, for brevity
            tempLine=[eD["Hdd0"],eD["Hdd1"],eD["Hdd2"],eD["Hpd0"],eD["Hpd1"],
                    eD["Hpp0"],eD["Hpp1"],eD["Hsd0"],eD["Hsp0"],eD["Hss0"],
                    eD["Sdd0"],eD["Sdd1"],eD["Sdd2"],eD["Spd0"],eD["Spd1"],
                    eD["Spp0"],eD["Spp1"],eD["Ssd0"],eD["Ssp0"],eD["Sss0"]]
        lineListList.append(tempLine)

    #spline
    lineListList.append(['Spline'])
    lineListList.append([1, pD["domainTB"][1]]) #nInt cutoff
    lineListList.append([0, 0, -1]) #a1 a2 a3  (for exp(-a1*r+a2)+a3)
    lineListList.append([pD["domainTB"][0],pD["domainTB"][1],0,0,0,0,0,0]) #start end c0 c1 c2 c3 c4 c5
    #(for c0+c1(r-r0)+c2(r-r0)^2+c3(r-r0)^3+c4(r-r0)^4+c5(r-r0)^5

    #convert numbers to strings and write lines
    for lineList in lineListList:
        lineListString=[str(elem) for elem in lineList]
        fileObject.write(' '.join(lineListString)+'\n')
    fileObject.close()
    

def plotSKF(fileName,domain):
    """
    Plots the elements of Hamiltonian and overlap matrix against distance r
    Only plots an element if lower 2/3 of respective .skf column has
    nonzero magnitude (since LATTE expects 1.0s below r_min).
    ---Inputs---
    fileName: name of a .skf file, string
    domain: domain (r values) on which to plot elements, list [r_min, r_max]
    ---Outputs---
    NONE: makes and shows plots
    """

    with open(fileName,'r') as f:
        linesAsterisksCommas=f.readlines()

    linesAsterisks=[line.replace(',','') for line in linesAsterisksCommas] #clean lines of commas
    lines=[]*len(linesAsterisks)
    for line in linesAsterisks: #write asterisk exapanded lines to lines variable
        if '*' in line:
            lineSplit=line.split() #split line on spaces
            for i_entry,entry in enumerate(lineSplit):
                if '*' in entry: #split *-containing entry on *
                    entrySplit=entry.split('*')
                    num=float(entrySplit[1]) #number to be repeated
                    timesRep=int(float(entrySplit[0])) #times to repeat
                    expandedEntries=' '.join([str(num)]*timesRep)
                    lineSplit[i_entry]=expandedEntries
            lines.append(' '.join(lineSplit))

        else:
            lines.append(line)
                
    if '@' in lines[0]:
        gridDist=float(lines[1].split()[0]) #distance between gridpoints, on second line for extended format
        integralTableLineLength=40 #.skf is in extended format
        HIntegralLabels=['Hff0', 'Hff1', 'Hff2', 'Hff3', 'Hdf0',
                         'Hdf1', 'Hdf2', 'Hdd0', 'Hdd1', 'Hdd2',
                         'Hpf0', 'Hpf1', 'Hpd0', 'Hpd1', 'Hpp0',
                         'Hpp1', 'Hsf0', 'Hsd0', 'Hsp0', 'Hss0']
        SIntegralLabels=['Sff0', 'Sff1', 'Sff2', 'Sff3', 'Sdf0',
                         'Sdf1', 'Sdf2', 'Sdd0', 'Sdd1', 'Sdd2',
                         'Spf0', 'Spf1', 'Spd0', 'Spd1', 'Spp0',
                         'Spp1', 'Ssf0', 'Ssd0', 'Ssp0', 'Sss0']
    else:
        gridDist=float(lines[0].split()[0]) #distance between gridpoints, on first line for simple format
        integralTableLineLength=20 #.skf is in simple format
        HIntegralLabels=['Hdd0', 'Hdd1', 'Hdd2', 'Hpd0', 'Hpd1',
                         'Hpp0', 'Hpp1', 'Hsd0', 'Hsp0', 'Hss0']
        SIntegralLabels=['Sdd0', 'Sdd1', 'Sdd2', 'Spd0', 'Spd1',
                         'Spp0', 'Spp1', 'Ssd0', 'Ssp0', 'Sss0']

    linesSplit=[line.split() for line in lines]
    #find index of integral table's first line
    firstLineIndex=False
    foundFirstIntegralLine=False
    i_line=0
    while (not foundFirstIntegralLine):
        splitLine=linesSplit[i_line]
        if len(splitLine)==integralTableLineLength:
            foundFirstIntegralLine=True
            firstLineIndex=i_line
            if (integralTableLineLength==20):
                #in simple format line BEFORE integral table has 20 entries
                #correct for this
                firstLineIndex+=1 
        i_line+=1
        
    #find index of integral table's last line
    afterLineIndex=False
    for i_line,splitLine in enumerate(linesSplit):
        if ('Spline' in splitLine):
            afterLineIndex=i_line #'Spline' occurs on line after integral table

    #make plucked integral table into numpy arrays for H and S
    integralTableLines=linesSplit[firstLineIndex:afterLineIndex]
    integralTable=np.array([[float(entry) for entry in splitLine] for splitLine in integralTableLines])
    HTable=integralTable[:,0:int(integralTableLineLength/2)] #first half of columns are for H
    STable=integralTable[:,int(integralTableLineLength/2):] #second half of columns are for S

    numPoints=HTable.shape[0] #number of points at which integrals are given
    rValues=gridDist*np.arange(numPoints)

    #plot nonzero H (Hamiltonian) integrals
    for i_integral,integralValues in enumerate(np.transpose(HTable)):
        if (sum(abs(integralValues[int(numPoints/3):]))>0.0):
               plt.plot(rValues,integralValues,label=HIntegralLabels[i_integral])
    plt.xlim(domain)
    plt.xlabel('internuclear distance, [Bohr radii]')
    plt.ylabel('energy, [Hartrees]')
    plt.legend()
    plt.show()

    #plot S (overlap) integrals
    for i_integral,integralValues in enumerate(np.transpose(STable)):
        if (sum(abs(integralValues[int(numPoints/3):]))>0.0):
               plt.plot(rValues,integralValues,label=SIntegralLabels[i_integral])
    plt.xlim(domain)
    plt.xlabel('internuclear distance, [Bohr radii]')
    plt.ylabel('energy, [Hartrees]')
    plt.legend()
    plt.show()


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
    #define unit conversion constants
    ang_per_bohr=0.529177 # [Anstroms/Bohr radius]
    eV_per_hart=27.2114 # [eV/Hartree]

    pD=parameterizationDictionary #for brevity
    lineList=[] #list to which complete strings will be appended 1 per line
    r_min=pD["domainPair"][0] #[Bohr radii]
    r_max=pD["domainPair"][1] #[Bohr radii]
    rValues=np.linspace(r_min,r_max,pD["nGridPoints"]) #[Bohr radii]

    lineList.append('# DATE: ' + datetime.today().strftime('%Y-%m-%d') + ' UNITS: metal CONTRIBUTOR: ' + pD["contributor"] + '\n')
    lineList.append('# ' + pD["pairDescription"] + '\n')
    lineList.append('\n')
    lineList.append(pD["pairKeyword"] + '\n')
    lineList.append('N ' + str(pD["nGridPoints"]) + '\n')
    lineList.append('\n')

    for i_r,r in enumerate(rValues): #r [Bohr radii]
        energy,force=pD["pairFunction"](r) #energy [Hartrees], force [Hartrees/Bohr radius]
        r_ang=r*ang_per_bohr
        energy_eV=energy*eV_per_hart
        force_eV_ang=force*(eV_per_hart/ang_per_bohr)
        tempValues=[str(i_r+1),str(r_ang),str(energy_eV),str(force_eV_ang)]
        lineList.append(' '.join(tempValues)+'\n')

    with open(outputFileName,'w') as f:
        f.writelines(lineList)
