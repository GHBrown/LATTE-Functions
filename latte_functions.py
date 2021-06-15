
import copy
import numpy as np
import scipy as sp
#import matplotlib.pyplot as plt

"""
A collection of functions useful for working with the LATTE code (written by scientists at Los Alamos National Lab).
This collection of functions was written by Gabriel Brown, University of Illinois Urbana-Champaign.
"""


def getLATTEProperty(prprty,fileName):
    """
    Searches LATTE output text file and returns the specified
    property as a number.
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
    for line in lines:
        if (line[0:propertyNameLength]==prprty):
            propertyString=line.split()[-1] #get last element in split (the value)
            return float(propertyString)
    return 'Specified property not found.'


def getNumberOfAtoms(simulationDictionary):
    """
    Determines the number of atoms in a simulation specified by simulationDictionary.
    """
    groups=simulationDictionary['groupDictionaries']
    numberOfAtoms=0
    for group in groups:
        numberOfAtoms+=group["qArray"].shape[0]
    return numberOfAtoms
        

def LATTEDatFile2Dict(fileName):
    #WILL NEED SOME FIXING TO ACCOMODATE CHANGE TO SIMULATION DICTIONARIES
    """
    Takes arbitrary LATTE .dat file, creates array of atomic coordinates
    for every atom type present in file and stores a list of unique atom types.
    Also stores periodicity vectors of simulation domain.
    ---Inputs---
    fileName: file name (or relative path to file), string
    ---Outputs---
    qDictionary: dictionary containing the list of coordinate arrays (one for
        each atom type), a list of periodicity vectors, and a (string) list of
        the unique atom types in the same order, dictionary
    """
    fileObject=open(fileName,'r')
    content=fileObject.readlines()
    fileObject.close()
    content=[x.strip() for x in content]
    N_atoms=int(content[0])
    existingAtoms=[] #allocate list for atom types which have already been found
    qsList=[] #allocate list to contain coordinate lists for corresponding atoms
    boxPeriodicityVectors=[0]*3 #allocate list to store periodicity vectors of box size
    for irow,row in enumerate(content[1:4]): #extract box periodicity vectors
        rowList=row.split()
        curPeriodicityVector=np.array([float(elem) for elem in rowList])
        boxPeriodicityVectors[irow]=curPeriodicityVector
    for irow,row in enumerate(content[4:4+N_atoms]):
        rowList=row.split()
        curAtomType=rowList[0] #get atom type of current row
        curCoordinatesFloat=[float(elem) for elem in rowList[1:4]] #convert coordinate string list into list of floats
        if (curAtomType in existingAtoms): #if current atom type exists
            atomTypeIndex=existingAtoms.index(curAtomType) #find index of coordinate list
            qsList[atomTypeIndex].append(curCoordinatesFloat) #append current coordinates to this existing list
        else: #if current atom type does not exist yet
            existingAtoms.append(curAtomType) #add current atom type to existing atoms types
            qsList.append([curCoordinatesFloat]) #add first line in a new coordinate list for the new atom type
                
    qArrayList=[np.array(elem) for elem in qsList] #convert list of coordinate lists to list of numpy coordinate arrays
    qDictionary= {
        "boxVectors": boxPeriodicityVectors,
        "atomTypes": existingAtoms,
        "qs": qArrayList
    }
    return qDictionary


def groupDictionaryToInputBlock(fileName,groupDictionary,twoDimensional=False):
    #WILL NEED SOME FIXING TO ACCOMODATE CHANGE TO SIMULATION DICTIONARIES
    """
    Takes a dictionary defining a set of atomic coordinates and
    simulation domain and writes to LATTE's .dat format.
    ---Inputs---
    fileName: the file to which the configuration will be written (including the file extension), string
    groupDictionary: the dictionary which holds the information about the configuration, dictionary
    twoDimensional: optional input which allows users to configure systems which are traditionally two dimensionsal (only two box vectors) and set a third box vector which will take the value of this variable, float (if used)
    ---Outputs---
    NONE, the information is just written to file defined by fileName
    """

    fileObject=open(fileName,'w+')
    if twoDimensional: #take special care to set z values of two dimensional
        zHeight=twoDimensional #change to a more sensible variable name
        groupDictionary=translateCoordinates(groupDictionary,np.array([0,0,zHeight]))
        zPeriodicityVector=np.array([0,0,zHeight*2]) #make artificially large lattice vector in z direction to ensure no 2-D periodicity
        groupDictionary['boxVectors']=np.vstack((groupDictionary['boxVectors'],zPeriodicityVector))
    qList=groupDictionary['qs']
    N_atoms=getNumberOfAtoms(groupDictionary)
    N_dim=qList[0].shape[1]

    #write number of atoms
    fileObject.write('       '+str(N_atoms)+'\n')
    #loop over rows of cellSizeArray and write in cell dimensions
    for i in range(N_dim):
        dataFloat=groupDictionary['boxVectors'][i]
        dataString=[f"{element:.4f}" for element in dataFloat]
        fileObject.write('   '+'   '.join(dataString)+'\n')
    #loop over each atom type, loop over each atom in each type and write coordinates to file
    for q in qList:  
        for i in range(q.shape[0]):
            dataFloat=q[i,:]
            dataString=[f"{element:.5f}" for element in dataFloat]
            fileObject.write('C'+'    '+'   '.join(dataString)+'\n')

    print('Wrote .dat file:',fileName)
    print('Number of atoms:',N_atoms)
    fileObject.close()


"""
SIMULATION SHOULD BE A COLLECTION OF GROUPS, EACH GROUP HAS ITS OWN DATA, AND ONLY CONSISTS OF ONE ATOM TYPE WHICH IS LISTED IN GROUP DATA
rewrite code so this is satisfied
"""

def refreshMetaData(simulationDictionary):
    """
    Updates the metadata (list containing number of atomtypes, and
    list containing corresponding masses) of a given simulation dictionary.
    """
    atomTypeListFull=[]
    massListFull=[]
    for groupNumber,group in enumerate(simulationDictionary["groupDictionaries"]):
        atomTypeListFull+=group["atomTypeList"] #assumes atom types are sorted
        massListFull+=group["massList"] 

    simulationDictionary["nAtoms"]=getNumberOfAtoms(simulationDictionary)
    simulationDictionary["atomTypeList"]=atomTypeListFull
    simulationDictionary["massList"]=massListFull
    simulationDictionary["nAtomTypes"]=len(atomTypeListFull)

def simulationDictionaryToLAMMPSData(fileName,groupDictionary,twoDimensional=False):
    """
    Takes a dictionary defining a set of atomic coordinates and
    simulation domain and writes to LATTE's .dat format.
    ---Inputs---
    fileName: the file to which the configuration will be written (including the file extension), string
    simulationDictionary: the dictionary which holds the information about the configuration, dictionary
    twoDimensional: optional input which allows users to configure systems which are traditionally two dimensionsal (only two box vectors) and set a third box vector which will take the value of this variable, float (if used)
    ---Outputs---
    NONE, the information is written to file defined by fileName
    """

    fileObject=open(fileName,'w+')
    if twoDimensional: #take special care to set z values of two dimensional
        zHeight=twoDimensional
        simulationDictionary=translateCoordinates(groupDictionary,np.array([0,0,zHeight]))
        zPeriodicityVector=np.array([0,0,zHeight*2]) #make artificially large lattice vector in z direction to ensure no 2-D periodicity
        simulationDictionary["boxVectors"]=np.vstack((groupDictionary["boxVectors"],zPeriodicityVector))
    groups=simulationDictionary["groupDictionaries"]
    N_atoms=simulationDictionary["nAtoms"]
    N_atomtypes=simulationDictionary["nAtomTypes"]
    N_dim=3 #can one do better than this?

    #header and cumulative data
    fileObject.write('LAMMPS Description\n \n')
    fileObject.write(str(N_atoms)+' atoms\n \n')
    fileObject.write(str(N_atomtypes)+' atom types\n \n')
    #box size
    dimLabels=['x','y','z']
    for i in range(N_dim):
        boxVector=simulationDictionary['boxVectors'][i]
        minMax=[str(np.min(boxVector)),str(np.max(boxVector))]
        extentLabels=[dimLabels[i]+'lo',dimLabels[i]+'hi\n']
        fileObject.write(' '.join(minMax)+' '+' '.join(extentLabels))
    fileObject.write('\n')
    fileObject.write('Masses\n\n')
    for atomTypeNumber,mass in enumerate(simulationDictionary["massList"]):
        fileObject.write(str(atomTypeNumber+1)+' '+str(mass)+'\n')
    fileObject.write('\n')
    #information for each atom
    fileObject.write('Atoms\n\n')
    for igroup,group in enumerate(groups):  
        atomicData=group['qArray']
        for i in range(atomicData.shape[0]):
            dataFloatMeta=atomicData[i,0:4]
            dataFloatCoordinates=atomicData[i,4:8]
            dataString=[str(int(element)) for element in dataFloatMeta]
            dataString+=[f"{element:.5f}" for element in dataFloatCoordinates]
            fileObject.write(' '.join(dataString)+'\n')

    print('Wrote .lmp file:',fileName)
    print('Number of atoms:',N_atoms)
    fileObject.close()


def translateCoordinates(simulationDictionary,displacementVector,groupControl=None):
    """
    Given a configuration of atoms, edits the groups so that all* coordinates have been
    shift by a given vector. No other changes are performed.
    *certain group may be targeted by specifying groupControl
    ---Inputs---
    groupDictionary: a group of atoms, dictionary (one field of dictionary must be 'q')
    displacementVector: vector by which all coordinates of old group should be displaced, 1D numpy array
    ---Outputs---
    groupDictionaryTranlated: a new group with translated coordinates, dictionary
    """
    if (not groupControl):
        groupControl=[1]*len(simulationDictionary["groupDictionaries"]) #default to translating all groups

    simulationDictionaryTranslated=copy.deepcopy(simulationDictionary)
    for groupNumber,group in enumerate(simulationDictionary['groupDictionaries']):
        if groupControl[groupNumber]:
            q=group["qArray"]
            qCoordinates=q[:,4:8] #extract only coordinate part of data array
            displacementVectorTuple=tuple([displacementVector]*qCoordinates.shape[0])
            displacementArray=np.vstack(displacementVectorTuple)
            qprime=qCoordinates+displacementArray
            simulationDictionaryTranslated["groupDictionaries"][groupNumber]["qArray"][:,4:8]=qprime #overwrite with translated coordinates
    return simulationDictionaryTranslated


def returnToBox(qDictionary):
    """
    Moves atoms which are outside of the box back inside the proper region using
    the periodicity vectors.
    NOTE: This (probably) only works for a rectangular simulation domain which has all 90 degree angles.
    NOTE: This function is not recursive in that if an atom is more than one box length outside of the box, it will only be moved closer to the box by one box length rather than the necessary integer > 2
    ---Inputs---
    qDictionary: dictionary containing the list of coordinate arrays (one for
        each atom type), a list of periodicity vectors, and a (string) list of
        the unique atom types in the same order, dictionary
    ---Outputs---
    NONE, the input dictionary itself is automatically altered
    """
    tol=1e-3 #used to allow atoms very close to zero to remain there
    #this is where the orthrhombic assumption happens, only the diagonal elements of the boxVector "array" are taken
    boxSizeVector=np.zeros(3)
    for iVector,vector in enumerate(qDictionary["boxVectors"]):
        boxSizeVector[iVector]=vector[iVector]

    for iAtom,curCoordinates in enumerate(qDictionary["qs"]): #loop over coordinate arrays (one per atom type)
        for iDimension in range(curCoordinates.shape[1]):
            curCoordinatesCurDim=curCoordinates[:,iDimension] #x,y, or z coordinates for the current atom type
            tooSmallMask=curCoordinatesCurDim < 0-tol
            tooLargeMask=curCoordinatesCurDim >= boxSizeVector[iDimension]
            curCoordinates[:,iDimension]=curCoordinatesCurDim+boxSizeVector[iDimension]*(tooSmallMask*np.ones(curCoordinatesCurDim.shape)-tooLargeMask*np.ones(curCoordinatesCurDim.shape))#overwrite existing coordinates with inbound coordinates
        qDictionary["qs"][iAtom]=curCoordinates #overwrite full coordinate array with fully in bounds coordinate array
            

def computeAverageDisplacement(qDictionary1,qDictionary2):
    """
    Computes the average distance between corresponding atoms in two arrangements (which both must have the same number of atoms for every atom type).
    ---Inputs---
    qDictionary1/2: dictionaries containing the list of coordinate arrays (one for
        each atom type), a list of periodicity vectors, and a (string) list of
        the unique atom types in the same order, dictionary
    ---Outputs---
    meanDisplacement: 2 element list of the maximum and minimum atomic coordinates in any
        direction (meaning each element of list is scalar), list
    """
    coordinateList1=qDictionary1["qs"]
    coordinateList2=qDictionary2["qs"]
    numerator=0 #numerator of the average displacement
    nTotal=0 #accumulator for total number of atoms
    for atomType in range(len(coordinateList1)):
        curCoordinates1=coordinateList1[atomType]
        curCoordinates2=coordinateList2[atomType]
        nCur=curCoordinates1.shape[0] #number of atoms of current atom type
        curVectorDisplacement=curCoordinates2-curCoordinates1
        displacementVector=np.zeros(curVectorDisplacement.shape[0])
        for curDim in range(curVectorDisplacement.shape[1]):
            displacementVector+=np.power(curVectorDisplacement[:,curDim],2)
        displacementVector=np.power(displacementVector,1/2)
        atomTypeTotalDisplacement=np.sum(displacementVector)
        nTotal+=nCur
        numerator+=atomTypeTotalDisplacement
    meanDisplacement=numerator/nTotal
    return meanDisplacement
                          
    
def getCoordinateBounds(qDictionary):
    """
    Gets the extreme coordinates of the atomic assembly and saves as list for
    specifying range when plotting.
    ---Inputs---
    qDictionary: dictionary containing the list of coordinate arrays (one for
        each atom type), a list of periodicity vectors, and a (string) list of
        the unique atom types in the same order, dictionary
    ---Outputs---
    bounds: 2 element list of the maximum and minimum atomic coordinates in any
        direction (meaning each element of list is scalar), list
    """
    for i_q,q in enumerate(qDictionary["qs"]):
        qMin=np.ndarray.min(q)
        qMax=np.ndarray.max(q)
        if i_q==0:
            runningMin=qMin
            runningMax=qMax
        else:
            if qMin<runningMin:
                runningMin=qMin
            if qMax>runningMax:
                runningMax=qMax
    bounds=[runningMin,runningMax]
    return bounds

    
def makeSKF(outputFileName,parameterizationDictionary):
    """
    Writes a .skf format corresponding to the given parameterization dictionary.
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
    #Generate lines
    lineListList.append([pD["gridDist"],pD["nGridPoints"]]) #line 1
    lineListList.append(pD["EVec"]+[pD["SPE"]]+pD["UVec"]+pD["fVec"])
    lineListList.append([pD["mass"]]+pD["cVec"]+[pD["domain"][1]]+dVec)
    for r in gridPointsVec:
        eD=pD["elementFunction"](r) #elementDict, for brevity
        tempLine=[eD["Hdd0"],eD["Hdd1"],eD["Hdd2"],eD["Hpd0"],eD["Hpd1"],
                  eD["Hpp0"],eD["Hpp1"],eD["Hsd0"],eD["Hsp0"],eD["Hss0"],
                  eD["Sdd0"],eD["Sdd1"],eD["Sdd2"],eD["Spd0"],eD["Spd1"],
                  eD["Spp0"],eD["Spp1"],eD["Ssd0"],eD["Ssp0"],eD["Sss0"]]
        lineListList.append(tempLine)

    #Convert numbers to strings and write lines
    for lineList in lineListList:
        lineListString=[str(elem) for elem in lineList]
        fileObject.write(' '.join(lineListString)+'\n')
    fileObject.close()
    
