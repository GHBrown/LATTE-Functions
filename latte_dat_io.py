# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:42:04 2021

@author: danpa
"""

import numpy as np 
import os
import re
from ase import Atoms
import ase.io
from ase.calculators.lammps import Prism, convert
from ase.utils import reader, writer
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, 'C:/Users/danpa/Documents/research/latte_tools/flat-graphene/flatgraphene')
import shift

def write_latte_dat(ase_obj,filename,electron_file=None):
    """write latte coordinate file from ase.atom.Atom object
    
    :param ase_obj: (ase.atom.Atom obj) ase atoms object where geometry is stored
    
    :param filename: (str) filename to write data file to
    """
    cell=np.array(ase_obj.get_cell())
    rx_=" ".join(map(str,cell[0,:]))
    ry_=" ".join(map(str,cell[1,:]))
    rz_=" ".join(map(str,cell[2,:]))
    
    xyz=ase_obj.get_positions()
    natom=np.shape(xyz)[0]
    
    #compare mass in ase object to mass in electrons file to find correct element
    if electron_file!=None:
        with open(electron_file) as f:
            lines=f.readlines()
            mass_dict={}
            for i,l in enumerate(lines):
                if i<2:
                    continue
                properties=l.split(" ")
                temp_dict={float(properties[7]) : properties[0]} # {Mass : Element}
                mass_dict.update(temp_dict)
                
    with open(filename,'w+') as f:
        f.write("      ")
        f.write('%d\n'%(natom))
        f.write(rx_+" \n")
        f.write(ry_+" \n")
        f.write(rz_+" \n")
        
        for a in ase_obj:
            if electron_file!=None:
                get_mass=a.mass
                for key in mass_dict.keys():
                    if np.isclose(float(key),get_mass,rtol=0.00001):
                        symbol=mass_dict[key]
            else:
                symbol=a.symbol
                
            f.write(symbol+" ")
            pos=np.array(a.position)
            str_pos=" ".join(map(str,pos))
            f.write(str_pos+" \n")
            
def read_latte_dat(filename,electron_file=None):
    """read latte data file into ase.Atoms object
    
    :param filename: (str) filename of latte data file to read
    
    :returns: (ase.Atoms) ase.Atoms object containing chemical symbols, positions
              and cell of system"""
              
    with open(filename,"r") as f:
        lines=f.readlines()
        
        natom=int(re.findall(r'[-+]?[.]?[:\.\d]+',lines[0])[0])
        pos=np.zeros((natom,3))
        symbols=np.empty(natom,dtype=np.unicode_)
        cell_x=re.findall(r'[-+]?[.]?[:\.\d]+',lines[1])
        cell_y=re.findall(r'[-+]?[.]?[:\.\d]+',lines[2])
        cell_z=re.findall(r'[-+]?[.]?[:\.\d]+',lines[3])
        
        for i,l in enumerate(lines[4:]):
            pos[i,:]=re.findall(r'[-+]?[.]?[:\.\d]+',l)
            sym=l.split(" ")[0]
            symbols[i]=sym
            
        #include masses in object
        if electron_file!=None:
            with open(electron_file) as f:
                lines=f.readlines()
                mass_dict={}
                for i,l in enumerate(lines):
                    if i<2:
                        continue
                    properties=l.split(" ")
                    temp_dict={properties[0] : float(properties[7])} # {Mass : Element}
                    mass_dict.update(temp_dict)
                    
            masses=np.zeros(natom)
            for k in mass_dict.keys():
                ind=np.where(symbols==k)
                masses[ind]=mass_dict[k]
            
            atom_obj=Atoms(positions=pos,\
                       cell=np.array([cell_x,cell_y,cell_z]))
            atom_obj.set_masses(masses)
        else:
            try:
                atom_obj=Atoms(symbols,positions=pos,\
                           cell=np.array([cell_x,cell_y,cell_z]))
            except:
                print("atomic labels in .dat file may not be atomic symbols. Try passing associated electrons.dat file")
        
    return atom_obj

    
if __name__=="__main__":
    import shift
    filename="test.dat"
    a_nn=2.529/np.sqrt(3)
    ase_obj=shift.make_graphene(stacking=['A','B'],cell_type='hex',n_layer=2,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=a_nn,sep=3.35,sym=['B','Ti']\
                    ,mass=[12.01,12.02],h_vac=3.0)
    efile="C:/Users/danpa/Documents/research/twisted-graphene-geometry-optimization/parameters_potentials/latte/Porezag_Popov_Van_Alsenoy/latte/electrons.dat"
    write_latte_dat(ase_obj,"test_coords.dat",electron_file=efile)
    obj=read_latte_dat("test_coords.dat",electron_file=efile)
    print(obj.get_masses())
    
    