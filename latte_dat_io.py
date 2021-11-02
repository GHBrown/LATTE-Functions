# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:42:04 2021

@author: danpa
"""

import numpy as np 
import os
import re
from ase import Atoms

def write_latte_dat(ase_obj,filename):
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
    with open(filename,'w+') as f:
        f.write("      ")
        f.write('%d\n'%(natom))
        f.write(rx_+" \n")
        f.write(ry_+" \n")
        f.write(rz_+" \n")
        
        for a in ase_obj:
            f.write(a.symbol+" ")
            pos=np.array(a.position)
            str_pos=" ".join(map(str,pos))
            f.write(str_pos+" \n")
            
def read_latte_dat(filename):
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
            
        
        atom_obj=Atoms(symbols,positions=pos,\
                       cell=np.array([cell_x,cell_y,cell_z]))
        
    return atom_obj

if __name__=="__main__":
    import shift
    filename="test.dat"
    ase_obj=shift.make_graphene(stacking=['A','B'],cell_type='hex',n_layer=2,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=3.0,sym=['C','C'])
    write_latte_dat(ase_obj,filename)
    ase_obj=read_latte_dat(filename)
    print(ase_obj)
    print(ase_obj.positions.shape)