import numpy as np 
import pandas as pd 
from ReadFile import Read
import astropy.units as u
from tabulate import tabulate

def ComponentMass(filename, ptype):
    """
    Calculate the total mass of a specific component of a galaxy.

    Imputs:
    filename (str): The path to the file containing the galaxy data.
    ptype (int): The particle type for which the mass is to be calculated.
                The particle types are: Halo type (1), Disk type (2), Bulge type (3).

    Returns:
    Final_mass (str): The total mass of the specified component formatted in scientific notation with units of solar masses.
    """
    #Read in the file
    time, total, data = Read(filename)
    index = np.where(data["type"] == ptype) #finding all the indexs in the data corresponding to particle type inputted.
    mass = data["m"][index]/100 #calculating the mass of the particles in 10^12 solar masses.
    total_mass = np.round(np.sum(mass), 3)  # total masses of the particles to get the total mass of the galaxy component.
    final_mass = float(total_mass)
    return final_mass

