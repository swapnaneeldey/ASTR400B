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
    print(final_mass)

local_grp = []  # List to store the total mass of each galaxy in the local group
disk = []  # List to store the disk mass of each galaxy
buldge = []  # List to store the bulge mass of each galaxy
halo = []  # List to store the halo mass of each galaxy

# iterating through each of the three galaxies to fill the above lists.
for i in range(3):
    if i == 0:
        filename = "/Users/swapnaneeldey/Desktop/ASTR400Bfiles/MW_000.txt" # Path to the file containing the MW data
        name = "MW"
    elif i == 1:
        filename = "/Users/swapnaneeldey/Desktop/ASTR400Bfiles/M31_000.txt" # Path to the file containing the M31 data
        name = "M31"
    else:
        filename = "/Users/swapnaneeldey/Desktop/ASTR400Bfiles/M33_000.txt" # Path to the file containing the M33 data
        name = "M33"
    
    # Calculate the mass of the halo component
    halo_mass = ComponentMass(filename, 1)
    print(f"{name} Halo mass: {halo_mass}")
    
    # Calculate the mass of the disk component
    disk_mass = ComponentMass(filename, 2)
    print(f"{name} Disk mass: {disk_mass}")
    
    # Calculate the mass of the bulge component 
    if i == 2:
        bulge_mass = "0.000 x 10^12 solMass" #uing if statement since M33 does not have a bulge component
    else:
        bulge_mass = ComponentMass(filename, 3)   # Calculate the mass of the bulge component for the rest galaxies
    print(f"{name} Bulge mass: {bulge_mass}")
    
    # Append the masses to their respective lists
    disk.append(float(disk_mass[0:5]))
    buldge.append(float(bulge_mass[0:5]))
    halo.append(float(halo_mass[0:5]))
    
    # Calculate the total mass of the galaxy
    total_galaxy_mass = float(halo_mass[0:5]) + float(disk_mass[0:5]) + float(bulge_mass[0:5])
    print(f"Final total mass of {name}:", total_galaxy_mass, "x 10^12 solMass")
    
    # Append the total mass to the local group list
    local_grp.append(total_galaxy_mass)
    
    # Calculate the baryon fraction of the galaxy
    f_b = ((float(disk_mass[0:5]) + float(bulge_mass[0:5])) / total_galaxy_mass)
    print(f"baryon fraction of {name}:", np.round(f_b, 3))

# Calculate and print the total mass of the local group
print("Local Group mass:", np.sum(local_grp), "x 10^12 solMass")

# Calculate and print the baryon fraction of the local group
f_b_tot = ((np.sum(disk) + np.sum(buldge)) / np.sum(local_grp))
print("baryon fraction of Local Group:", np.round(f_b_tot, 3))