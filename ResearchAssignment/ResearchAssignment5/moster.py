import numpy as np

class AbundanceMatching:
    """ Class to define the abundance matching relations from 
    Moster et al. 2013, which relate the stellar mass of a galaxy
    to the expected dark matter halo mass, according to 
    Lambda Cold Dark Matter (LCDM) theory """
    
    
    def __init__(self, mhalo, z):
        """ Initialize the class
        
        PARAMETERS
        ----------
            mhalo: float
                Halo mass in Msun
            z: float
                redshift
        """
        
        #initializing the parameters:
        self.mhalo = np.array(mhalo, dtype=float)  # Convert to NumPy array
        self.z = float(z)  # Redshift
        
        
    def logM1(self):
        """eq. 11 of Moster 2013
        OUTPUT: 
            M1: float 
                characteristic mass in log(Msun)
        """
        M10      = 11.59
        M11      = 1.195 
        return M10 + M11*(self.z/(1+self.z))  
    
    
    def N(self):
        """eq. 12 of Moster 2013
        OUTPUT: 
            Normalization for eq. 2
        """
        N10      = 0.0351
        N11      = -0.0247
    
        return N10 + N11*(self.z/(1+self.z))
    
    
    def Beta(self):
        """eq. 13 of Moster 2013
        OUTPUT:  power of the low mass slope"""
        beta10      = 1.376
        beta11      = -0.826
    
        return beta10 + beta11*(self.z/(1+self.z))
    
    def Gamma(self):
        """eq. 14 of Moster 2013
        OUTPUT: power of the high mass slope """
        gamma10      = 0.608
        gamma11      = 0.329
    
        return gamma10 + gamma11*(self.z/(1+self.z))
    
    
    def SHMratio(self):
        """ 
        eq. 2 of Moster + 2013
        The ratio of the stellar mass to the halo mass
        
        OUTPUT: 
            SHMratio float
                Stellar mass to halo mass ratio
        """
        M1 = 10**self.logM1() # Converting characteristic mass 
        # to Msun from Log(Msun)
        
        A = (self.mhalo/M1)**(-self.Beta())  # Low mass end
        
        B = (self.mhalo/M1)**(self.Gamma())   # High mass end
        
        Norm = 2*self.N() # Normalization
    
        SHMratio = Norm*(A+B)**(-1)
    
        return SHMratio 
    
 # Q1: add a function to the class that takes the SHM ratio and returns 
# The stellar mass 
    def StellarMass(self):
        """Method to compute the stellar mass
        using eq 2 of Moster+2013 (stellar / halo mass ratio)
        
        Output;
            starMass: float, stellar mass in Msun"""
        
        starMass = self.mhalo*self.SHMratio()

        return starMass
    
    def error(self):
        """Method to compute the error in the stellar mass
        using eq 2 of Moster+2013 (stellar / halo mass ratio)
        
        Output;
            starMass: float, stellar mass in Msun"""
        M1_up = 10**(11.59+0.236) # Converting characteristic mass 
        # to Msun from Log(Msun)
        M1_low = 10**(11.59-0.236) # Converting characteristic mass
        A = (self.mhalo/M1_up)**(-(1.376+0.153))  # Low mass end
        
        B = (self.mhalo/M1_low)**(0.608+0.059)   # High mass end
        
        Norm = 2*(0.0351+0.0058) # Normalization
        print(Norm*(A+B)**(-1))
        sm_upperlimit = self.mhalo*Norm*(A+B)**(-1)

        M1 = 10**(11.59-0.236) # Converting characteristic mass 
        # to Msun from Log(Msun)
        
        A = (self.mhalo/M1_low)**(-(1.376-0.153))  # Low mass end
        
        B = (self.mhalo/M1_up)**(0.608-0.059)   # High mass end
        
        Norm = 2*(0.0351-0.0058) # Normalization
        print(Norm*(A+B)**(-1))
        sm_lower_lim = self.mhalo*Norm*(A+B)**(-1)

        return sm_upperlimit, sm_lower_lim

    def compute_sigma_m(self, unc_N = 0.0058, unc_M1 = 0.236, unc_beta = 0.153, unc_gamma = 0.059):
        M = self.mhalo                          # Halo mass
        logM1 = self.logM1()                    # Log base 10 of characteristic mass
        M1 = 10**logM1                          # Characteristic mass
        N = self.N()                            # Normalization factor
        beta = self.Beta()                      # Slope at low mass end
        gamma = self.Gamma()                    # Slope at high mass end

        x = M / M1
        denom = (x**(-beta) + x**gamma)
        
        # Partial derivatives
        dm_dN = (2 * M) / denom

        dm_dlogM1 = 2 * M * N * np.log(10) * (gamma * x**gamma - beta * x**(-beta)) / denom**2

        dm_dbeta = 2 * M * N * (x**(-beta) * np.log(x)) / denom**2

        dm_dgamma = -2 * M * N * (x**gamma * np.log(x)) / denom**2

        # Total variance from each parameter's uncertainty
        sigma_m_squared = (
            (dm_dN * unc_N)**2 +
            (dm_dlogM1 * unc_M1)**2 +  # Uncertainty in M1 is already in log10 units
            (dm_dbeta * unc_beta)**2 +
            (dm_dgamma * unc_gamma)**2
        )

        sigma_m = np.sqrt(sigma_m_squared)
        return sigma_m
