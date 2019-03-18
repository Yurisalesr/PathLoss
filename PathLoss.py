import numpy as np

class RMa_LOS:
    def __init__(self,fc,hbs,hut,h,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.h = h
        self.fc = fc/(10**9)
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        self.dbp = 2*np.pi*self.hbs*self.fc*10**9/(3*10**8)
        
    def calcPathloss(self):
        if self.d2d>=10 and self.d2d <= self.dbp:
            return (20*np.log10(40*np.pi*self.d3d*self.fc/3) + \
                           min(0.03*self.h**1.72,10)*np.log10(self.d3d) - \
                           min(0.044*self.h**1.72,14.77) + 0.002*np.log10(self.h)*self.d3d) +\
                           self.sf*np.random.randn(1) 
                
        
      
        if self.d2d>self.dbp and self.d2d<=10**3:
            return (20*np.log10(40*np.pi*self.dbp*self.fc/3) +\
                       min(0.03*self.h**1.72,10)*np.log10(self.dbp) - \
                       min(0.044*self.h**1.72,14.77) + 0.002*np.log10(self.h)*self.dbp\
                       + 40*np.log10(self.d3d/self.dbp)) + self.sf*np.random.randn(1)
            
            
class RMa_NLOS:
    def __init__(self,fc,hbs,hut,h,W,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.W = W
        self.h = h
        self.fc = fc/(10**9) 
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        self.dbp = 2*np.pi*self.hbs*self.fc*10**9/(3*10**8)
        
    def calcPathloss(self):
        self.PL_RMa_NLOS = 161.04 - 7.1*np.log10(self.W) +\
                      7.5*np.log10(self.h) -\
                      (24.37 - 3.7*(self.h/self.hbs)**2)*np.log10(self.hbs) +\
                      (43.42-3.1*np.log10(self.hbs))*(np.log10(self.d3d) -3) +\
                      20*np.log10(self.fc) -\
                      (3.2*(np.log10(11.75*self.hut))**2 - 4.97)
        if self.d2d>=10 and self.d2d <= self.dbp:
                self.PL1 = 20*np.log10(40*np.pi*self.d3d*self.fc/3) + \
                           min(0.03*self.h**1.72,10)*np.log10(self.d3d) - \
                           min(0.044*self.h**1.72,14.77) + 0.002*np.log10(self.h)*self.d3d
                self.PL_RMa_LOS = self.PL1
     
        if self.d2d>self.dbp and self.d2d<=10**3:
                self.PL2 = 20*np.log10(40*np.pi*self.dbp*self.fc/3) +\
                       min(0.03*self.h**1.72,10)*np.log10(self.dbp) - \
                       min(0.044*self.h**1.72,14.77) + 0.002*np.log10(self.h)*self.dbp\
                       + 40*np.log10(self.d3d/self.dbp)
                self.PL_RMa_LOS = self.PL2
                
        return max(self.PL_RMa_NLOS,self.PL_RMa_LOS) + self.sf*np.random.randn(1)
        

class UMa_LOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9) 
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        
    def calcPathloss(self):
            
        if self.hut <13:
            self.he = 1
            self.dbp = 4*np.pi*(self.hbs - self.he)*(self.hut - self.he)*self.fc*10**9/(3*10**8)
        
        if self.hut>=13 and self.hut<=23:
            if self.d2d <= 18:
                self.g = 0
            
            if self.d2d > 18:
                self.g = 5/4*(self.d2d/100)**3*np.exp(-self.d2d/150)
            
                self.C = ((self.hut - 13)/10)**1.5*self.g
                self.rand = np.random.rand(1)
            
            if self.rand < self.C:
                self.he = 1
            
            if self.rand >= self.C:
                self.he = np.random.choice(np.arange(12,self.hut-1.5,3))
            
        self.dbp = 4*np.pi*(self.hbs - self.he)*(self.hut - self.he)*self.fc*10**9/(3*10**8)
        
        if self.d2d >= 10 and self.d2d <= self.dbp:
            return (28 + 22*np.log10(self.d3d) + 20*np.log10(self.fc)) +\
                    self.sf*np.random.randn(1)
                
        if self.d2d > self.dbp and self.d2d<=(5*10**3):
            return (28 + 40*np.log10(self.d3d) + 20*np.log10(self.fc) -\
                       9*np.log10((self.dbp)**2 + (self.hbs - self.hut)**2)) \
                    + self.sf*np.random.randn(1)


class UMa_NLOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9)
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        
    def calcPathloss(self):
        self.PL_UMa_NLOS = 13.54 + 39.08*np.log10(self.d3d) + 20*np.log10(self.fc) -\
                      0.6*(self.hut - 1.5)
        if self.hut <13:
            self.he = 1
            self.dbp = 4*np.pi*(self.hbs - self.he)*(self.hut - self.he)*self.fc*10**9/(3*10**8)
        
        if self.hut>=13 and self.hut<=23:
                       
            if self.d2d <= 18:
                self.g = 0
            
            if self.d2d > 18:
                self.g = 5/4*(self.d2d/100)**3*np.exp(-self.d2d/150)
            
                self.C = ((self.hut - 13)/10)**1.5*self.g
                self.rand = np.random.rand(1)
            
            if self.rand < self.C:
                self.he = 1
            
            if self.rand >= self.c:
                self.he = np.random.choice(np.arange(12,self.hut-1.5,3))
            
        self.dbp = 4*np.pi*(self.hbs - self.he)*(self.hut - self.he)*self.fc*10**9/(3*10**8)
        
        if self.d2d >= 10 and self.d2d <= self.dbp:
            self.PL1 = 28 + 22*np.log10(self.d3d) + 20*np.log10(self.fc)
            self.PL_UMa_LOS = self.PL1
        if self.d2d > self.dbp and self.d2d<=5*10**3:
            self.PL2 = 28 + 40*np.log10(self.d3d) + 20*np.log10(self.fc) -\
                       9*np.log10((self.dbp)**2 + (self.hbs - self.hut)**2)
            self.PL_UMa_LOS = self.PL2
                
        return max(self.PL_UMa_LOS,self.PL_UMa_NLOS) + self.sf*np.random.randn(1)

class UMi_LOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9)
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        self.dbp = 4*np.pi*(self.hbs - 1)*(self.hut - 1)*self.fc/(3*10**8)
        
    def calcPathloss(self):
        if self.d2d <= self.dbp and self.d2d >= 10:
            return 32.4 + 21*np.log10(self.d3d) + 20*np.log10(self.fc) +\
                    self.sf*np.random.randn(1)
                        
        if self.d2d > self.dbp and self.d2d<= 5*10**3:   
            return  32.4 + 40*np.log10(self.d3d) + 20*np.log10(self.fc) \
                       -9.5*np.log10((self.dbp)**2 + (self.hbs - self.hut)**2) + \
                     self.sf*np.random.randn(1)
                
                        
class UMi_NLOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9) 
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        self.dbp = 4*np.pi*(self.hbs - 1)*(self.hut - 1)*self.fc/(3*10**8)
        
    def calcPathloss(self):
        self.PL_UMi_NLOS = 35.3*np.log10(self.d3d) + 22.4 + 21.3*np.log10(self.fc) -\
                               0.3*(self.hut-1.5)
        if self.d2d <= self.dbp and self.d2d >= 10:
            self.PL1 = 32.4 + 21*np.log10(self.d3d) + 20*np.log10(self.fc)
            self.PL_UMi_LOS = self.PL1
        
        if self.d2d > self.dbp and self.d2d<= 5*10**3:   
            self.PL2 = 32.4 + 40*np.log10(self.d3d) + 20*np.log10(self.fc) \
                       -9.5*np.log10((self.dbp)**2 + (self.hbs - self.hut)**2)
            self.PL_UMi_LOS = self.PL2
                
        return max(self.PL_UMi_NLOS,self.PL_UMi_LOS) + self.sf*np.random.randn(1)    

                        
class InH_LOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9)
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        
    def calcPathloss(self):
        return 32.4 + 17.3*np.log10(self.d3d) + 20*np.log10(self.fc) + self.sf*np.random.randn(1) 

class InH_NLOS:
    def __init__(self,fc,hbs,hut,sf):
        self.hbs = hbs[2]
        self.hut = hut[2]
        self.fc = fc/(10**9)
        self.sf = sf
        self.d2d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2)
        self.d3d = np.sqrt((hbs[0] - hut[0])**2 + (hbs[1] - hut[1])**2 + (hbs[2] - hut[2])**2)
        
    def calcPathloss(self):
        self.PL_InH_LOS = 32.4 + 17.3*np.log10(self.d3d) + 20*np.log10(self.fc)
        self.PL_InH_NLOS = 38.3*np.log10(self.d3d) + 17.3 + 24.9*np.log10(self.fc)
        return  max(self.PL_InH_LOS,self.PL_InH_NLOS) + self.sf*np.random.randn(1)