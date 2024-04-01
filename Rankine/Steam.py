import numpy as np
from scipy.interpolate import griddata
class satProps():
    def __init__(self, vals=None):
        self.Tsat=None
        self.hf=None
        self.hg=None
        self.sf=None
        self.sg=None
        self.vf=None
        self.vg=None

    def set(self, vals):
        self.Tsat, self.hf, self.hg, self.sf, self.sg,self.vf, self.vg=vals

    def get(self):
        return (self.Tsat, self.hf, self.hg, self.sf, self.sg, self.vf, self.vg)

class steam():
    def __init__(self, pressure, T=None, x=None, v=None, h=None, s=None, name=None):
        self.p = pressure #pressure - kPa
        self.T = T #Temperature - degrees C
        self.x = x #quality
        self.v = v #specific volume - m^3/kg
        self.h = h #specific enthalpy - kj/kg
        self.s = s #entropy - kj/(kg*K)
        self.satProp=satProps()
        self.name = name #a useful identifier
        self.region = None #'superheated' or 'saturated' or 'two-phase'
        if T==None and x==None and v==None and h==None and s==None: return
        else: self.calc()

    def calc(self):
        '''
        The Rankine cycle operates between two isobars (i.e., p_high (Turbine inlet state 1) & p_low (Turbine exit state 2)
        So, given a pressure, we need to determine if the other given property puts
        us in the saturated or superheated regime.
        :return: nothing returned, just set the properties
        '''
        #1. need to determine which second property is known
        #2. determine if two-phase/saturated or superheated
        #3. find all unknown thermodynamic properties by interpolation from appropriate steam table

        #read in the thermodynamic data from files
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs= np.loadtxt('sat_water_table.txt',skiprows=1, unpack=True)
        tcol, hcol, scol, pcol = np.loadtxt('superheated_water_table.txt', skiprows=1, unpack=True)

        R=8.314/(18/1000) #ideal gas constant for water [J/(mol K)]/[kg/mol]
        Pbar=self.p/100 #pressure in bar

        #get saturated properties
        #note: these are local variables only.  Their scope is only in this function
        Tsat = float(griddata((ps), ts, (Pbar)))
        hf=float(griddata((ps),hfs,(Pbar)))
        hg=float(griddata((ps),hgs,(Pbar)))
        sf=float(griddata((ps),sfs,(Pbar)))
        sg=float(griddata((ps),sgs,(Pbar)))
        vf=float(griddata((ps),vfs,(Pbar)))
        vg=float(griddata((ps),vgs,(Pbar)))
        self.satProp.set((Tsat, hf, hg, sf, sg, vf, vg))
        #self.hf=hf #this creates a member variable for the class that can be accessed from an object

        if self.T!=None:
            if self.T>Tsat:
                self.region='Superheated'
                self.h = float(griddata((tcol, pcol), hcol, (self.T, self.p)))
                self.s = float(griddata((tcol, pcol), scol, (self.T, self.p)))
                self.x=1.0
                TK = self.T + 273.14  # temperature conversion to Kelvin
                self.v=R*TK/(self.p*1000)  #ideal gas approximation for volume
        elif self.x!=None:
            self.region='Saturated'
            self.T=Tsat
            self.h=hf+self.x*(hg-hf)
            self.s=sf+self.x*(sg-sf)
            self.v=vf+self.x*(vg-vf)
        elif self.h!=None:
            self.x=(self.h-hf)/(hg-hf)
            if self.x<=1.0:
                self.region='Saturated'
                self.T=Tsat
                self.s=sf+self.x*(sg-sf)
                self.v=vf+self.x*(vg-vf)
            else:
                self.region='Superheated'
                self.T = float(griddata((hcol, pcol), tcol, (self.h, self.p)))
                self.s = float(griddata((hcol, pcol), scol, (self.h, self.p)))
        elif self.s!=None:
            self.x=(self.s-sf)/(sg-sf)
            if self.x<=1.0:
                self.region='Saturated'
                self.T=Tsat
                self.h=hf+self.x*(hg-hf)
                self.v=vf+self.x*(vg-vf)
            else:
                self.region = 'Superheated'
                self.T = float(griddata((scol, pcol), tcol, (self.s, self.p)))
                self.h = float(griddata((scol, pcol), scol, (self.s, self.p)))
        pass

    def getVaporDome_TS(self, points=500):
        '''
        This function gets the TS data for the vapor dome using a cubic interpolation from saturated steam tables.
        :param points: how many data points I want
        :return: a numpy array with three columns of T, sf, sg
        '''
        #read the saturated water table into columns
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs= np.loadtxt('sat_water_table.txt',skiprows=1, unpack=True)
        T_min=ts.min()
        T_max=ts.max()
        TSData=np.empty((0,3)) #create the empty array for return
        T=np.linspace(T_min,T_max, points)  #create the linspace with then number of points I want
        for t in T:
            dat=np.array([[t,float(griddata((ts),sfs,(t), method='cubic')),float(griddata((ts),sgs,(t),method='cubic'))]])
            TSData=np.append(TSData,dat,axis=0)
        return TSData

    def getVaporDome_TS2(self):
        '''
        This function just returns the tuple of (T, sf, sg) directly from the saturated steam table.
        :return:
        '''
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs= np.loadtxt('sat_water_table.txt',skiprows=1, unpack=True)
        return (ts, sfs, sgs)

    def getSatProp(self, P_KPa=None, P_Bar=None):
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)
        Pbar = P_Bar if P_Bar!=None else P_KPa / 100  # get pressure in bar

        # get saturated properties
        Tsat = float(griddata((ps), ts, (Pbar)))
        hf = float(griddata((ps), hfs, (Pbar)))
        hg = float(griddata((ps), hgs, (Pbar)))
        sf = float(griddata((ps), sfs, (Pbar)))
        sg = float(griddata((ps), sgs, (Pbar)))
        vf = float(griddata((ps), vfs, (Pbar)))
        vg = float(griddata((ps), vgs, (Pbar)))
        self.satProp.set((Tsat, hf, hg, sf, sg, vf, vg))
        return self.satProp.get()

    def print(self):
        print('Name: ', self.name)
        if self.x<0.0: print('Region: compressed liquid')
        else: print('Region: ', self.region)
        print('p = {:0.2f} kPa'.format(self.p))
        if self.x >= 0.0: print('T = {:0.1f} degrees C'.format(self.T))
        print('h = {:0.2f} kJ/kg'.format(self.h))
        if self.x >= 0.0:
            print('s = {:0.4f} kJ/(kg K)'.format(self.s))
            if self.region == 'Saturated': print('v = {:0.6f} m^3/kg'.format(self.v))
            if self.region == 'Saturated': print('x = {:0.4f}'.format(self.x))
        print()

def main():
    inlet=steam(7350,name='Turbine Inlet') #not enough information to calculate
    inlet.x=0.9 #90 percent quality
    inlet.calc()
    inlet.print()

    h1=inlet.h
    s1=inlet.s
    print(h1,s1,'\n')

    outlet=steam(100, s=inlet.s, name='Turbine Exit')
    outlet.print()

    another=steam(8575, h=2050, name='State 3')
    another.print()

    yetanother = steam(8575, h=3125, name='State 4')
    yetanother.print()

#the following if statement causes main() to run
#only if this file is being run explicitly, not if it is
#being imported into another Python program as a module
if __name__=="__main__":
    main()
