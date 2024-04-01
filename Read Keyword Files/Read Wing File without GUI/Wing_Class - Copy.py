import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Wing:

    def __init__(self,):
        self.title = None
        self.dist_unit = None
        self.force_unit = None
        self.St = None
        self.Ss =  None
        self.modulus = None
        self.length = None
        self.taper = None
        self.load_factor = None
        self.margin = None
        self.sparH = None
        self.sparW = None
        self.up_dist = []
        self.aft_dist = []
        self.up_point = []
        self.aft_point = []

    
    def processWingData(self,data):
        #from the array of strings, fill the wing dictionary


        for line in data: #loop over all the lines
            line = line.strip().replace("'","").replace('"','').replace("(","").replace(")","")
            cells=line.split(',')
            keyword=cells[0].lower()
    
            if keyword=='title': self.title=cells[1].replace("'","")
            if keyword=='distance_unit': self.dist_unit=cells[1].replace("'","")
            if keyword=='force_unit': self.force_unit=cells[1].replace("'","")
            if keyword=='tensile_strength': self.St=float(cells[1])
            if keyword=='shear_strength': self.Ss=float(cells[1])
            if keyword=='modulus': self.modulus=float(cells[1])
            if keyword=='wing_length': self.length=float(cells[1])
            if keyword=='wing_taper': self.taper=float(cells[1])
            if keyword=='load_factor': self.load_factor=float(cells[1])
            if keyword=='margin_of_safety': self.margin=float(cells[1])
    
            if keyword=='spar_envelope':
                self.sparH=float(cells[1])    
                self.sparW=float(cells[2])
    
            if keyword=='dist_up_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.up_dist.append(this_load)
    
            if keyword=='dist_aft_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.aft_dist.append(this_load)
                    
            if keyword=='point_up_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.up_point.append(this_load)
    
            if keyword=='point_aft_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.aft_point.append(this_load)
