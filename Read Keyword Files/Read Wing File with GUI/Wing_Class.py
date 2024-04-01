import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class Wing:
    def __init__(self,):
        '''
        Constructor for the Wing class.
        '''
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
        '''
        data is a list of strings (i.e., the lines read from a file)
        :param data:
        :return:
        '''

        for line in data: #loop over all the lines
            line=line.strip() # first, strip all leading and trailing blank spaces
            cells=line.split(',') # then split using the ',' delimiter
            keyword=cells[0].lower().strip()
            if keyword=='title':
                self.title=cells[1].replace("'","")
            elif keyword=='distance_unit':
                self.dist_unit=cells[1].replace("'","")
            elif keyword=='force_unit':
                self.force_unit=cells[1].replace("'","")
            elif keyword=='tensile_strength':
                self.St=float(cells[1])
            elif keyword=='shear_strength':
                self.Ss=float(cells[1])
            elif keyword=='modulus':
                self.modulus=float(cells[1])
            elif keyword=='wing_length':
                self.length=float(cells[1])
            elif keyword=='wing_taper':
                self.taper=float(cells[1])
            elif keyword=='load_factor':
                self.load_factor=float(cells[1])
            elif keyword=='margin_of_safety':
                self.margin=float(cells[1])
            elif keyword=='spar_envelope':
                self.sparH=float(cells[1])    
                self.sparW=float(cells[2])
            elif keyword=='dist_up_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    cell=cell.replace('(','').replace(')','')
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.up_dist.append(this_load)
            elif keyword=='dist_aft_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.aft_dist.append(this_load)
            elif keyword=='point_up_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.up_point.append(this_load)
            elif keyword=='point_aft_load':
                this_load=[cells[1].replace("'","")]
                ncells=len(cells)
                for cell in cells[2:]:
                    value=float(cell.replace("(","").replace(")",""))
                    this_load.append(value)
                self.aft_point.append(this_load)

    def generate_report(self):
        report = '\tTitle: '+ self.title
        report += '\n\tWing length: '+ str(self.length)
        report += '\n\tLoad factor: '+ str(self.load_factor)

        report +='\n\n'
        report += '\n\t\t Distributed Loads - Upwards'

        for set in self.up_dist:
            report += '\n\nLoad Set: ' + set[0]
            for i in range(1,len(set),2):
                report += '\n' + str(set[i]) + "   " + str(set[i+1])

        return report






