from Rankine_GUI import Ui_Form
from PyQt5 import uic
import sys
from PyQt5 import QtWidgets as qtw
from Rankine import rankine
from Steam import steam

#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """MainWindow constructor"""
        super().__init__()
        self.setupUi(self)
        # Main UI code goes here
        self.le_TurbineInletCondition.setEnabled(False)
        #creating a canvas to draw a figure for the rankine cycle
        self.figure=Figure(figsize=(3,8),tight_layout=True, frameon=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.gb_Output.layout().addWidget(self.canvas,6,0,1,6)

        #setting up some signals and slots
        self.le_PHigh.editingFinished.connect(self.setPHigh)  #triggered by hitting enter or leaving the line edit
        self.le_PLow.editingFinished.connect(self.setPLow)   #triggered by hitting enter or leaving the line edit
        self.le_TurbineEff.editingFinished.connect(self.checkTurbineEffRange)
        self.rdo_Quality.toggled.connect(self.setQualityOrTHigh) #triggered when the state of the radio button changes
        self.btn_Calculate.clicked.connect(self.calcRankine)
        # End main ui code

        #create a rankine object to work with later
        self.RC=rankine(8,8000,name='Default Rankine Cycle')
        #create a steam object to help with retrieving saturated properties
        self.satsteam=steam(8000,x=1.0)

        self.TSatHigh=0
        #call the functions to set the saturation properties during construction of this class
        self.setPHigh()
        self.setPLow()

        #show the form
        self.show()

    def clamp(self, val, low, high):
        if self.isfloat(val):
            val=float(val)
            if val>high:
                return float(high)
            if val <low:
                return float(low)
            return val
        return float(low)

    def isfloat(self,value):
        '''
        This function is a check to verify that a string can be converted to a float
        :return:
        '''
        if value=='NaN':return False
        try:
            float(value)
            return True
        except ValueError:
            return False

    def setPHigh(self):
        #make sure it is a number
        ph=self.le_PHigh.text()
        if not self.isfloat(ph):
            ph='80'
            self.le_PHigh.setText(ph)

        PHigh = float(ph) #convert text to number
        sat_high=self.satsteam.getSatProp(P_Bar=PHigh)
        #(Tsat, hf, hg, sf, sg, vf, vg)
        self.TSatHigh=sat_high[0]
        st_high='PSat = {:0.2f} bar, TSat = {:0.2f} C'.format(PHigh,sat_high[0])
        st_high +='\nhf = {:0.2f} kJ/kg, hg = {:0.2f} kJ/kg'.format(sat_high[1], sat_high[2])
        st_high +='\nsf = {:0.2f} kJ/kg*K, sg = {:0.2f} kJ/kg*k'.format(sat_high[3], sat_high[4])
        st_high +='\nvf = {:0.4f} m^3/kg, vg = {:0.2f} m^3/kg'.format(sat_high[5], sat_high[6])
        self.lbl_SatPropHigh.setText(st_high)

    def setPLow(self):
        #make sure it is a number
        pl=self.le_PLow.text()
        if not self.isfloat(pl):
            pl='0.08'
            self.le_PLow.setText(pl)

        PLow=float(self.le_PLow.text()) #convert text to number
        sat_low=self.satsteam.getSatProp(P_Bar=PLow)
        #(Tsat, hf, hg, sf, sg, vf, vg)
        st_low='PSat = {:0.2f} bar, TSat = {:0.2f} C'.format(PLow,sat_low[0])
        st_low +='\nhf = {:0.2f} kJ/kg, hg = {:0.2f} kJ/kg'.format(sat_low[1], sat_low[2])
        st_low +='\nsf = {:0.2f} kJ/kg*K, sg = {:0.2f} kJ/kg*k'.format(sat_low[3], sat_low[4])
        st_low +='\nvf = {:0.4f} m^3/kg, vg = {:0.2f} m^3/kg'.format(sat_low[5], sat_low[6])
        self.lbl_SatPropLow.setText(st_low)

    def checkTurbineEffRange(self):
        '''
        Makes sure turbine efficiency is in the range from 0 to 1
        :return:
        '''
        e=self.clamp(self.le_TurbineEff.text(),0.0,1.0)
        self.le_TurbineEff.setText(str(e))


    def setQualityOrTHigh(self):
        TF=self.rdo_Quality.isChecked()
        if TF:
            self.lbl_TurbineInletCondition.setText('Turbine Inlet: x=')
            self.le_TurbineInletCondition.setText(str(1.0))
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            self.lbl_TurbineInletCondition.setText('Turbine Inlet: T High =')
            self.le_TurbineInletCondition.setText('{:0.2f}'.format(self.TSatHigh+1))
            self.le_TurbineInletCondition.setEnabled(True)

    def calcRankine(self):
        '''
        This is called when the calculate button is clicked
        :return: nothing
        '''
        #read the high and low pressure isobar values.  no range checking.
        PHigh=float(self.le_PHigh.text())
        PLow=float(self.le_PLow.text())

        #create a new rankine object with values depending on which radio buttton checked
        if(self.rdo_Quality.isChecked()):
            self.RC.set(p_low=PLow*100, p_high=PHigh*100, eff_turbine=float(self.le_TurbineEff.text()))
        else:
            self.RC.set(p_low=PLow*100, p_high=PHigh*100,eff_turbine=float(self.le_TurbineEff.text()), t_high=float(self.le_TurbineInletCondition.text()))
        #calculate the cycle efficiency (and states 1,2,3,4)
        self.RC.calc_efficiency()

        #fill out the enthalpy values
        self.le_H1.setText('{:0.2f}'.format(self.RC.state1.h))
        self.le_H2.setText('{:0.2f}'.format(self.RC.state2.h))
        self.le_H3.setText('{:0.2f}'.format(self.RC.state3.h))
        self.le_H4.setText('{:0.2f}'.format(self.RC.state4.h))

        #fill out the other properties for the rankine cycle
        self.le_Efficiency.setText('{:0.2f}'.format(self.RC.efficiency))
        self.le_TurbineWork.setText('{:0.2f}'.format(self.RC.turbine_work))
        self.le_PumpWork.setText('{:0.2f}'.format(self.RC.pump_work))
        self.le_HeatAdded.setText('{:0.2f}'.format(self.RC.heat_added))

        self.ax.clear()
        self.RC.plot_cycle_TS(self.ax)

        self.canvas.draw()

#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine Cycle Calculator')
    sys.exit(app.exec())