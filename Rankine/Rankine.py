from Steam import steam
from matplotlib import pyplot as plt
import numpy as np

class rankine():
    def __init__(self, p_low=8, p_high=8000, t_high=None, eff_turbine=1.0, name='Rankine Cycle'):
        '''
        Constructor for rankine power cycle.  If t_high is not specified, the State 1
        is assigned x=1 (saturated steam @ p_high).  Otherwise, use t_high to find State 1.
        :param p_low: the low pressure isobar for the cycle in kPa
        :param p_high: the high pressure isobar for the cycle in kPa
        :param t_high: optional temperature for State1 (turbine inlet) in degrees C
        :param eff_turbine: the turbine efficiency eta=(h1-h2)/(h1-h2s)<=1.0
        :param name: a convenient name
        '''
        self.set(p_low=p_low, p_high=p_high, t_high=t_high, eff_turbine=eff_turbine, name=name)
        self.steam=steam(8000)
        self.tsdata=self.steam.getVaporDome_TS(100)

    def set(self, p_low=8, p_high=8000, t_high=None, eff_turbine=1.0, name='Rankine Cycle'):
        self.p_low = p_low
        self.p_high = p_high
        self.t_high = t_high
        self.name = name
        self.efficiency = None
        self.eff_turbine = eff_turbine
        self.turbine_work = 0
        self.pump_work = 0
        self.heat_added = 0
        self.state1 = None
        self.state2 = None
        self.state2s = None
        self.state3 = None
        self.state4 = None

    def calc_efficiency(self):
        # calculate the 4 states
        # state 1: turbine inlet (p_high, t_high) superheated or saturated vapor
        if (self.t_high == None):
            self.state1 = steam(self.p_high, x=1.0, name='Turbine Inlet')
        else:
            self.state1 = steam(self.p_high, T=self.t_high, name='Turbine Inlet')
        # state 2: turbine exit (p_low, s=s_turbine inlet) two-phase
        # create state 2s for 100% efficient turbine
        self.state2s = steam(self.p_low, s=self.state1.s, name="Turbine Exit")
        # use turbine efficiency to calculate h2
        h2 = self.state1.h - (self.state1.h - self.state2s.h) * self.eff_turbine
        # finally, find state 2
        self.state2 = steam(self.p_low, h=h2, name="Turbine Exit")
        # state 3: pump inlet (p_low, x=0) saturated liquid
        self.state3 = steam(self.p_low, x=0, name='Pump Inlet')
        # state 4: pump exit (p_high,s=s_pump_inlet) typically sub-cooled, but estimate as saturated liquid
        self.state4 = steam(self.p_high, s=self.state3.s, name='Pump Exit')
        self.state4.h = self.state3.h + self.state3.v * (self.p_high - self.p_low)

        self.turbine_work = self.state1.h - self.state2.h
        self.pump_work = self.state4.h - self.state3.h
        self.heat_added = self.state1.h - self.state4.h
        self.efficiency = 100.0 * (self.turbine_work - self.pump_work) / self.heat_added
        return self.efficiency

    def print_summary(self):
        if self.efficiency == None:
            self.calc_efficiency()
        print('Cycle Summary for: ', self.name)
        print('\tEfficiency: {:0.3f}%'.format(self.efficiency))
        print('\tTurbine Work: {:0.3f} kJ/kg'.format(self.turbine_work))
        print('\tPump Work: {:0.3f} kJ/kg'.format(self.pump_work))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(self.heat_added))
        self.state1.print()
        self.state2.print()
        self.state3.print()
        self.state4.print()

    def get_summary(self):
        '''
        This returns a formatted string to put on the plot of the rankine cycle.
        :return:
        '''
        if self.efficiency == None:
            self.calc_efficiency()
        s = r'Summary:'
        s += '\n$\eta$: {:0.1f}% '.format(self.efficiency)
        s += '\n$\eta_{turbine}$: ' + '{:0.2f}'.format(self.eff_turbine) if self.eff_turbine < 1.0 else ''
        s += '\n$W_{turbine}$: ' + '{:0.1f} kJ/k'.format(self.turbine_work)
        s += '\n$W_{pump}$: ' + '{:0.1f} kJ/kg'.format(self.pump_work)
        s += '\n$Q_{boiler}$: ' + '{:0.1f} kJ/kg'.format(self.heat_added)
        return s

    def plot_cycle_TS(self, ax=None):
        '''
        1. generate the data needed for mapping the vapor dome (i.e., saturated liquid and saturated vapor)
        2. generate (s, T) data for the p_high isobar from sat liq. to t_high.
        3. generate (s, T) data for the p_low isobar from sat liq to s_2.
        4. generate (s, T) data to linearly interpolate from (s_sat_f, p_low) to (s_sat_f, p_high).
        5. generate (s, T) data to linearly interpolate from (s_1, p_high) to (s_2, p_low)
        6. splice together data from 4+2+5
        :return:
        '''
        QTPlotting=True #assumes we are plotting onto a QT GUI form
        if ax==None:
            ax=plt.subplot()
            QTPlotting=False  #actually, we are just using CLI and showing the plot
        # 1 generate the data needed for mapping the vapor dome
        # one way to do it
        tsdata = self.tsdata  #(self.p_high).getVaporDome_TS(points=100)
        T = tsdata[:, 0]
        SF = tsdata[:, 1]
        SG = tsdata[:, 2]
        # another way to do it
        #T,SF,SG=steam().getVaporDome_TS2()
        # actually plot the vapor dome
        ax.plot(SF, T, color='b')
        ax.plot(SG, T, color='r')

        # 2 generate generate (s, T) data for the p_high isobar from sat liq. to t_high.
        p_highdata = np.empty((0, 2))  # create an empty, single axis numpy array with two columns and zero rows to start
        satsteam=self.state1 #create a steam object
        sf = satsteam.satProp.sf # get sf
        tsat=satsteam.satProp.Tsat # get tsat
        sg = satsteam.satProp.sg # get sg
        #build up numpy array of (tsat, s)
        for i in range(100):
            p_highdata = np.append(p_highdata, [[tsat, sf + i / 100.0 * (sg - sf)]], axis=0)
        #make sure to include value for saturated liquid
        p_highdata = np.append(p_highdata, [[tsat, sf]], axis=0)
        #build data for superheated if necessary
        if self.state1.T > tsat:
            for i in range(20):
                #interpolate between tsat+1 to state1.T
                t = tsat + 1 + i / 20.0 * (self.state1.T - tsat - 1)
                p_highdata = np.append(p_highdata, [[t, steam(pressure=self.p_high, T=t).s]], axis=0)
            #make sure it includes state 1
            p_highdata = np.append(p_highdata, [[self.state1.T, self.state1.s]], axis=0)

        # 3 generate generate (s, T) data for the p_low isobar from sat liq. to s_2.
        p_lowdata = np.empty((0, 2))
        sf = self.state3.satProp.sf  #already known, so use state 3
        tsat = self.state3.satProp.Tsat
        for i in range(100):
            p_lowdata = np.append(p_lowdata, [[tsat, sf + i / 100.0 * (self.state2.s - sf)]], axis=0)
        p_lowdata = np.append(p_lowdata, [[tsat, self.state2.s]], axis=0)
        # add low pressure isobar T-s data to the plot
        ax.plot(p_lowdata[:, 1], p_lowdata[:, 0], color='k')

        # 4 generate (s, T) data to linearly interpolate from (s_sat_f, p_low) to (s_sat_f, p_high).
        shigh = self.state1.satProp.sf #p_highdata[0, 1]
        slow = sf #from above for #3
        thigh = self.state1.satProp.Tsat #p_highdata[0, 0]
        tlow = tsat #from above for #3
        scl = np.empty((0, 2))
        for i in range(100):
            scl = np.append(scl, [[tlow + i / 100.0 * (thigh - tlow), sf + i / 100.0 * (shigh - slow)]], axis=0)
        scl = np.append(scl, [[thigh, shigh]], axis=0)
        p_highdata = np.append(scl, p_highdata, axis=0)

        # 5. generate (s, T) data to linearly interpolate from (s_1, p_high) to (s_2, p_low)
        TH, SH, TL, SL =self.state1.T, self.state1.s, self.state2.T, self.state2.s
        for i in range(100):
            p_highdata = np.append(p_highdata, [[TH + i / 100.0 * (TL - TH), SH + i / 100.0 * (SL - SH)]], axis=0)
        p_highdata = np.append(p_highdata, [[TL, SL]], axis=0)
        # temperature values for a straight line in T-s space across the vapor dome connecting state 3 with state 2

        tdatalow = np.array([tlow for t in p_highdata[:, 0]])
        ax.plot(p_highdata[:, 1], p_highdata[:, 0], color='g')
        ax.fill_between(p_highdata[:, 1], p_highdata[:, 0], tdatalow, color='grey', alpha=0.2)

        #add axis labels
        ax.set_ylabel(r'T ($^oC$)', fontsize='large' if QTPlotting else 'medium')
        ax.set_xlabel(r'S $\left(\frac{kJ}{kg\cdot K}\right)$', fontsize='large' if QTPlotting else 'medium')
        #put a title on the plot
        self.name='Rankine Cycle - ' + self.state1.region + ' at Turbine Inlet'
        ax.set_title(self.name, fontsize='large' if QTPlotting else 'medium')
        #add the summary text to the plot
        ax.text(0.5, 350, self.get_summary())
        # #get reference to the current axes of the plot with gca()
        # ax = plt.gca()
        #modify the tick marks (not required for the exam)
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize='large' if QTPlotting else 'medium')  # format tick marks

        #plot the circles for states 1, 2, and 3
        ax.plot(self.state1.s, self.state1.T, marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(self.state2.s, self.state2.T, marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(self.state3.s, self.state3.T, marker='o', markerfacecolor='w', markeredgecolor='k')
        #set limits on x and y
        ax.set_xlim(tsdata[:, 1].min(), max(tsdata[:, 2].max(), self.state1.s))
        ax.set_ylim(0, 550)
        #show the plot
        if QTPlotting==False:
            plt.show()

def main():
    rankine1 = rankine(8, 8000, t_high=500, eff_turbine=0.95, name='Rankine Cycle - Superheated at turbine inlet')
    # t_high is specified
    # if t_high were not specified, then x_high = 1 is assumed
    eff = rankine1.calc_efficiency()
    print(eff)
    rankine1.state3.print()
    rankine1.print_summary()
    rankine1.plot_cycle_TS()
    # hf=rankine1.state1.hf
    # hg=rankine1.state1.hg
    rankine2 = rankine(8, 8000, eff_turbine=0.95, name='Rankine Cycle - Saturated at turbine inlet')
    eff2 = rankine2.calc_efficiency()
    rankine2.plot_cycle_TS()
    print(eff2)

    rankine2.print_summary()


if __name__ == "__main__":
    main()
