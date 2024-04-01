#region imports
import numpy as np
import math
from copy import copy, deepcopy
from scipy.optimize import fsolve
import PyQt5.QtWidgets as qtw
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
#endregion

#region class definitions
class Fluid:
    def __init__(self, mu=0.00089, rho=1000.0, SI=True):
        """
        default properties are for water
        :param mu: dynamic viscosity in Pa*s -> (kg*m/s^2)*(s/m^2) -> kg/(m*s) or (lb*s/ft^2)
        :param rho: density in kg/m^3
        :param SI: tells constructor if unit conversion is needed from english to si if SI==False
        """
        #(lb*s)/ft^2*((3.3ft/m)^2)*(1kg/2.2lb)*(9.81m/s^2)->(Pa*s)
        self.mu=mu if SI==True else mu*(3.3**2)*(1/2.2)*9.81
        #(lb/ft^3)*((3.3ft/m)^3)*(1kg/2.2lb) -> kg/m^3
        self.rho=rho if SI==True else rho*(3.3**3)/(2.2)
        self.nu=self.mu/self.rho  # kinematic viscosity in units of m^2/s

    def m_to_psi(self, p):
        #convert m of water to psi
        #(m)*(3.3*12in/m)*rho(kg/m^3)*(2.2lb/kg)*(1m/(3.3*12in))^3
        psi=p*self.rho*2.2/((3.3*12)**2)
        return psi

    def psi_to_m(self, p):
        #convert psi to m of water
        #(lb/in^2)*(1kg/2.2lb)*(3.3*12in/m)^2*(1/rho)(m^3/kg)
        m=p*(1/2.2)*((3.3*12)**2)*(1/self.rho)
        return m

class Node():
    def __init__(self, name=None, pipes=None, ext_flow=None, z=None, specifiedP=None, min_ph=None, oldnode=None):
        """
        A node in a pipe network.
        :param name: name of the node
        :param pipes: a list/array of pipe objects connected to this node
        :param ext_flow: any external flow into (+) or out (-) of this node in L/s
        :param z: elevation of the node in m
        """
        self.Name = name if oldnode is None else oldnode.getName()
        self.Pipes = pipes if oldnode is None else oldnode.Pipes
        self.E = 0.0 if oldnode is None else oldnode.E
        self.E = ext_flow if ext_flow is not None else self.E
        self.P = 0.0  # the pressure head at the node in m of fluid
        self.SpecifiedP = 0.0 if oldnode is None else oldnode.SpecifiedP
        self.SpecifiedP = specifiedP if specifiedP is not None else self.SpecifiedP
        self.Z = 0.0 if oldnode is None else oldnode.Z
        self.Z = z if z is not None else self.Z #elevation of the node in m
        self.IsSprinkler = False
        self.K=1.0
        self.MinPH = 0.0 if oldnode is None else oldnode.MinPH
        self.MinPH = min_ph if min_ph is not None else self.MinPH

    def getNetFlowRate(self):
        """
        Calculates the net flow rate into this node in L/s
        :return:
        """
        Qtot=self.E  #count the external flow first
        for p in self.Pipes:
            #retrieves the pipe flow rate (+) if into node (-) if out of node.  see class for pipe.
            Qtot+=p.getFlowIntoNode(self.Name)
        return Qtot

    def getName(self):
        return self.Name

    def modifyNode(self, name=None, pipes=None, ext_flow=None, z=None, specifiedP=None, min_ph=None):
        if not name==None:
            self.Name=name
        if not pipes==None:
            self.Pipes=pipes
        if not ext_flow==None:
            self.E=ext_flow
        if not z==None:
            self.Z=z
        if not specifiedP==None:
            self.SpecifiedP=specifiedP
        if not min_ph==None:
            self.MinPH=min_ph

    def setExtFlow(self, E, SI=True):
        #(ft^3/s)*(1m/3.3ft)^3*(1000L/m^3)->L/s
        self.E=E if SI else 1000*E*(1/3.3)**3

class SprinklerHead(Node):
    def __init__(self, name=None, pipes=None, ext_flow=None, z=None, specifiedP=None, min_ph=None, k=None, oldnode=None, SI=True):
        """
        SprinklerHead inherits from node (i.e., it is a special kind of node).
        :param name: inherits
        :param pipes: inherits
        :param ext_flow: inherits
        :param z: inherits
        :param min_ph: minimum pressure head for sprinkler
        :param k: discharge coefficient k=Q/sqrt(P)
        :param oldnode: a node object that I want to copy properties from if it is specified
        """
        if oldnode!=None:
            #if I pass an object from which to inherit properties.
            #run parent constructor
            super().__init__(oldnode=oldnode)
        else:
            #run parent constructor
            super().__init__(name=name, pipes=pipes,ext_flow=ext_flow, z=z, specifiedP=specifiedP, min_ph=min_ph)
        self.E=ext_flow if ext_flow is not None else self.E
        self.setExtFlow(self.E, SI)
        self.Z=z if z!=None else self.Z
        self.minPH = min_ph if min_ph is not None else 2
        self.minPH=self.minPH if SI else self.minPH/3.3
        self.K=1.0 if oldnode is None else oldnode.K
        self.K = k if k is not None else self.K
        self.IsSprinkler = True

    def calc_k(self):  # calculate sprinkler discharge coefficient
        self.K=abs(self.E/math.sqrt(self.P))
        return self.K

    def calc_q (self):  # calculate outflow of sprinkler.  Warning: P should be positive to avoid error.
        self.E = self.K*math.sqrt(self.P)
        return self.E

class Pipe:
    def __init__(self, start='A', end='B', length=100.0, dia=200.0, roughness=0.00025, fluid=Fluid(), SI=True):
        """
        Defines a generic pipe with orientation from lowest letter to highest.
        :param start: the start node
        :param end: the end node
        :param length: the pipe length in m
        :param dia: the pipe diameter in mm
        :param roughness: the pipe roughness in m
        :param SI: if SI==False, need to convert len, roughness from ft to m and dia from in to m
        """
        self.modifyPipe(start=start,end=end,length=length, dia=dia, roughness=roughness,fluid=fluid, SI=SI)

    def modifyPipe(self, start='A', end='B', length=100.0, dia=200.0, roughness=0.00025, fluid=Fluid(), SI=True):
        self.startNode = min(start, end)  # makes sure to use the lowest letter for startNode
        self.endNode = max(start, end)  # makes sure to use the highest letter for the endNode
        self.setLength(length,SI)
        self.setDiam(dia, SI)
        self.setRough(roughness,SI)
        self.Q = 10  # working in units of L/s.  Signed flow rate.
        self.fluid = fluid  # the fluid in the pipe
        self.vel = self.v()  # calculate the initial velocity of the fluid
        self.reynolds = self.Re()  # calculate the initial reynolds number

    def setDiam(self, dia, SI=True):
        self.diam = dia / 1000.0 if SI else dia / (12 * 3.3)  # pipe diameter in (m)

    def setRough(self, roughness, SI=True):
        self.rough = roughness if SI else roughness / 3.3  # pipe roughness in (m)

    def setLength(self, length, SI=True):
        self.length = length if SI else length / 3.3  # pipe length in (m)

    def v(self):
        """
        Calculate average velocity for self.Q
        :return: average velocity
        """
        A = math.pi / 4.0 * self.diam ** 2
        self.vel=(abs(self.Q)/1000.0)/A # need to convert Q to m^3/s
        return self.vel

    def Re(self):
        """
        Calculate the reynolds number under current conditions.
        :return:
        """
        self.reynolds= self.v() * self.diam / self.fluid.nu
        return self.reynolds

    def FrictionFactor(self):
        """
        Use the Colebrook equation to find the friction factor.
        NOTE:  math.log is natural log, math.log10 is base 10 log
        """
        relrough=self.rough/self.diam
        def ffc(ff):  # friction factor calculator
            LHS = 1 / (ff ** 0.5)
            RHS = -2.0 * math.log10(relrough / 3.7 + 2.51 / (self.Re() * ff ** 0.5))
            return LHS - RHS
        f = fsolve(ffc, np.array([0.008]))  # use fsolve to find friction factor
        return f[0]

    def HeadLoss(self):  # calculate headloss through a section of pipe in m of fluid
        """
        Use the Darcy-Weisbach equation to find the head loss through a section of pipe.
        """
        g = 9.81  # m/s^2
        ff = self.FrictionFactor()
        hl = ff * (self.length / self.diam) * (self.v() ** 2) / (2 * g)  # m of water
        return hl

    def getHeadLoss(self, s):
        """
        Calculate the head loss for the pipe in direction of loop traversal.
        :param s: the node i'm starting with in a traversal of the pipe
        :return: the signed headloss through the pipe in m of fluid
        """
        #while traversing a loop, if s = startNode I'm traversing in same direction as positive pipe flow
        nTraverse= 1 if s==self.startNode else -1
        #if flow is positive sense, scalar =1 else =-1
        nFlow=1 if self.Q >= 0 else -1
        return nTraverse*nFlow*self.HeadLoss()

    def getName(self):
        """
        Gets the pipe name.
        :return: pipe name (e.g., 'a-b')
        """
        return self.startNode+'-'+self.endNode

    def oContainsNode(self, node):
        #does the pipe connect to the node?
        return self.startNode==node or self.endNode==node

    def printPipeFlowRate(self, SI=True):
        q_units = 'L/s' if SI else 'cfs'
        q=self.Q if SI else self.Q/1000*(3.3**3)
        print('Q for {} = {:0.2f} {}'.format(self.getName(), q, q_units))

    def getPipeFlowRateOutput(self, SI=True, tabs=''):
        q_units = 'L/s' if SI else 'cfs'
        q=self.Q if SI else self.Q/1000*(3.3**3)
        return tabs+'Q in {} = {:0.2f} {}'.format(self.getName(), q, q_units)

    def getFlowIntoNode(self, n):
        """
        determines the flow rate into node n
        :param n: a node object
        :return: +/-Q
        """
        return -self.Q if n==self.startNode else self.Q

class Loop():
    def __init__(self, Name='A', loopPipes = None):
        """
        Defines a loop in a pipe network.
        :param Name: name of the loop
        :param loopPipes: a list/array of pipe objects in this loop
        """
        self.name=Name
        self.pipes=loopPipes if not loopPipes==None else []

    def getName(self):
        return self.name

    def getPipes(self):
        return self.pipes

    def getNode(self, nodes, name):
        #get a node object from list of nodes by name
        return next(n for n in nodes if n.getName() == name)

    def getLoopHeadLoss(self, nodes):
        """
        Calculates the net head loss as I traverse around the loop, in m of fluid.
        Also, calculates pressure head at each node.
        :param nodes: pass along the list of nodes for the pipe network so I can find pressures at nodes
        :return:
        """
        deltaP = 0 # initialize loop pressure drop to zero
        sn_n = self.pipes[0].startNode # name of node at the start of the first pipe
        sn_o = self.getNode(nodes, sn_n) # get the stating node object
        ph = sn_o.P  # pressure head at startNode
        for p in self.pipes:
            # calculates the flow head loss in the pipe considering loop traversal and flow directions
            phl=p.getHeadLoss(sn_n)
            # determine next node in the loop
            nn_n = p.endNode if p.startNode == sn_n else p.startNode
            # gets the node object at other end of current pipe
            nn_o = self.getNode(nodes, nn_n)
            #calculate head loss due to elevation change
            deltaZ = nn_o.Z-sn_o.Z
            # update the loop head loss
            deltaP += (phl+deltaZ)
            #calc pressure head at end node
            ph -= (phl + deltaZ)
            #set pressure head at end node
            nn_o.P = ph
            # set ns node object to ne node object
            sn_o = nn_o
            # set start node (name) to the next node (name)
            sn_n = nn_n
        return deltaP

class PipeNetwork_model():
    def __init__(self, Pipes=None, Loops=None, Nodes=None, fluid=Fluid(), SI=True):
        """
        The pipe network is built from pipe, node, loop, and fluid objects.
        :param Pipes: a list of pipe objects
        :param Loops: a list of loop objects
        :param Nodes: a list of node objects
        :param fluid: a fluid object
        """
        self.pipes=Pipes if not Pipes==None else []
        self.loops=Loops if not Loops==None else []
        self.nodes=Nodes if not Nodes==None else []
        self.Fluid=fluid
        self.SI=SI

    def getNodePipes(self, node):
        #returns a list of pipe objects that are connected to the node object
        l=[p for p in self.pipes if p.oContainsNode(node)] #a list comprehension
        return l

    def buildNodes(self):
        #automatically create the node objects by looking at the pipe ends
        for p in self.pipes:
            if self.hasNode(p.startNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.startNode))
            if self.hasNode(p.endNode)==False:
                #instantiate a node object and append it to the list of nodes
                self.nodes.append(Node(p.endNode))
        for n in self.nodes:
            n.Pipes=self.getNodePipes(n.getName())

    #region Functions for evaluating the pipe network
    def findFlowRates(self):
        """
        A method to analyze the pipe network and find the flow rates in each pipe
        given the constraints of no net flow into a node and no net pressure drops in the loops.
        Here, the sprinkler min pressure is set to 2m and the output of each sprinkler head is fixed
        so that we can calculate the discharge coefficient for each sprinkler.
        :return: a list of flow rates in the pipes
        """
        #see how many nodes and loops there are
        N=len(self.nodes)+len(self.loops)
        #build an initial guess for flow rates
        Q0=np.full(N,5)
        self.iterations=0
        def fn(q): #the callback for fsolve
            #update the flow rate in each pipe object
            for i in range(len(self.pipes)):
                self.pipes[i].Q=q[i]
            # calculate the net head loss for the loop objects
            LHL = self.getLoopHeadLosses()
            #calculate the net flow rate for the node objects
            NFR = self.getNodeFlowRates()
            self.iterations+=1
            return LHL+NFR
        #using fsolve to find the flow rates
        FR=fsolve(fn,Q0)
        return FR

    def findFlowRates2(self):
        """
        A method to analyze the pipe network and find the flow rates in each pipe.
        Given the constraints of: i) no net flow into a node, ii) no net pressure drop in the loops and
        iii) no mass accumulation in the pipe network.
        I used the last variable in q for the minimum sprinkler pressure.
        I found that setting this value at 2m, didn't allow system to converge, so I let fsolve
        vary the min sprinkler pressure and got convergence.
        :return: a list of flow rates in the pipes + other variables
        """
        #see how many nodes and loops there are +1 for network mass balance
        N=len(self.nodes)+len(self.loops)+1
        #build an initial guess for flow fsolve variables
        Q0=np.full(N,1)
        self.iterations=0
        def fn(q): #the callback for fsolve
            #update the flow rate in each pipe object
            for i in range(len(self.pipes)):
                self.pipes[i].Q=q[i]

            #calculate the net head loss for the loop objects
            LHL = self.getLoopHeadLosses()
            #update the sprinkler flow rates and calculate net outflow for the network
            minP=abs(q[len(q)-1]) #allows fsolve to set min sprinkler pressure
            NFO = self.getNetOutFlow(minP)
            #calculate the net flow rate for the node objects
            NNF = self.getNodeFlowRates()
            self.iterations+=1
            return NFO+LHL+NNF
        #using fsolve to find the flow rates
        FR=fsolve(fn,Q0,factor=1000)
        NFO=self.getNetOutFlow(setmin=False)  #debugging check
        return FR

    def getNetOutFlow(self, minP=2.0, setmin=True):
        """
        Calculate the mass balance for the pipe network.
        :param minP: the minimum pressure head at a sprinkler in (m).
        :param setmin: if True, set min pressure for sprinkler by increasing pump pressure.
        :return: [net mass balance]
        """
        if setmin:
            self.setMinNodePressureHead(minP, calcK=False)
        netOut=0
        for n in self.nodes:
            if n.IsSprinkler:
                n.E = -1.0*n.K*math.sqrt(n.P)
            netOut+=n.E
        return [netOut]

    def getNodeFlowRates(self):
        #each node object is responsible for calculating its own net flow rate
        fr=[n.getNetFlowRate() for n in self.nodes]
        return fr

    def getLoopHeadLosses(self):
        #each loop object is responsible for calculating its own net head loss
        for N in self.nodes:
            N.P = 0
        lhl=[l.getLoopHeadLoss(self.nodes) for l in self.loops]
        return lhl

    def setMinNodePressureHead(self, minPH, calcK=True):
        """
        Assuming we have calculated the pressure head at each node in the calculation of loop pressure drops,
        where we set ph at node 'a'=0, we can find the minimum sprinkler pressure and set it to minPH by adding
        delta to each node PH.  This is same as saying increase the pressure output of the pump at node 'a'
        :param minPH:
        :param calcK:
        :return:
        """
        minph = None
        for n in self.nodes:  #look for minimum sprinkler pressure head
            if n.IsSprinkler and (minph == None or n.P < minph):
                minph = n.P
        delta = minPH - minph  #calculate how much we need to up the pressure head at the pump
        #dial up the pump pressure/node pressure
        for n in self.nodes:
            n.P += delta
            if n.IsSprinkler and calcK:  #if called for, calculate discharge coefficient of sprinklers
                n.calc_k()

    def setRefPressure(self, RefP=0.0, nodeName='a'):
        deltaP=RefP-self.getNode(nodeName).P
        for n in self.nodes:
            n.P += deltaP

    def recalcPipeLengths(self, Riser=False):
        """
        When altering the elevation of the nodes in the pipe network, I make the assumption that
        the pipe connecting two nodes is a straight piece of pipe, now at an angle relative
        to the horizontal.  So calculate the new length of that pipe.
        The other option would be to assume a vertical riser and simply add on the length of the
        riser to the length of pipe.
        :return: length of the new pipe
        """
        for p in self.pipes:
            p.dZ=self.getNode(p.endNode).Z-self.getNode(p.startNode).Z
            p.length= p.length+p.dZ if Riser else math.sqrt(p.dZ**2+p.length**2)
    #endregion

    #region Functions to determine existance of Node, Pipe, or Loop by passing name
    def hasPipe(self, name):
        #determines if this pipe exists in the list of pipes
        return any(x.getName() == name for x in self.pipes)

    def hasLoop(self,name):
        #determines if this pipe exists in the list of pipes
        return any(x.getName() == name for x in self.loops)

    def hasNode(self, name):
        #determines if I have already constructed this node object (by name)
        return any(x.getName() == name for x in self.nodes)
    #endregion

    #region functions to get Node, Pipe, or Loop objects by passing name
    def getNode(self, name):
        return next(x for x in self.nodes if x.getName() == name)

    def getPipe(self, name):
        #returns a pipe object by its name
        return next(p for p in self.pipes if p.getName() == name)

    def getLoop(self, name):
        #returns a loop object by its name
        return next(l for l in self.loops if l.getName() == name)
    #endregion

    #region functions to get Node, Pipe or Loop index by passing name
    def getNodeIndex(self, name):
        """
        This way of searching a list of objects uses enumerate to return
        both the index in the list and the object from the list.  We can
        compare the object name property to the target name and return the
        index value for the object with the desired property.
        :param name:
        :return:
        """
        I=next(i for i, x in enumerate(self.nodes) if x.getName() == name)
        return I

    def getPipeIndex(self, name):
        """
        This way of searching a list of objects uses enumerate to return
        both the index in the list and the object from the list.  We can
        compare the object name property to the target name and return the
        index value for the object with the desired property.
        :param name:
        :return:
        """
        I=next(i for i, x in enumerate(self.pipes) if x.getName() == name)
        return I

    def getLoopIndex(self, name):
        """
        This way of searching a list of objects uses enumerate to return
        both the index in the list and the object from the list.  We can
        compare the object name property to the target name and return the
        index value for the object with the desired property.
        :param name:
        :return:
        """
        I=next(i for i, x in enumerate(self.loops) if x.getName() == name)
        return I
    #endregion

class PipeNetwork_controller():
    def __init__(self):
        self.model = PipeNetwork_model()
        self.view = PipeNetwork_view()

    # region Functions to add, delete or modify Nodes, Pipes or Loops by passing QTreeWidget item
    def addPipe(self, item, PN=None):
        if PN is None:
            PN=self.model
        if type(item) is qtw.QTreeWidgetItem:
            self.addPipe_TreeItem(item, PN)
        elif type(item) is Pipe:
            PN.pipes.append(item)

    def addPipe_TreeItem(self, item, PN=None):
        if PN is None:
            PN = self.model
        self.modifyPipe(item, 0, PN)

    def addLoop(self, item, PN=None):
        if PN is None:
            PN = self.model
        if type(item) is qtw.QTreeWidgetItem:
            self.addLoop_TreeItem(item, PN)
        elif type(item) is Loop:
            PN.loops.append(item)

    def addLoop_TreeItem(self, item, PN=None):
        if PN is None:
            PN = self.model
        self.modifyLoop(item, 0, PN)

    def deletePipe(self, item, PN=None):
        if PN is None:
            PN = self.model
        name = item.text()[0]
        if PN.hasPipe(name):
            PN.pipes.pop(PN.getPipeIndex(name))

    def deleteLoop(self, item, PN=None):
        if PN is None:
            PN = self.model
        name = item.text()[0]
        if PN.hasLoop(name):
            PN.loops.pop(PN.getLoopIndex(name))

    def deleteNode(self, item, PN=None):
        if PN is None:
            PN = self.model
        name = item.text()[0]
        if PN.hasNode(name):
            PN.nodes.pop(PN.getNodeIndex(name))

    def modifyPipe(self, item, col=0, PN=None):
        if PN is None:
            PN = self.model
        i = item
        cols = i.columnCount()
        p = [i.text(c) for c in range(cols)]
        name, length, diam, rough = p
        length = float(length) * (1 if PN.SI else 1.0 / 3.3)
        diam = float(diam) * (1000.0 if PN.SI else 12 * 3.3)
        rough = float(rough) * (1 if PN.SI else 1.0 / 3.3)
        nodes = name.split('-')
        stNode = nodes[0]
        enNode = nodes[1]
        TF = PN.hasPipe(name)
        if TF:
            p = PN.getPipe(name)
            p.length = float(length)
            p.diam = float(diam)
            p.modifyPipe(start=stNode, end=enNode, length=float(length), dia=float(diam), roughness=float(rough),
                         fluid=Fluid(), SI=self.SI)
        else:
            PN.pipes.append(Pipe(start=stNode, end=enNode, length=float(length), dia=float(diam), roughness=float(rough), fluid=Fluid(), SI=PN.SI))
        pass

    def modifyNode(self, item=None, col=0, PN=None, **kwargs):
        if PN is None:
            PN = self.model
        if type(item) is qtw.QTreeWidgetItem:
            self.modifyNode_TreeItem(item, col, PN)
        else:
            name=kwargs.get('name')
            TF = PN.hasNode(name)
            if TF:
                oldnode=PN.getNode(name)
                index=PN.getNodeIndex(name)
                if kwargs.get('makesprinkler'):
                    extFlow=kwargs.get('ext_flow')
                    z=kwargs.get('z')
                    minPH=kwargs.get('min_ph')
                    PN.nodes[index] = SprinklerHead(oldnode=oldnode, ext_flow=extFlow, z=z, min_ph=minPH)
                else:
                    extFlow=kwargs.get('ext_flow')
                    z=kwargs.get('z')
                    specPH=kwargs.get('specified_PH')
                    PN.nodes[index] = Node(oldnode=oldnode, ext_flow=extFlow, z=z,specifiedP=specPH)
            else:
                if kwargs.get('makesprinkler'):
                    extFlow = kwargs.get('ext_flow')
                    z = kwargs.get('z')
                    minPH = kwargs.get('min_ph')
                    PN.nodes.append(name=name, ext_flow=extFlow, z=z, min_ph=minPH)
                else:
                    extFlow = kwargs.get('ext_flow')
                    z = kwargs.get('z')
                    specPH = kwargs.get('specified_PH')
                    PN.nodes.append(Node(name=name, ext_flow=extFlow, z=z, specifiedP=specPH))

    def modifyNode_TreeItem(self, item, col=0, PN=None):
        if PN is None:
            PN = self.model
        i = item
        cols = i.columnCount()
        p = [i.text(c) for c in range(cols)]
        name, extflow, issprinkler, minph = p
        TF = PN.hasNode(name)
        if issprinkler.lower().strip() == 'true':
            N = SprinklerHead(name=name, ext_flow=float(extflow),
                              min_ph=float(minph) if not minph.lower().strip() == 'none' else None)
        else:
            N = Node(name=name, ext_flow=float(extflow),
                     specifiedP=float(minph) if not minph.lower().strip() == 'none' else None)
        N.Pipes = PN.getNodePipes(N.getName())
        if TF:
            i = PN.getNodeIndex(N.getName())
            PN.nodes[i] = N
        else:
            PN.nodes.append(N)

    def modifyLoop(self, item, col=0, PN=None):
        """
        Generally, the item passed is a tree with a trunk of a node name and leaves of pipe names
        :param item: a tree widget item
        :param col: the col that was modified
        :return: nothing
        """
        if PN is None:
            PN = self.model
        i = item
        loopname = item.text()[0]
        pipes = []
        for c in range(item.childCount()):
            pipes.append(PN.getPipe(item.child(i).text()[0]))
        L = Loop(Name=loopname, loopPipes=pipes)
        TF = PN.hasLoop(loopname)
        if TF:
            i = PN.getLoopIndex(loopname)
            PN.loops[i] = L
        else:
            PN.loops.append(L)
    # endregion

    #region Functions responsible for importing/building objects from a data file
    def importPipeNetwork(self, data, PN=None):
        """
        Build a pipe network from a data file
        :param data: a list of lines from a file.  Lines starting with # are comments to be ignored
        :return:
        """
        if PN is None:
            PN = self.model
        i=0
        PN.pipes.clear()
        PN.nodes.clear()
        PN.loops.clear()
        while i <len(data):
            L=''+data[i]
            L=L.lower().strip()
            #read in a pipe
            if L.find('#')==0:
                n=1
                pass
            elif L.find('pipe')>=0:
                #read in a pipe and increase i to end of pipe section
                i=self.importPipe(data, i,PN)
                pass
            elif L.find('loop')>=0:
                i=self.importLoop(data, i, PN)
                pass
            elif L.find('node')>=0:
                i=self.importNode(data,i,PN)
            elif L.find('fluid')>=0:
                i=self.importFluid(data, i, PN)
                pass
            elif L.find('network')>=0:
                i=self.importUnits(data, i, PN)
                pass
            i+=1
        PN.buildNodes()
        pass

    def importPipe(self, data, j, PN=None):
        """
        Read in a pipe from the data that comes from a file
        :param data: the list of strings from reading a file
        :param j: an index value for the list where the pipe definition starts
        :return: an index value for where the pipe definition ends
        """
        if PN is None:
            PN = self.model
        j+=1
        L=data[j].lower()
        P=Pipe(SI=PN.SI, fluid=PN.Fluid)
        while L.find('/pipe')==-1:
            cells=data[j].split(':')
            if cells[0].lower().find('nodes')>=0:
                nodes=cells[1].replace('(','').replace(')','').strip().split(',')
                P.startNode=nodes[0]
                P.endNode=nodes[1]
            elif cells[0].lower().find('length')>=0:
                P.setLength(float(cells[1].strip()),self.SI)
            elif cells[0].lower().find('dia')>=0:
                P.setDiam(float(cells[1].strip()), self.SI)
            elif cells[0].lower().find('rough')>=0:
                P.setRough(float(cells[1].strip()), self.SI)
            j+=1
            L=data[j].lower()
        if not PN.hasPipe(P.getName()):
            PN.pipes.append(P)
        return j

    def importLoop(self, data, j, PN=None):
        """
        Read in a loop from the data that comes from a file
        :param data: the list of strings from reading a file
        :param j: an index value for the list where the loop definition starts
        :return: an index value for where the loop definition ends
        """
        if PN is None:
            PN = self.model
        j += 1
        L = data[j].lower()
        Lp=Loop()
        while L.find('/loop') == -1:
            cells = data[j].split(':')
            if cells[0].lower().find('pipes') >= 0:
                cells=cells[1].split(',')
                for c in cells:
                    Lp.pipes.append(PN.getPipe(c.strip().replace("'",'')))
            elif cells[0].lower().find('name') >= 0:
                Lp.name=cells[1].strip()
            j += 1
            L = data[j].lower()
        if not PN.hasLoop(Lp.getName()):
            PN.loops.append(copy(Lp))
        return j

    def importFluid(self, data, j, PN=None):
        """
        Read in a fluid from the data that comes from a file
        :param data: the list of strings from reading a file
        :param j: an index value for the list where the fluid definition starts
        :return: an index value for where the fluid definition ends
        """
        if PN is None:
            PN=self.model
        j += 1
        L = data[j].lower()
        f = Fluid(SI=PN.SI)
        while L.find('/fluid') == -1:
            cells = data[j].split(':')
            if cells[0].lower().find('mu') >= 0:
                f.mu=float(cells[1].strip())
            elif cells[0].lower().find('rho') >= 0:
                f.rho=float(cells[1].strip())
            j += 1
            L = data[j].lower()
        PN.Fluid=f
        return j

    def importUnits(self, data, j, PN=None):
        """
        Read in a loop from the data that comes from a file
        :param data: the list of strings from reading a file
        :param j: an index value for the list where the loop definition starts
        :return: an index value for where the loop definition ends
        """
        if PN is None:
            PN=self.model
        j += 1
        L = data[j].lower()
        while L.find('/network') == -1:
            cells = data[j].split(':')
            if cells[0].lower().find('units') >= 0:
                self.SI=cells[1].strip().lower()=='si'
            j += 1
            L = data[j].lower()
        return j

    def importNode(self, data, j, PN=None):
        """
        Read in a node from the data that comes from a file
        :param data: the list of strings from reading a file
        :param j: an index value for the list where the node definition starts
        :return: an index value for where the node definition ends
        """
        if PN is None:
            PN=self.model
        j+=1
        L=data[j].lower()
        N=Node()
        isSprinkler = False #assume not a sprinkler at first
        while L.find('/node')==-1:
            cells=data[j].split(':')
            if cells[0].lower().find('name')>=0:
                N.Name=cells[1].strip()
            elif cells[0].lower().find('sprinkler')>=0:
                isSprinkler = cells[1].strip().lower()=='true'
            elif cells[0].lower().find('ext')>=0:
                N.E=float(cells[1].strip())
            elif cells[0].lower().find('minp')>=0:
                N.MinPH=float(cells[1].strip())
            elif cells[0].lower().find('specified')>=0:
                N.SpecifiedP=float(cells[1].strip())
            j+=1
            L=data[j].lower()
        if not PN.hasNode(N.Name):
            N.Pipes=PN.getNodePipes(N.Name)
            PN.nodes.append(N)
        else:
            PN.getNode(N.Name).modifyNode(ext_flow=N.E, z=N.Z, specifiedP=N.SpecifiedP, min_ph=N.MinPH)
        if isSprinkler:  #may need to convert a node to a sprinkler head
            if PN.getNode(N.Name).IsSprinkler==False:
                n=PN.getNodeIndex(N.Name)
                PN.nodes[n]=SprinklerHead(oldnode=PN.getNode(N.Name))
        return j
    #endregion

class PipeNetwork_view():
    def __init__(self):
        pass

    def setWidgets(self,args):
        self.tree_Pipes, self.tree_Nodes, self.tree_Loops, self.lbl_FlowRates, self.lbl_NodePressures, self.lbl_PipeHeadLosses=args
        pass

    #region Functions to display on GUI
    def updatePipeTree(self, tree, PN=None):
        if tree is None:
            tree=self.tree_Pipes
        tree.clear()
        for p in PN.pipes:
            itm = qtw.QTreeWidgetItem((p.getName(), str(p.length), str(p.diam), str(p.rough)))
            itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
            tree.addTopLevelItem(itm)

    def updateNodeTree(self, tree, PN=None):
        if tree is None:
            tree=self.tree_Nodes;
        tree.clear()
        for n in PN.nodes:
            itm = qtw.QTreeWidgetItem(
                (n.getName(), str(n.E), str(n.IsSprinkler), str(n.MinPH if n.IsSprinkler else n.SpecifiedP)))
            itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
            tree.addTopLevelItem(itm)

    def updateLoopTree(self, tree, PN=None):
        if tree is None:
            tree=self.tree_Loops
        tree.clear()
        for l in PN.loops:
            loopName = [l.getName()]
            pipes = l.getPipes()
            loop = qtw.QTreeWidgetItem(loopName)
            loop.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
            for p in pipes:
                itm = qtw.QTreeWidgetItem([p.getName()])
                itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
                loop.addChild(itm)
            tree.addTopLevelItem(loop)
    #endregion

    #region Functions to print output to CLI
    def printPipeFlowRates(self, SI=True, PN=None):
        for p in PN.pipes:
            p.printPipeFlowRate(SI=SI)

    def printNetNodeFlows(self, SI=True, PN=None):
        for n in PN.nodes:
            Q = n.getNetFlowRate()
            Q = Q if SI else Q / 1000 * (3.3 ** 3)
            q_units = 'L/s' if SI else 'cfs'
            print('Q into {} = {:0.2f} {}'.format(n.Name, n.getNetFlowRate(), q_units))

    def printLoopHeadLoss(self, SI=True, PN=None):
        for l in PN.loops:
            hl=l.getLoopHeadLoss( PN.nodes)
            hl=hl if SI else PN.pipes[0].fluid.m_to_psi(hl)
            p_units='m of water' if SI else 'psi'

            print('HL for {} = {:0.2f} {}'.format(l.name, l.getLoopHeadLoss(PN.nodes),p_units))

    def printPipeHeadLoss(self, SI=True, PN=None):
        for p in PN.pipes:
            hl = p.HeadLoss()
            hl = hl if SI else hl*12.0*3.3
            l = p.length if SI else p.length * 3.3
            d = p.diam if SI else p.diam * 3.3 * 12
            p_units = 'm of water' if SI else 'in of water'

            print('HL for {} (L={:0.2f}, d={:0.3f}) = {:0.2f} {}'.format(p.getName(), l, d, hl, p_units))

    def printNodeHead(self, SI=True, PN=None):
        for N in PN.nodes:
            p=N.P if SI else PN.pipes[0].fluid.m_to_psi(N.P)
            p_units='m of water' if SI else 'psi'
            if N.IsSprinkler:
                print('PH at sprinkler {} (Z={:0.2f}m) = {:0.2f} {} (K={:0.2f}, Q={:0.2f} L/S)'.format(N.Name, N.Z, p, p_units, N.K, abs(N.E)))
            else:
                print('PH at node {} (Z={:0.2f}m) = {:0.2f} {}'.format(N.Name, N.Z, p, p_units))
    #endregion

    # region Functions to output strings to be displayed on GUI
    def getPipeFlowRatesOutput(self, SI=True, tabs='', PN=None, lbl=None):
        stOut = 'Flow Results:\n'
        for p in PN.pipes:
            stOut += p.getPipeFlowRateOutput(SI=SI, tabs=tabs) + '\n'
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getNetNodeFlowsOutput(self, SI=True, tabs='', PN=None, lbl=None):
        stOut = 'Net Flow Into Nodes:\n'
        for n in PN.nodes:
            Q = n.getNetFlowRate()
            Q = Q if SI else Q / 1000 * (3.3 ** 3)
            q_units = 'L/s' if SI else 'cfs'
            stOut += tabs + 'Q into {} = {:0.2f} {}'.format(n.getName(), n.getNetFlowRate(), q_units) + '\n'
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getLoopHeadLossOutput(self, SI=True, tabs='', PN=None, lbl=None):
        stOut = 'Loop Head Losses:\n'
        for l in PN.loops:
            hl = l.getLoopHeadLoss( PN.nodes)
            hl = hl if SI else PN.pipes[0].fluid.m_to_psi(hl)
            p_units = 'm of water' if SI else 'psi'
            stOut += tabs + 'HL for {} = {:0.2f} {}'.format(l.getName(), l.getLoopHeadLoss( PN.nodes), p_units) + '\n'
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getPipeHeadLossesOutput(self, SI=True, tabs='', PN=None, lbl=None):
        stOut = 'Pipe Head Losses:\n'
        for p in PN.pipes:
            hl = p.HeadLoss()
            hl = hl if SI else p.fluid.m_to_psi(hl)
            l = p.length if SI else p.length * 3.3
            d = p.diam if SI else p.diam * 3.3 * 12
            p_units = 'm of water' if SI else 'psi'
            stOut += tabs + 'HL in {} (L={:0.2f}, d={:0.3f}) is {:0.2f} {}'.format(p.getName(), l, d, hl, p_units) + '\n'
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getNodeHeadOutput(self, SI=True, tabs='', PN=None, lbl=None):
        stOut = 'Node Pressures\n'
        for N in PN.nodes:
            p = N.P if SI else PN.pipes[0].fluid.m_to_psi(N.P)
            p_units = 'm of water' if SI else 'psi'
            if N.IsSprinkler:
                stOut += tabs + 'PH at sprinkler {} (Z={:0.2f}m) = {:0.2f} {} (K={:0.2f}, Q={:0.2f} L/S)'.format(
                    N.getName(), N.Z, p, p_units, N.K, abs(N.E)) + '\n'
            else:
                stOut += tabs + 'PH at node {} (Z={:0.2f}m) = {:0.2f} {}'.format(N.getName(), N.Z, p, p_units) + '\n'
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getRealityCheckOutput(self, SI=True, tabs='', PN=None, lbl=None):
        lblA=qtw.QLabel()
        lblB=qtw.QLabel()
        self.getLoopHeadLossOutput(SI=SI,tabs='\t', PN=PN, lbl=lblA)
        stOut=lblA.text()
        stOut += '\n\n'
        self.getNetNodeFlowsOutput(SI=SI,tabs='\t', PN=PN, lbl=lblB)
        stOut+=lblB.text()
        if not lbl is None:
            lbl.setText(stOut)
        return stOut

    def getPipeNetworkOutputForFile(self, PN=None):
        stTmp = ''
        stTmp += '<Fluid>\n'
        stTmp += '\tmu: {:0.5f}\n'.format( PN.Fluid.mu)
        stTmp += '\trho: {:0.1f}\n'.format( PN.Fluid.rho)
        stTmp += '</Fluid>\n'
        for p in PN.pipes:
            stTmp += '<Pipe>\n'
            stTmp += '\tnodes: ({},{})\n'.format(p.startNode, p.endNode)
            stTmp += '\tlength: {:0.1f}\n'.format(p.length)
            stTmp += '\tdiam: {:0.1f}\n'.format(p.diam * 1000.0)
            stTmp += '\trough: {:0.5f}\n'.format(p.rough)
            stTmp += '</Pipe>\n'
        for l in PN.loops:
            stTmp += '<Loop>\n'
            stTmp += '\tName: {}\n'.format(l.getName())
            stTmp += '\tPipes: '
            for i in range(len(l.pipes)):
                p = l.pipes[i]
                stTmp += (', ' if i > 0 else '') + '{}'.format(p.getName())
            stTmp += '\n'
            stTmp += '</Loop>\n'
        for n in PN.nodes:
            if n.E > 0.0 or n.E < 0.0:
                stTmp += '<Node>\n'
                stTmp += '\tName: {}\n'.format(n.getName())
                stTmp += '\tExternal Flow: {:0.1f}\n'.format(n.E)
                stTmp += '\tMinP: {}\n'.format(n.MinPH) if n.IsSprinkler else ''
                stTmp += '\tSprinkler: True\n' if n.IsSprinkler else ''
                stTmp += '</Node>\n'
        return stTmp
    # endregion
#endregion

#region function definitions
def mainX2():
    """
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc.
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    """
    #instantiate a Fluid object to define the working fluid as water
    #water=Fluid(mu=0.00089,rho=1000.0)
    SIUnits=False
    water=Fluid(mu=20.50E-6, rho=62.3, SI=SIUnits)
    #roughness = 0.00025  # in m
    r_CI=0.00085 #ft
    r_CN=0.003 #ft

    #instantiate a new PipeNetwork object
    PN=PipeNetwork_model() # this is the model
    PNC=PipeNetwork_controller() # this is the controller that modifies the model and reads the view
    PNV=PipeNetwork_view()  # this is the view
    #add Pipe objects to the pipe network
    PNC.addPipe(Pipe('a','b',1000, 18, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('a','h',1600, 24, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('b','c',500, 18, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('b','e',800, 16, r_CI, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('c','d',500, 18, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('c','f',800, 16, r_CI, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('d','g',800, 16, r_CI, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('e','f',500, 12, r_CI, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('e','i',800, 18, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('f','g',500, 12, r_CI, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('g','j',800, 18, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('h','i',1000, 24, r_CN, water, SI=SIUnits),PN)
    PNC.addPipe(Pipe('i','j',1000, 24, r_CN, water, SI=SIUnits),PN)

    #add Node objects to the pipe network
    PN.buildNodes()
    #update the external flow of inlet node
    PNC.modifyNode(PN=PN, name='h', ext_flow=10)
    PNC.modifyNode(PN=PN, name='e', ext_flow=-3)
    PNC.modifyNode(PN=PN, name='f', ext_flow=-5)
    PNC.modifyNode(PN=PN, name='h', ext_flow=-2)
    PN.getNode('h').setExtFlow(10,SI=SIUnits)
    #set external flows at required nodes
    PN.getNode('e').setExtFlow(-3, SI=SIUnits)
    PN.getNode('f').setExtFlow(-5, SI=SIUnits)
    PN.getNode('d').setExtFlow(-2, SI=SIUnits)

    #add Loop objects to the pipe network
    PNC.addLoop(Loop('A',[PN.getPipe('a-b'),PN.getPipe('b-e'), PN.getPipe('e-i'),PN.getPipe('h-i'), PN.getPipe('a-h')]), PN=PN)
    PNC.addLoop(Loop('B',[PN.getPipe('b-c'), PN.getPipe('c-f'),PN.getPipe('e-f'), PN.getPipe('b-e')]), PN=PN)
    PNC.addLoop(Loop('C',[PN.getPipe('c-d'), PN.getPipe('d-g'),PN.getPipe('f-g'), PN.getPipe('c-f')]), PN=PN)
    PNC.addLoop(Loop('D',[PN.getPipe('e-f'),PN.getPipe('f-g'), PN.getPipe('g-j'),PN.getPipe('i-j'), PN.getPipe('e-i')]), PN=PN)

    #call the findFlowRates method of the PN (a PipeNetwork object)
    PN.findFlowRates()
    pRef=water.psi_to_m(80) #reference pressure in m of water
    PN.setRefPressure(pRef,'h')  #set reference pressure at node 'h' and increase all other nodes too.

    #get output for zero height nodes
    SIUnits=False
    PNV.printPipeFlowRates(SI=SIUnits, PN=PN)
    print('')
    print('Check node flows:')
    PNV.printNetNodeFlows(SI=SIUnits, PN=PN)
    print('')
    print('Check loop head loss:')
    PNV.printLoopHeadLoss(SI=SIUnits, PN=PN)
    print('')
    PNV.printPipeHeadLoss(SI=SIUnits, PN=PN)
    print('')
    PNV.printNodeHead(SI=SIUnits, PN=PN)

    #new node elevations
    #now change node elevations (which also changes pipe lengths by geometry assuming straight pipe runs between nodes)
    # sf=1.0
    # PN.getNode('a').Z= 0.0
    # PN.getNode('b').Z= 2.5/sf #sprinkler
    # PN.getNode('c').Z= 5.0/sf
    # PN.getNode('d').Z= 5.0/sf #sprinkler
    # PN.getNode('e').Z= 5.0/sf
    # PN.getNode('f').Z= 4.0/sf #sprinkler
    # PN.getNode('g').Z= 4.5/sf
    # PN.getNode('h').Z= 4.0/sf #sprinkler
    #
    # PN.recalcPipeLengths()
    #
    # #calculate the flow rates for the new pipe network layout without fixed sprinkler flow rates
    # PN.findFlowRates2()
    #
    # print()
    # print('Second pipe network')
    # PN.printPipeFlowRates()
    # print('')
    # print('Check node flows:')
    # PN.printNetNodeFlows()
    # print('')
    # print('Check loop head loss:')
    # PN.printLoopHeadLoss()
    # print('')
    # PN.printPipeHeadLoss()
    # print('')
    # PN.printNodeHead()
    # print('Mass Balance Network: {:0.3f}'.format(PN.getNetOutFlow(setmin=False)[0]))

def mainHW6():
    '''
    This program analyzes flows in a given pipe network based on the following:
    1. The pipe segments are named by their endpoint node names:  e.g., a-b, b-e, etc.
    2. Flow from the lower letter to the higher letter of a pipe is considered positive.
    3. Pressure decreases in the direction of flow through a pipe.
    4. At each node in the pipe network, mass is conserved.
    5. For any loop in the pipe network, the pressure loss is zero
    Approach to analyzing the pipe network:
    Step 1: build a pipe network object that contains pipe, node, loop and fluid objects
    Step 2: calculate the flow rates in each pipe using fsolve
    Step 3: output results
    Step 4: check results against expected properties of zero head loss around a loop and mass conservation at nodes.
    :return:
    '''
    #instantiate a Fluid object to define the working fluid as water
    water=Fluid(mu=0.00089,rho=1000.0)
    roughness = 0.00025  # in m

    #instantiate a new PipeNetwork object
    PN = PipeNetwork_model()
    #instantiate a new PipeNetworkController
    PNC = PipeNetwork_controller()
    #instantiate a new PipeNetworkView
    PNV = PipeNetwork_view()
    #add Pipe objects to the pipe network through PNC
    PNC.addPipe(Pipe('a','b',185.0, 300.0, roughness, water),PN)
    PNC.addPipe(Pipe('a','c',100.0, 200.0, roughness, water),PN)
    PNC.addPipe(Pipe('b','e',165.0, 200.0, roughness, water),PN)
    PNC.addPipe(Pipe('c','d',125.0, 200.0, roughness, water),PN)
    PNC.addPipe(Pipe('c','f',100.0, 150.0, roughness, water),PN)
    PNC.addPipe(Pipe('d','e',125.0, 200.0, roughness, water),PN)
    PNC.addPipe(Pipe('d','g',100.0, 150.0, roughness, water),PN)
    PNC.addPipe(Pipe('e','h',100.0, 150.0, roughness, water),PN)
    PNC.addPipe(Pipe('f','g',125.0, 250.0, roughness, water),PN)
    PNC.addPipe(Pipe('g','h',125.0, 250.0, roughness, water),PN)
    #add Node objects to the pipe network
    PN.buildNodes()
    #update the external flow of inlet node
    PNC.modifyNode(item=None, PN=PN, name='a', ext_flow=60)
    #PN.getNode('a').E=float(input("External flow rate at node 'a'?"))
    #place the sprinkler heads at required nodes
    PNC.modifyNode(item=None, PN=PN,name='b', makesprinkler=True, ext_flow=-10, z=0, min_ph=2.0)
    PNC.modifyNode(item=None, PN=PN,name='d', makesprinkler=True, ext_flow=-20, z=0, min_ph=2.0)
    PNC.modifyNode(item=None, PN=PN,name='f', makesprinkler=True, ext_flow=-15, z=0, min_ph=2.0)
    PNC.modifyNode(item=None, PN=PN,name='h', makesprinkler=True, ext_flow=-15, z=0, min_ph=2.0)

    #add Loop objects to the pipe network
    PNC.addLoop(Loop('A',[PN.getPipe('a-b'), PN.getPipe('b-e'),PN.getPipe('d-e'), PN.getPipe('c-d'), PN.getPipe('a-c')]), PN=PN)
    PNC.addLoop(Loop('B',[PN.getPipe('c-d'), PN.getPipe('d-g'),PN.getPipe('f-g'), PN.getPipe('c-f')]), PN=PN)
    PNC.addLoop(Loop('C',[PN.getPipe('d-e'), PN.getPipe('e-h'),PN.getPipe('g-h'), PN.getPipe('d-g')]), PN=PN)

    #call the findFlowRates method of the PN (a PipeNetwork object)
    PN.findFlowRates()
    PN.setMinNodePressureHead(2)

    #get output for zero height nodes
    PNV.printPipeFlowRates(PN=PN)
    print('')
    print('Check node flows:')
    PNV.printNetNodeFlows(PN=PN)
    print('')
    print('Check loop head loss:')
    PNV.printLoopHeadLoss(PN=PN)
    print('')
    PNV.printPipeHeadLoss(PN=PN)
    print('')
    PNV.printNodeHead(PN=PN)

    #new node elevations
    #now change node elevations (which also changes pipe lengths by geometry assuming straight pipe runs between nodes)
    sf=1.0
    PNC.modifyNode(PN=PN, name='a', z=0.0)
    PNC.modifyNode(PN=PN, name='b', z=2.5, makesprinkler=True)
    PNC.modifyNode(PN=PN, name='c', z=5.0)
    PNC.modifyNode(PN=PN, name='d', z=5.0, makesprinkler=True)
    PNC.modifyNode(PN=PN, name='e', z=5.0)
    PNC.modifyNode(PN=PN, name='f', z=4.0, makesprinkler=True)
    PNC.modifyNode(PN=PN, name='g', z=4.5)
    PNC.modifyNode(PN=PN, name='h', z=4.0, makesprinkler=True)
    # PN.getNode('a').Z= 0.0
    # PN.getNode('b').Z= 2.5/sf #sprinkler
    # PN.getNode('c').Z= 5.0/sf
    # PN.getNode('d').Z= 5.0/sf #sprinkler
    # PN.getNode('e').Z= 5.0/sf
    # PN.getNode('f').Z= 4.0/sf #sprinkler
    # PN.getNode('g').Z= 4.5/sf
    # PN.getNode('h').Z= 4.0/sf #sprinkler

    PN.recalcPipeLengths()

    #calculate the flow rates for the new pipe network layout without fixed sprinkler flow rates
    PN.findFlowRates2()

    print()
    print('Second pipe network')
    PNV.printPipeFlowRates(PN=PN)
    print('')
    print('Check node flows:')
    PNV.printNetNodeFlows(PN=PN)
    print('')
    print('Check loop head loss:')
    PNV.printLoopHeadLoss(PN=PN)
    print('')
    PNV.printPipeHeadLoss(PN=PN)
    print('')
    PNV.printNodeHead(PN=PN)
    print('Mass Balance Network: {:0.3f}'.format(PN.getNetOutFlow(setmin=False)[0]))#
#endregion

if __name__=='__main__':
    #mainHW6()
    mainX2()
