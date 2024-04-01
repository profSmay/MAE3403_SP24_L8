#region imports
from PipeNetwork_GUI import Ui_Form
import PyQt5.QtWidgets as qtw
import PyQt5.QtCore as qtc
import PyQt5.QtGui as qtg
import sys
from PipeNetwork_classes import PipeNetwork_controller
#endregion

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        Main window constructor.
        """
        super().__init__()
        self.setupUi(self)

        self.Controller = PipeNetwork_controller()  # controller for modifying the Pipe Network #$CONTROLLER$#
        self.Model = self.Controller.model  # start with an empty pipe network #$MODEL$#
        self.View = self.Controller.view  # view for Pipe Network #$VIEW$#
        self.View.setWidgets((self.tree_Pipes,self.tree_Nodes, self.tree_Loops, self.lbl_FlowRates, self.lbl_NodePressures, self.lbl_PipeHeadLosses))
        # Main UI code goes here
        self.setupSignalsSlotsEventFilter()
        # end Main UI code
        self.show()  # show the main widget

    def setupSignalsSlotsEventFilter(self):
        """
        This is called from the constructor to wire up the signals and slots and install event filter if needed.
        :return: nothing
        """
        # CONNECTING SIGNALS AND SLOTS
        # these available signals can be easily connected to a slot
        self.btn_AddPipe.clicked.connect(self.addPipeToPipeList)
        self.btn_DeletePipe.clicked.connect(self.deletePipe)
        self.btn_AddPipeToLoop.clicked.connect(self.addPipeToLoop)
        self.btn_CreateLoop.clicked.connect(self.createLoop)
        self.btn_Evaluate.clicked.connect(self.EvaluatePipeNetwork)
        self.tree_Pipes.itemChanged.connect(self.modifyPipe)
        self.tree_Nodes.itemChanged.connect(self.modifyNode)
        self.tree_Loops.itemChanged.connect(self.modifyLoop)
        self.btn_OpenPipeNetworkFile.clicked.connect(self.readPipeNetworkFile)

        # INSTALLING AN EVENT FILTER ON THE TREE WIDGETS
        # I do this if there is no easy signal to connect to do what I want.
        # The event filter allows the widget to analyze events from the windows event loop when
        # they occur on the widget.
        # if a signal is not available, I can use the event filter to take action
        self.tree_Pipes.installEventFilter(self)
        self.tree_LoopPipes.installEventFilter(self)
        self.tree_Loops.installEventFilter(self)

    def eventFilter(self, obj, event):
        """
        This overrides the default eventFilter of the widget.  It takes action on events and then passes the event
        along to the parent widget.
        :param obj: The object on which the event happened
        :param event: The event itself
        :return: boolean from the parent widget
        """
        # allow tree_pipes, tree_nodes, tree_LoopPipes, and tree_Loops to respond to delete key
        if event.type() == qtc.QEvent.KeyPress:
            if event.key() == qtc.Qt.Key_Delete:
                if obj == self.tree_Pipes:
                    self.deletePipe()
                elif obj == self.tree_Nodes:
                    self.deleteNode()
                elif obj == self.tree_LoopPipes:
                    self.deleteLoopPipe()
                elif obj == self.tree_Loops:
                    self.deleteLoop()

        # pass the event along to the parent widget if there is one.
        return super(MainWindow, self).eventFilter(obj, event)

    #region Functions that act as Slots
    def readPipeNetworkFile(self):
        """
        Read the information from a pipe network file.
        :return:
        """
        # open the file dialog box to search for the file I want
        filename = qtw.QFileDialog.getOpenFileName()[0]
        if len(filename) == 0:  # no file selected
            return
        self.le_FileName.setText(filename)  # echo the filename on the GUI
        file = open(filename, 'r')  # open the file
        data = file.readlines()  # read all the lines of the file into a list of strings
        self.Controller.importPipeNetwork(data, PN=self.Model)  # import the pipe network information
        self.updateView()  # update the view of the model
        pass

    def updateNodesTree(self):
        nodesInUse = []
        for i in range(self.tree_Pipes.topLevelItemCount()):
            pipeName = self.tree_Pipes.topLevelItem(i).text(0)
            pipeName = pipeName.split('-')
            stNode = pipeName[0]
            enNode = pipeName[1]
            if nodesInUse.__contains__(stNode) is False:
                nodesInUse.append(stNode)
            if nodesInUse.__contains__(enNode) is False:
                nodesInUse.append(enNode)
        createdNodes = []
        for i in range(self.tree_Nodes.topLevelItemCount()):
            createdNodes.append(self.tree_Nodes.topLevelItem(i).text(0))
        # now, create a new node for any nodeInUse that is not in createdNodes
        for n in nodesInUse:
            if createdNodes.__contains__(n) == False:
                itm = qtw.QTreeWidgetItem([n, '0', 'False'])
                itm.setFlags(
                    qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
                self.tree_Nodes.addTopLevelItem(itm)

    def addPipeToPipeList(self):
        """
        I use this as the slot for the clicked signal of the Add Pipe button.  It reads from the
        line edit boxes and creates a top level QTreeWidgetItem and places it in the tree_Pipes widget.
        :return: none
        """
        node1 = self.le_StartNode.text()  # read from GUI
        node2 = self.le_EndNodeName.text()  # read from GUI
        name = '{}-{}'.format(min(node1, node2), max(node1, node2))  # follow alphabetical naming convention
        length = self.le_PipeLength.text()  # read from GUI
        diam = self.le_Diam.text()  # read from GUI
        rough = self.le_Roughness.text()  # read from GUI
        itm = qtw.QTreeWidgetItem((name, length, diam, rough))
        self.Controller.addPipe(itm, PN=self.Model)  # part of the controller that updates the model with changes on the view
        itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
        self.tree_Pipes.addTopLevelItem(itm)
        self.updateNodesTree()

    def modifyPipe(self, item, col):
        self.Controller.modifyPipe(item, col, PN=self.Model)

    def modifyNode(self, item, col):
        self.Controller.modifyNode(item, col, PN=self.Model)

    def modifyLoop(self, item, col):
        self.Controller.modifyLoop(item, col, PN=self.Model)

    def deletePipe(self):
        """
        This is the slot for the clicked signal of the Delete Pipe button as well as the response to the delete
        key being pressed when in the tree_Pipe widget.  This latter behavior is implemented with the eventFilter
        method that is installed on the tree_Pipe widget.
        :return: none
        """
        index = self.tree_Pipes.currentIndex().row()
        self.Controller.deletePipe(self.tree_Pipes.takeTopLevelItem(index), PN=self.Model)

    def deleteNode(self):
        index = self.tree_Nodes.currentIndex().row()
        self.Controller.deleteNode(self.tree_Nodes.takeTopLevelItem(index), PN=self.Model)

    def deleteLoopPipe(self):
        """
        This is the response to the delete
        key being pressed when in the tree_LoopPipes widget.  It implemented with the eventFilter
        method that is installed on the tree_LoopPipe widget.
        :return:
        """
        index = self.tree_LoopPipes.currentIndex().row()
        self.tree_LoopPipes.takeTopLevelItem(index)

    def deleteLoop(self):
        """
        This is the response to the delete
        key being pressed when in the tree_Loops widget.  It implemented with the eventFilter
        method that is installed on the tree_Loops widget.
        :return:
        """
        isParent = self.tree_Loops.currentItem().parent() is None
        itm = self.tree_Loops.currentItem()
        index = self.tree_Loops.currentIndex().row()
        if isParent:
            self.Controller.deleteLoop(self.tree_Loops.takeTopLevelItem(index), PN=self.Model)
        else:
            parent = self.tree_Loops.currentItem().parent()
            self.tree_Loops.currentItem().parent().removeChild(itm)
            self.PND.modifyLoop(parent, 0, PN=self.Model)

    def addPipeToLoop(self):
        """
        I use this as the slot for the clicked signal of the Add Pipe to loop button.  It reads from the
        tree_Pipes widget and creates a top level QTreeWidgetItem and places it in the tree_LoopPipes widget.
        :return: none
        """
        for p in self.tree_Pipes.selectedItems():
            name = p.text(0)
            rows = str(self.tree_LoopPipes.topLevelItemCount())
            itm = qtw.QTreeWidgetItem((rows, name))
            itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
            self.tree_LoopPipes.addTopLevelItem(itm)
        rows = self.tree_LoopPipes.topLevelItemCount()

    def createLoop(self):
        """
        This is the slot for the clicked signal of btn_CreateLoop.  It reads from the tree_LoopPipes and builds
        the hierarchy of a loop object and adds it to the tree_Loops widget.
        :return:
        """
        loopName = [self.le_LoopName.text()]
        pipes = []
        while self.tree_LoopPipes.topLevelItemCount() > 0:
            pipe = self.tree_LoopPipes.takeTopLevelItem(0)
            pipes.append([pipe.text(1)])
        loop = qtw.QTreeWidgetItem(loopName)
        loop.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
        for p in pipes:
            itm = qtw.QTreeWidgetItem(p)
            itm.setFlags(qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable | qtc.Qt.ItemIsDragEnabled | qtc.Qt.ItemIsUserCheckable | qtc.Qt.ItemIsEnabled)
            loop.addChild(itm)
        self.tree_Loops.addTopLevelItem(loop)
        self.Controller.addLoop(loop, PN=self.Model)

    def EvaluatePipeNetwork(self):
        """
        Use the Pipe Network object to evaluate the pipe network and display output
        :return:
        """
        self.Model.findFlowRates()
        self.Model.setMinNodePressureHead(minPH=2, calcK=True)
        # show output
        self.View.getPipeFlowRatesOutput(tabs='\t', PN=self.Model, lbl=self.lbl_FlowRates)
        self.View.getPipeHeadLossesOutput(tabs='\t', PN=self.Model, lbl=self.lbl_PipeHeadLosses)
        self.View.getNodeHeadOutput(tabs='\t', PN=self.Model, lbl=self.lbl_NodePressures)
        self.View.getRealityCheckOutput(tabs='\t', PN=self.Model, lbl=self.lbl_PressureAndFlowChecks)
        return None
    #endregion

    def updateView(self):
        """
        Update the representation of the model in the tree Widgets
        :return:
        """
        self.tree_LoopPipes.clear()
        # update the pipe tree
        self.View.updatePipeTree(self.tree_Pipes, PN=self.Model)
        #update the node tree
        self.View.updateNodeTree(self.tree_Nodes, PN=self.Model)
        #update the loop tree
        self.View.updateLoopTree(self.tree_Loops, PN=self.Model)

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Pipe Network Designer')
    sys.exit(app.exec())
