
import numpy as np
import matplotlib.pyplot as plt

import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog,QMessageBox
from PyQt5.QtGui import QCursor
from PyQt5.QtCore import Qt

from Aero_Demo_ui import Ui_Dialog
from Wing_Class import Wing

class main_window(QDialog, Ui_Dialog):
    def __init__(self):
        '''
        Constructor for main_window.  This is a little different than inheriting from
        QWidgit, though QDialog inherits from QWidgit.
        '''
        super(main_window,self).__init__()
        self.setupUi(self)
        #signals and slots are assigned in this function
        self.assign_widgets()

        self.wing=None  # the primary data element for this program

        self.show()

    def assign_widgets(self):
        '''
        connect signals and slots
        :return:
        '''
        self.pushButton_Exit.clicked.connect(self.ExitApp)
        self.pushButton_GetWing.clicked.connect(self.GetWing)
        self.comboBox.currentTextChanged.connect(self.combobox_changed)
        #add some choices to the comboBox
        self.comboBox.addItem("Birch")
        self.comboBox.addItem("Oak")
        self.comboBox.addItem("Sycamore")

    def combobox_changed(self):
        '''
        Slot for comboBox currentTextChanged
        :return:
        '''
        text = self.comboBox.currentText()
        self.lineEdit_fromCombo.setText(text)
        pass

    def GetWing(self):
        # get the filename using the OPEN dialog
        filename=QFileDialog.getOpenFileName()[0]
        if len(filename)==0: 
            no_file()
            return
        self.textEdit_filename.setText(filename)
        #we do this in case it takes a long time to read the file
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))

        # Read the file
        f1 = open(filename, 'r')  # open the file for reading
        data = f1.readlines()  # read the entire file as a list of strings
        f1.close()  # close the file  ... very important

        self.wing=Wing()  # create a wing instance (object)

        try:  #an example of handling an error
            self.wing.processWingData(data)
            self.lineEdit_Height.setText('{:.2f}'.format(self.wing.sparH))
            self.lineEdit_Width.setText('{:.2f}'.format(self.wing.sparW))
            report = self.wing.generate_report()
            self.plainTextEdit_Report.setPlainText(report)

            QApplication.restoreOverrideCursor()
        except:
            QApplication.restoreOverrideCursor()
            bad_file()

    def PlotSomething(self):
        x=np.linspace(0,6*np.pi,300)
        y=np.zeros_like(x)
        for i in range(300):
            y[i]=np.exp(-x[i]/5)*np.sin(x[i])
        plt.plot(x,y)
        plt.show()
        return

    def ExitApp(self):
        app.exit()

def no_file():
    msg = QMessageBox()
    msg.setText('There was no file selected')
    msg.setWindowTitle("No File")
    retval = msg.exec_()
    return None

def bad_file():
    msg = QMessageBox()
    msg.setText('Unable to process the selected file')
    msg.setWindowTitle("Bad File")
    retval = msg.exec_()
    return

if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())
    
 





