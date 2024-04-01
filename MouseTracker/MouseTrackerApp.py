from MouseTracker import Ui_MouseTracker
from PyQt5 import uic
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg

class MainWindow(qtw.QWidget, Ui_MouseTracker):
    def __init__(self):
        '''
        Main window constructor.
        '''
        super().__init__()
        self.setupUi(self)
        # Main UI code goes here
        #I'm installing an event filter because there is no signal for keypress
        self.spin_MouseTracker.installEventFilter(self)
        self.spin_MouseTracker.setMouseTracking(True)
        self.pushButton.pressed.connect(self.colorButton)
        self.pushButton.installEventFilter(self)
        self.pushButton_2.installEventFilter(self)
        self.pushButton.setMouseTracking(True)
        self.pushButton_2.setMouseTracking(True)
        # end Main UI code
        self.show()  #show the main widget

    #since my main window is a widget, I can customize its events by overriding the default event
    def mouseMoveEvent(self, event):
        self.setWindowTitle('x:{}, y:{}'.format(event.x(), event.y()))

    def colorButton(self):
        self.pushButton.setAutoFillBackground(False)

    #the event filter allows me to respond to specific events for a specific object
    def eventFilter(self, obj, event):
        if event.type()==qtc.QEvent.KeyPress:
            et=event.type()
            if event.key()==qtc.Qt.Key_Up:
                self.setWindowTitle('up pushed')
            if event.key() == qtc.Qt.Key_Down:
                self.setWindowTitle('down pushed')
            if event.key() == qtc.Qt.Key_O:
                self.setWindowTitle('Go Pokes!!')
            pass
        if event.type() == event.Enter:
            if obj == self.pushButton:
                self.pushButton.setText("Over me")
            if obj == self.pushButton_2:
                self.pushButton_2.setText("Over me")
        if event.type() == event.Leave:
            if obj == self.pushButton:
                self.pushButton.setText("pushButton")
            if obj == self.pushButton_2:
                self.pushButton_2.setText("pushButton2")
        return super(MainWindow, self).eventFilter(obj, event)

if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Mouse Tracker')
    sys.exit(app.exec())