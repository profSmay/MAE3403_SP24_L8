# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MouseTracker.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MouseTracker(object):
    def setupUi(self, MouseTracker):
        MouseTracker.setObjectName("MouseTracker")
        MouseTracker.resize(710, 698)
        MouseTracker.setMouseTracking(True)
        self.spin_MouseTracker = QtWidgets.QSpinBox(MouseTracker)
        self.spin_MouseTracker.setGeometry(QtCore.QRect(270, 250, 91, 25))
        self.spin_MouseTracker.setMouseTracking(True)
        self.spin_MouseTracker.setToolTipDuration(3000)
        self.spin_MouseTracker.setObjectName("spin_MouseTracker")
        self.pushButton = QtWidgets.QPushButton(MouseTracker)
        self.pushButton.setGeometry(QtCore.QRect(260, 80, 112, 34))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(MouseTracker)
        self.pushButton_2.setGeometry(QtCore.QRect(50, 450, 112, 34))
        self.pushButton_2.setObjectName("pushButton_2")

        self.retranslateUi(MouseTracker)
        QtCore.QMetaObject.connectSlotsByName(MouseTracker)

    def retranslateUi(self, MouseTracker):
        _translate = QtCore.QCoreApplication.translate
        MouseTracker.setWindowTitle(_translate("MouseTracker", "Form"))
        self.spin_MouseTracker.setToolTip(_translate("MouseTracker", "This object has an event filter installed"))
        self.pushButton.setText(_translate("MouseTracker", "PushButton"))
        self.pushButton_2.setText(_translate("MouseTracker", "PushButton"))
