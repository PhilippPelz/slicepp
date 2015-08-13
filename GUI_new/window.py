
#https://github.com/MapQuest/avecado/blob/master/src/python_module.cpp
#http://renatogarcia.blog.br/en/posts/boost-property-tree-extensions.html

# https://docs.python.org/2/library/json.html
# http://www.boost.org/doc/libs/1_46_1/doc/html/boost_propertytree/parsers.html#boost_propertytree.parsers.info_parser
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.Qt as Qt
import numpy as np
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
from PyQt4.Qt import QTreeWidgetItem
#from pyE17 import io
from pyE17.io.h5rw import h5read, h5write
import os
import json

projectPath = os.path.dirname(os.path.realpath(__file__)) + '/../'
projectPath1 = projectPath + '/validation/si3n4/'
defaultResultsPath= projectPath1 + '7.11.h5'
defaultLoadConfigPath= projectPath1 + '7.11.json'
defaultSaveConfigPath= projectPath1 + '7.11.json'
defaultExePath= projectPath + 'build/release/bin/stem3'

class SliceGUI(QtGui.QMainWindow):
    def __init__(self):
        super(SliceGUI, self).__init__()
        
        self.labelText = [['Abs', 'Phase'], ['Real', 'Imag']]

        self.resultsRadio1 = QtGui.QRadioButton()
        self.resultsRadio2 = QtGui.QRadioButton()
        self.view1 = pg.ImageView(view=pg.PlotItem())
        self.view2 = pg.ImageView(view=pg.PlotItem())
        self.trResults = pg.DataTreeWidget()
        self.v1Label = QtGui.QLabel(self.labelText[0][0])
        self.v2Label = QtGui.QLabel(self.labelText[0][1])
        self.resultssplitter = QtGui.QSplitter(1)
        self.consolelayout = QtGui.QVBoxLayout()     
        
        self.btnLoadConfigFile = QtGui.QPushButton('Load')
        self.btnSaveConfigFile = QtGui.QPushButton('Save')
        self.btnLoadResults = QtGui.QPushButton('Load')
        self.loadFileButton = QtGui.QPushButton('Load Configuration File')   
        
        self.process = QtCore.QProcess(self)
        self.consoleOutput = QtGui.QTextEdit()  
        
        self.selected = np.array([1])
        
        self.addTopLevel()
        self.addLeftSide()
        self.addWaveItems()
        self.addResultItems()
        self.addConsoleItems()



    def addWaveItems(self):
        pass
    def addResultItems(self):
        d = {
            'list1': [1, 2, 3, 4, 5, 6, {'nested1': 'aaaaa', 'nested2': 'bbbbb'}, "seven"],
            'dict1': {
                'x': 1,
                'y': 2,
                'z': 'three'
            },
            'array1 (20x20)': np.ones((2, 2))
        }
        resultsViewLayout = QtGui.QVBoxLayout()
        rblayout = QtGui.QHBoxLayout()
        resultsViewWidget = QtGui.QWidget()
        resultsViewWidget.setLayout(resultsViewLayout)
        gb = QtGui.QGroupBox()

        self.resultsRadio1.setText("Abs/Phase")
        self.resultsRadio1.setChecked(True)
        self.resultsRadio2.setText("Real/Imag")

        gb.setLayout(rblayout)
        rblayout.addWidget(self.resultsRadio1)
        rblayout.addWidget(self.resultsRadio2)
        self.resultsRadio1.clicked.connect(self.updateResultsView)
        self.resultsRadio2.clicked.connect(self.updateResultsView)

        self.trResults.setMaximumWidth(800)
        self.trResults.connect(self.trResults, QtCore.SIGNAL('itemClicked(QTreeWidgetItem*, int)'), self.resultsTreeItemClick)
        
        
        resultsViewLayout.addWidget(self.v1Label)
        resultsViewLayout.addWidget(self.view1)
        resultsViewLayout.addWidget(self.v2Label)
        resultsViewLayout.addWidget(self.view2)
        resultsViewLayout.addWidget(gb)
        self.resultssplitter.addWidget(resultsViewWidget)
        self.resultssplitter.addWidget(self.trResults)

    def dataReady(self):

        cursor = self.consoleOutput.textCursor()
        cursor.movePosition(cursor.End)
        bytes1 = str(self.process.readAll())
#         print bytes1
        cursor.insertText(bytes1)
        self.consoleOutput.ensureCursorVisible()
        
    def runFinished(self):
        self.btnRun.setEnabled(True)
        self.tabs.setCurrentIndex(2)
        self.loadResultFileClicked()
        
    def addConsoleItems(self):

        self.process.started.connect(lambda: self.btnRun.setEnabled(False))
        self.process.finished.connect(self.runFinished)
        self.process.readyReadStandardOutput.connect(self.dataReady)
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels);

        self.consolelayout.addWidget(self.consoleOutput)
        self.consoleOutput.setFontPointSize(9.0)
        font = QtGui.QFont("Monospace")
        font.setPointSize(6)
        self.consoleOutput.setFont(font)

    def addTopLevel(self):
        self.splitter = QtGui.QSplitter(1)
        self.setCentralWidget(self.splitter)
        self.vboxwidget = QtGui.QWidget()
        self.vbox = QtGui.QVBoxLayout()
        self.tabs = QtGui.QTabWidget()
        self.widg1 = QtGui.QWidget()
        self.widg2 = QtGui.QWidget()
        self.widg3 = QtGui.QWidget()
        self.wavelayout = QtGui.QGridLayout()

        self.widg1.setLayout(self.wavelayout)
        self.widg3.setLayout(self.consolelayout)
        self.tabs.addTab(self.widg3, 'Console')
        self.tabs.addTab(self.resultssplitter, 'Results')

        self.vboxwidget.setLayout(self.vbox)
        self.splitter.addWidget(self.vboxwidget)
        self.splitter.addWidget(self.tabs)
        self.splitter.setSizes([300, 724])

    def jsonDictToParamDict(self, d):
        params = []
        # {'name': 'Basic parameter data types', 'type': 'group', 'children': []}
        for key, value in d.iteritems():
            if isinstance(value, dict):
                children = self.jsonDictToParamDict(value)
                params.append({'name': key, 'type':'group', 'expanded':False, 'children': children})
            else:
                params.append({'name': key, 'type':'str', 'value':value })
        return params

    def paramDictToJsonDict(self, d):
        params = {}
        for entry in d['children'].items():
            if entry[1]['value'] is not None:
                params[entry[0]] = entry[1]['value']
            else:
                params[entry[0]] = self.paramDictToJsonDict(entry[1])
            pass
        return params

    def addLeftSide(self):
        self.loadconfiglayout = QtGui.QHBoxLayout()
        self.saveconfiglayout = QtGui.QHBoxLayout()
        self.loadresultlayout = QtGui.QHBoxLayout()
        self.loadconfigwidget = QtGui.QWidget()
        self.saveconfigwidget = QtGui.QWidget()
        self.loadresultwidget = QtGui.QWidget()
        self.loadconfigwidget.setLayout(self.loadconfiglayout)
        self.saveconfigwidget.setLayout(self.saveconfiglayout)
        self.loadresultwidget.setLayout(self.loadresultlayout)

        self.btnRun = QtGui.QPushButton('Save & Run')
        self.loadPathTb = QtGui.QLineEdit()
        self.loadResultsPathTb = QtGui.QLineEdit()
        self.loadResultsPathTb.setText(defaultResultsPath)
        self.loadPathTb.setText(defaultLoadConfigPath)
        self.savePathTb = QtGui.QLineEdit()
        self.savePathTb.setText(defaultSaveConfigPath)


        self.btnLoadConfigFile.clicked.connect(self.loadConfigFileClicked)
        self.btnSaveConfigFile.clicked.connect(self.saveConfigFileClicked)
        self.btnLoadResults.clicked.connect(self.loadResultFileClicked)
        self.btnRun.clicked.connect(self.runBtnCLicked)

        self.loadconfiglayout.addWidget(self.loadPathTb)
        self.loadconfiglayout.addWidget(self.btnLoadConfigFile)
        self.saveconfiglayout.addWidget(self.savePathTb)
        self.saveconfiglayout.addWidget(self.btnSaveConfigFile)
        self.loadresultlayout.addWidget(self.loadResultsPathTb)
        self.loadresultlayout.addWidget(self.btnLoadResults)

        self.paramTree = ParameterTree()
        self.vbox.addWidget(QtGui.QLabel('load config path'))
        self.vbox.addWidget(self.loadconfigwidget)
        self.vbox.addWidget(QtGui.QLabel('load result file'))
        self.vbox.addWidget(self.loadresultwidget)
        self.vbox.addWidget(QtGui.QLabel('save config path'))
        self.vbox.addWidget(self.saveconfigwidget)
        self.vbox.addWidget(self.paramTree)
        self.vbox.addWidget(self.btnRun)
        self.loadConfigFileClicked()

    def runBtnCLicked(self):
        command = defaultExePath
        self.consoleOutput.clear()
        self.tabs.setCurrentIndex(1)
        self.saveConfigFileClicked()
        self.process.start(command, [self.savePathTb.text()])

    def loadConfigFileClicked(self):
        fn = self.loadPathTb.text()
        with open(fn, 'r') as f:
            d = json.load(f)
        p = self.jsonDictToParamDict(d)
        self.par = Parameter.create(name='params', type='group', children=p)
        self.paramTree.setParameters(self.par, showTop=False)

    def saveConfigFileClicked(self):
        d = self.paramDictToJsonDict(self.par.saveState())
        fn = self.savePathTb.text()
        with open(fn, 'w') as f:
            json.dump(d, f)

    def loadResultFileClicked(self):
        fn = self.loadResultsPathTb.text()
        s = QtCore.QString()
        fn = str(fn)
        with pg.BusyCursor():
            self.data = h5read(fn)
        d = self.data.copy()
        d1 = {key: str(d[key].shape) for key, value in d.items() if key != 'config'}

        self.trResults.setData(d1, hideRoot=True)

    def resultsTreeItemClick(self, item, column):
        self.selected = self.data[str(item.text(0))]
        self.updateResultsView()

    def updateResultsView(self):  
        if self.resultsRadio1.isChecked():
                self.v1Label.setText(self.labelText[0][0])
                self.v2Label.setText(self.labelText[0][1])
        else:
                self.v1Label.setText(self.labelText[1][0])
                self.v2Label.setText(self.labelText[1][1])
        if self.selected.shape[0] > 1:
            if self.resultsRadio1.isChecked():
                self.view1.setImage(np.abs(self.selected), autoRange=True, autoLevels=True, levels=None, axes=None, \
                                        xvals=None, pos=None, scale=None, transform=None, autoHistogramRange=True)
                self.view2.setImage(np.angle(self.selected), autoRange=True, autoLevels=True, levels=None, axes=None, \
                                        xvals=None, pos=None, scale=None, transform=None, autoHistogramRange=True) 
            else:
                self.view1.setImage(np.real(self.selected), autoRange=True, autoLevels=True, levels=None, axes=None, \
                                        xvals=None, pos=None, scale=None, transform=None, autoHistogramRange=True)
                self.view2.setImage(np.imag(self.selected), autoRange=True, autoLevels=True, levels=None, axes=None, \
                                        xvals=None, pos=None, scale=None, transform=None, autoHistogramRange=True)  

if __name__ == '__main__':


    app = QtGui.QApplication([])
    gui = SliceGUI()
    gui.resize(1920, 1080)
    gui.show()
    gui.setWindowTitle('slice++')
    app.exec_()


