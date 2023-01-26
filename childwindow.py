from PySide2.QtWidgets import *
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import *
from PySide2.QtGui import *
import FDMundergroundwater.onedimensionflow as fo


class One_dimension_confined_aquifer_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/odcasf.ui")
        # self.ui.setWindowIcon(QIcon("water.ico"))
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Confined_aquifer_SF()
        print("fuck you bitch")

    def return_main(self):
        self.ui.close()
