import FDMundergroundwater.onedimensionflow as fo
from PySide2.QtWidgets import QApplication, QMessageBox
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import QFile


class MainWindow:
    def __init__(self):
        # 从文件中加载ui格式
        qfile_mainwindow = QFile("ui/MainWindow - untitled.ui")
        qfile_mainwindow.open(QFile.ReadOnly)
        qfile_mainwindow.close()
        self.ui = QUiLoader().load(qfile_mainwindow)


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    cc = fo.Stableflow()
