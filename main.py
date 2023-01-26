from PySide2.QtWidgets import QApplication, QMainWindow
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import *
import childwindow as cw
from PySide2.QtGui import QIcon
from multiprocessing import Process
import sys


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/MainWindow - untitled.ui")
        # 加载图标
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 为按钮添加点击动作
        self.ui.pushButton.clicked.connect(self.next)
        # 实例化具象：承压含水层一维稳定流
        self.one_dimension_confined_aquifer_flow_window = cw.One_dimension_confined_aquifer_flow()

    def next(self):  # 该函数用于打开每一种水流模式所对应的主窗口
        # 获取radioButton的数据
        flow_dimension = self.ui.choose_flow_dimension.checkedButton()
        aquifer_type = self.ui.choose_aquifer_type.checkedButton()
        stable_status = self.ui.choose_stable_status.checkedButton()
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流承压含水层稳定流":
            p1 = Process(target=self.one_dimension_confined_aquifer_flow_window.ui.exec())
            p1.start()


def main():
    # QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    window = MainWindow()
    window.ui.show()
    sys.exit(app.exec_())


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    p0 = Process(target=main())
    p0.start()
