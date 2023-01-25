import FDMundergroundwater.onedimensionflow as fo
from PySide2.QtWidgets import QApplication, QMessageBox, QMainWindow
from PySide2.QtUiTools import QUiLoader
from multiprocessing import Process


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/MainWindow - untitled.ui")
        self.ui.pushButton.clicked.connect(self.hello)

    def hello(self):
        flow_dimension = self.ui.choose_flow_dimension.checkedButton()
        aquifer_type = self.ui.choose_aquifer_type.checkedButton()
        stable_status = self.ui.choose_stable_status.checkedButton()
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流承压含水层稳定流":

            print("fuck you")
            p = Process(target=one_dimension_confined_aquifer_flow())
            p.start()

def one_dimension_confined_aquifer_flow():
    class One_dimension_confined_aquifer_flow:
        def __init__(self):
            # 从文件中加载ui格式
            self.ui = QUiLoader().load("ui/MainWindow - untitled.ui")
            # 获取自编库中的类的用法
            self.flow = fo.Confined_aquifer_SF()

    app_ = QApplication([])
    window_ = One_dimension_confined_aquifer_flow()
    window_.ui.show()
    app_.exec_()


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    window.ui.show()
    app.exec_()
