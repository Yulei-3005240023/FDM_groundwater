import gevent
from gevent import monkey

# 给程序打补丁使其变成异步模式
monkey.patch_all()
from PySide2.QtWidgets import QApplication, QMainWindow
from PySide2.QtUiTools import QUiLoader
import childwindow as cw
from PySide2.QtGui import QIcon
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
        self.one_dimension_confined_aquifer_stable_flow_window = cw.One_dimension_confined_aquifer_stable_flow()
        # 实例化具象：承压含水层一维非稳定流
        self.one_dimension_confined_aquifer_unstable_flow_window = cw.One_dimension_confined_aquifer_unstable_flow()
        # 实例化具象：潜水含水层一维稳定流
        self.one_dimension_unconfined_aquifer_stable_flow_window = cw.One_dimension_unconfined_aquifer_stable_flow()
        # 实例化具象：潜水含水层一维非稳定流
        self.one_dimension_unconfined_aquifer_unstable_flow_window = cw.One_dimension_unconfined_aquifer_unstable_flow()
        # 实例化具象：承压含水层二维稳定流
        self.two_dimension_confined_aquifer_stable_flow_window = cw.Two_dimension_confined_aquifer_stable_flow()
        # 实例化具象：潜水含水层二维稳定流
        self.two_dimension_unconfined_aquifer_stable_flow_window = cw.Two_dimension_unconfined_aquifer_stable_flow()

    def next(self):  # 该函数用于打开每一种水流模式所对应的主窗口
        # 获取radioButton的数据
        flow_dimension = self.ui.choose_flow_dimension.checkedButton()
        aquifer_type = self.ui.choose_aquifer_type.checkedButton()
        stable_status = self.ui.choose_stable_status.checkedButton()
        # 水流模式的条件判断
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流承压含水层稳定流":
            task1 = gevent.spawn(self.one_dimension_confined_aquifer_stable_flow_window.ui.show())  # 创建多进程任务
            task_list.append(task1)  # 把该进程加入到进程列表
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流承压含水层非稳定流":
            task2 = gevent.spawn(self.one_dimension_confined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task2)  # 把该进程加入到进程列表
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流潜水含水层稳定流":
            task3 = gevent.spawn(self.one_dimension_unconfined_aquifer_stable_flow_window.ui.show())
            task_list.append(task3)  # 把该进程加入到进程列表
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "一维流潜水含水层非稳定流":
            task4 = gevent.spawn(self.one_dimension_unconfined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task4)  # 把该进程加入到进程列表
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "二维流承压含水层稳定流":
            task5 = gevent.spawn(self.two_dimension_confined_aquifer_stable_flow_window.ui.show())
            task_list.append(task5)  # 把该进程加入到进程列表
        if flow_dimension.text() + aquifer_type.text() + stable_status.text() == "二维流潜水含水层稳定流":
            task6 = gevent.spawn(self.two_dimension_unconfined_aquifer_stable_flow_window.ui.show())
            task_list.append(task6)  # 把该进程加入到进程列表


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    # QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication()
    window = MainWindow()

    task_list = []  # 创建多进程列表
    task0 = gevent.spawn(window.ui.show())
    task_list.append(task0)
    gevent.joinall(task_list)
    sys.exit(app.exec_())
