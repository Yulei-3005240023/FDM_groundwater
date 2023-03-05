import gevent
from PySide2.QtWidgets import QApplication, QMainWindow, QTreeWidget
from PySide2.QtUiTools import QUiLoader
import childwindow as cw
from PySide2.QtGui import QIcon
import sys
import multiprocessing


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/MainWindow - untitled.ui")
        # 加载图标
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 为按钮添加点击动作
        self.ui.toth.clicked.connect(self.toth)
        self.ui.treeWidget.clicked.connect(self.next)
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
        # 实例化具象：承压含水层二维非稳定流
        self.two_dimension_confined_aquifer_unstable_flow_window = cw.Two_dimension_confined_aquifer_unstable_flow()
        # 实例化具象：潜水含水层二维非稳定流
        self.two_dimension_unconfined_aquifer_unstable_flow_window = cw.Two_dimension_unconfined_aquifer_unstable_flow()
        # 实例化具象：Tóth复杂盆地
        self.toth_difficult_basin_window = cw.Two_dimension_Toth_difficult_baisn()

    def next(self):  # 该函数用于打开每一种水流模式所对应的主窗口
        item = self.ui.treeWidget.currentItem()
        # 水流模式的条件判断
        if item.whatsThis(0) == "一维流承压含水层稳定流":
            task1 = gevent.spawn(self.one_dimension_confined_aquifer_stable_flow_window.ui.show())  # 创建多进程任务
            task_list.append(task1)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "一维流承压含水层非稳定流":
            task2 = gevent.spawn(self.one_dimension_confined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task2)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "一维流潜水含水层稳定流":
            task3 = gevent.spawn(self.one_dimension_unconfined_aquifer_stable_flow_window.ui.show())
            task_list.append(task3)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "一维流潜水含水层非稳定流":
            task4 = gevent.spawn(self.one_dimension_unconfined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task4)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "二维流承压含水层稳定流":
            task5 = gevent.spawn(self.two_dimension_confined_aquifer_stable_flow_window.ui.show())
            task_list.append(task5)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "二维流潜水含水层稳定流":
            task6 = gevent.spawn(self.two_dimension_unconfined_aquifer_stable_flow_window.ui.show())
            task_list.append(task6)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "二维流承压含水层非稳定流":
            task7 = gevent.spawn(self.two_dimension_confined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task7)  # 把该任务加入到协程列表
        if item.whatsThis(0) == "二维流潜水含水层非稳定流":
            task8 = gevent.spawn(self.two_dimension_unconfined_aquifer_unstable_flow_window.ui.show())
            task_list.append(task8)  # 把该任务加入到协程列表

    def toth(self):  # 进入Tóth复杂盆地的部分
        task9 = gevent.spawn(self.toth_difficult_basin_window.ui.show())
        task_list.append(task9)  # 把该任务加入到协程列表


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    #QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    multiprocessing.freeze_support()  # 预防多进程打包后运行出错
    app = QApplication()
    window = MainWindow()
    task_list = []  # 创建协程列表
    task0 = gevent.spawn(window.ui.show())
    task_list.append(task0)
    gevent.joinall(task_list)
    sys.exit(app.exec_())
