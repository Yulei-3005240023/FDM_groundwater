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
        #self.ui.toth.clicked.connect(self.toth)
        # 树状水流模式选择
        self.ui.treeWidget.clicked.connect(self.next)
        # 关于程序与算法的列表选择
        self.ui.listWidget_0.clicked.connect(self.algorithm_information)
        # 关于使用的地下水概念的列表选择
        self.ui.listWidget_1.clicked.connect(self.undergroundwater_information)
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

    def algorithm_information(self):
        item = self.ui.listWidget_0.currentItem()
        # 清空文本框
        self.ui.textBrowser_0.clear()
        if item.whatsThis() == '有限差分法':
            self.ui.textBrowser_0.setText('有限差分法(Finite Difference Method)')
            self.ui.textBrowser_0.append('  是求解偏微分方程边值问题和初值问题的一种数值方法，其实质是利用导数的差分近似形式代替偏微分方程形成的差分方程组，通过求解方程组得到离散点的待求变量作为连续场的一种近似结果。')
            self.ui.textBrowser_0.append('  有限差分法中的“有限”，是指网格中的节点、单元或块体数目是有限的。(《地下水运动方程》王旭升等 p99)')
            self.ui.textBrowser_0.append('  在本程序中，使用该方法对地下水含水层的水头进行数值解计算并且绘图。')
        if item.whatsThis() == '差分格式':
            self.ui.textBrowser_0.setText('  在本程序中主要使用向后隐式差分，因为其简单且绝对收敛，该格式的截断误差为o=Δt+Δx^2')
            self.ui.textBrowser_0.append('  在承压含水层一维非稳定流中，还提供了Crank-Nicolson中心差分的方法，该格式的截断误差为o=Δt^2+Δx^2')
        if item.whatsThis() == '差分方程组求解':
            self.ui.textBrowser_0.setText('  对于数值法求解偏微分方程而言，如何使用程序求解方程组是一个重要的问题，在线性代数计算中有几种方法可以求解大型线性方程组：\n')
            self.ui.textBrowser_0.append('1.直接法，其中包括：高斯消元法，LU分解法，追赶法')
            self.ui.textBrowser_0.append('2.迭代法：常见的迭代法有：Jacobi迭代法，Gauss-Seidel迭代法，超松弛迭代法，预条件迭代法')
            self.ui.textBrowser_0.append('  在本程序中使用的是LU分解法，由python扩展库numpy.linalg中的solve函数提供，使用上个世纪90年代编写的lapack程序包的_gesv例程。\n')
            self.ui.textBrowser_0.append('注意：本程序在进行大规模矩阵求解运算时会全额占用计算机的CPU和内存，程序页面暂时卡死是正常现象。')
        if item.whatsThis() == '矩阵赋值':
            self.ui.textBrowser_0.setText('  numpy.linalg.solve(a, b)函数的适用方法为aX=b:')
            self.ui.textBrowser_0.append('  其中a为系数矩阵,b为常数矩阵，对solve函数输入上述两个矩阵会返回求解的X矩阵。因此按照对应的差分格式对矩阵a,b进行赋值后即可带入函数进行求解。')

    def undergroundwater_information(self):
        item = self.ui.listWidget_1.currentItem()
        # 清空文本框
        self.ui.textBrowser_1.clear()
        if item.whatsThis() == '渗透系数':
            self.ui.textBrowser_1.setText('渗透系数(hydraulic conductivity)')
            self.ui.textBrowser_1.append('  渗透系数是一个及其重要的水文地质参数，是表征多空介质透水能力的参数，常用单位为m/d。')
            self.ui.textBrowser_1.append('  渗透系数既与多孔介质的空隙性质有关，也与渗透液体的物理性质有关。(《地下水动力学》陈崇希等 p10,《地下水科学概论》周训等 p44)')
        if item.whatsThis() == '导水系数':
            self.ui.textBrowser_1.setText('导水系数(transmissivity)')
            self.ui.textBrowser_1.append('  虽然渗透系数(K)可以说明岩层的透水能力，但不能单独说明含水层的出水能力。对于承压含水层，由于其厚度(M)是定值，则T=KM也是定制。T称为导水系数，它指的是在水力梯度等于1水流经整个含水层厚度上的单宽流量，常用单位是m2/d。导水系数是表征承压含水层导水能力的参数，只使用于二维流，对于三维流则没有意义。')

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
