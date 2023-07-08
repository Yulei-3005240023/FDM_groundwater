from PySide2.QtWidgets import *
from PySide2.QtUiTools import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import FDMgroundwater.randomflow as ra
import sys
import time


class Set_hydrogeological_parameter(QDialog):
    # 定义信号
    signal_hydrogeological_parameter = Signal(list)

    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/Set_hydrogeological_parameter.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《确认》
        self.ui.sure.clicked.connect(self.sure)
        # 监测按钮《返回》
        self.ui.back.clicked.connect(self.back)

    def sure(self):
        hydraulic_conductivity = self.ui.doubleSpinBox_hydraulic_conductivity.value()
        specific_yield = self.ui.doubleSpinBox_specific_yield.value()
        hydrogeological_parameter = [hydraulic_conductivity, specific_yield]
        self.signal_hydrogeological_parameter.emit(hydrogeological_parameter)
        self.ui.close()

    def back(self):
        self.ui.close()


class Set_FDM(QDialog):
    # 定义信号
    signal_FDM_parameter = Signal(list)

    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/Set_FDM.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《确认》
        self.ui.sure.clicked.connect(self.sure)
        # 监测按钮《返回》
        self.ui.back.clicked.connect(self.back)

    def sure(self):
        step_length = self.ui.spinBox_step_length.value()
        step_time = self.ui.spinBox_step_time.value()
        FDM_parameter = [step_length, step_time]
        self.signal_FDM_parameter.emit(FDM_parameter)
        self.ui.close()

    def back(self):
        self.ui.close()


class Random_one_dimension_boussinesq_window(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/Random_flow.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 实例化具象：潜水含水层一维随机流
        self.flow = ra.Random_one_dimension_boussinesq()
        # 监测按钮《数值解计算》
        self.ui.solve_FDM.clicked.connect(self.solve_FDM)
        # 实例化具象：设置水文地质学参数
        self.Set_hydrogeological_parameter_window = Set_hydrogeological_parameter()
        # 设置水文地质学参数
        self.ui.actionSet_hydrogeological_parameter.triggered.connect(self.actionSet_hydrogeological_parameter)
        # 默认渗透系数为10，给水度为0.1
        self.flow.hydraulic_conductivity(10.0)
        self.flow.specific_yield(0.1)
        # 实例化具象：设置数值解求解参数
        self.Set_FDM_parameter_window = Set_FDM()
        # 默认空间差分步长为1，时间差分步长为1
        self.flow.step_length(1)
        self.flow.step_time(1)
        # 设置数值解求解参数
        self.ui.actionSet_FDM.triggered.connect(self.actionSet_FDM_parameter)
        # 计算得到的水头存储
        self.solve_fdm = None
        # 获取当前系统时间戳
        t = time.localtime()
        # 日志时间
        self.time = str(time.strftime("%Y-%m-%d_%H时%M分%S秒", t))
        self.ui.textBrowser.append('程序运行时间（北京时间）：' + str(time.strftime("%Y-%m-%d %H:%M:%S", t)))
        self.ui.textBrowser.append(self.time)

    def actionSet_hydrogeological_parameter(self):
        self.Set_hydrogeological_parameter_window.ui.show()
        self.Set_hydrogeological_parameter_window.signal_hydrogeological_parameter.connect(
            self.get_hydrogeological_parameter)

    def get_hydrogeological_parameter(self, hydrogeological_parameter):  # 主窗口获得水文地质学参数的槽函数
        self.flow.hydraulic_conductivity(hydrogeological_parameter[0])
        self.flow.specific_yield(hydrogeological_parameter[1])

    def actionSet_FDM_parameter(self):
        self.Set_FDM_parameter_window.ui.show()
        self.Set_FDM_parameter_window.signal_FDM_parameter.connect(self.get_FDM_parameter)

    def get_FDM_parameter(self, FDM_parameter):  # 主窗口获得数值解求解参数的槽函数
        self.flow.step_length(FDM_parameter[0])
        self.flow.step_time(FDM_parameter[1])
        print(self.flow.st)

    def left_boundary(self):
        if self.ui.comboBox_left_boundary.currentText() == '一类边界（给定水头）':
            self.flow.l_boundary(self.ui.doubleSpinBox_left_boundary.value(), Dirichlet=True)
            self.ui.textBrowser.append(
                '左边界（给定水头边界）：\n水头值为：' + str(self.ui.doubleSpinBox_left_boundary.value()))
        else:
            self.flow.l_boundary(self.ui.doubleSpinBox_left_boundary.value(), Neumann=True)
            self.ui.textBrowser.append(
                '左边界（给定通量边界）：\n通量值为：' + str(self.ui.doubleSpinBox_left_boundary.value()))

    def right_boundary(self):
        if self.ui.comboBox_right_boundary.currentText() == '一类边界（给定水头）':
            self.flow.r_boundary(self.ui.doubleSpinBox_right_boundary.value(), Dirichlet=True)
            self.ui.textBrowser.append(
                '右边界（给定水头边界）：\n水头值为：' + str(self.ui.doubleSpinBox_right_boundary.value()))
        else:
            self.flow.r_boundary(self.ui.doubleSpinBox_right_boundary.value(), Neumann=True)
            self.ui.textBrowser.append(
                '右边界（给定通量边界）：\n通量值为：' + str(self.ui.doubleSpinBox_right_boundary.value()))

    def solve_FDM(self):
        self.left_boundary()
        self.right_boundary()
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.ui.textBrowser.append('正在进行数值解求解')
        start_time = time.perf_counter()
        self.solve_fdm = self.flow.solve()
        end_time = time.perf_counter()
        self.ui.textBrowser.append('计算完毕，用时' + str(end_time - start_time) + '秒')

    def set_time_choose_box(self):  # 设置时刻选择条
        time_all = int(self.flow.tl / self.flow.st) + 1
        self.ui.spinBox_time.setMaximum(time_all - 1)
        self.ui.textBrowser_time.setPlainText(
            '计算时长为：' + str(self.flow.tl) + '天，时间分割步长为：' + str(self.flow.st) + '天，共计有时刻' + str(
                time_all) + '个。')


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    # QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication()
    window = Random_one_dimension_boussinesq_window()
    window.ui.show()
    sys.exit(app.exec_())
