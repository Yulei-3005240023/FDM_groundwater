from PySide2.QtWidgets import *
from PySide2.QtUiTools import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import numpy as np
from numpy import sin
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
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


class Set_new_wave(QDialog):
    # 定义信号
    signal_wave_parameter = Signal(list)

    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/Set_new_wave.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 平均每天的降雨量
        self.we = None
        # 监测按钮《保存》
        self.ui.save.clicked.connect(self.save)
        # 监测按钮《返回》
        self.ui.back.clicked.connect(self.back)

    def save(self):
        amplitude = self.ui.doubleSpinBox_amplitude.value()
        cycle = self.ui.doubleSpinBox_cycle.value()
        wave_parameter = [amplitude, cycle]
        self.signal_wave_parameter.emit(wave_parameter)

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
        # 绘图选定时刻存放
        self.time_location = None
        # 监测降雨量期望变化
        self.ui.doubleSpinBox_rain.valueChanged.connect(self.Set_rain_expectation)
        # 监测按钮《绘制数值解线性图》
        self.ui.draw_solve_line.clicked.connect(self.draw_solve_line)

        # 监测按钮《新建一个降雨量波动》
        self.ui.new_wave.clicked.connect(self.Set_new_wave)
        # 监测按钮《随机新建一个降雨量波动》
        self.ui.random_new_wave.clicked.connect(self.random_new_wave)
        # 实例化具象：新建降雨量波动
        self.Set_new_wave_window = Set_new_wave()
        # 降雨量波动函数
        self.rain_function = None
        # 降雨量波动存放
        self.wave_list = []
        # 监测按钮《删除上一个降雨量波动》
        self.ui.delete_wave.clicked.connect(self.delete_wave)
        # 监测按钮《查看时域图像》
        self.ui.time_space_figure.clicked.connect(self.draw_rain_function)

        # 创建Matplotlib绘图的Figure对象和Canvas对象(水头)
        # self.figure_head = plt.figure()
        # self.canvas_head = FigureCanvas(self.figure_head)

        # 创建Matplotlib工具栏并添加到主窗口
        # self.ui.toolbar = NavigationToolbar(self.canvas_head, self)
        # self.ui.addToolBar(self.ui.toolbar)

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

    def Set_rain_expectation(self):
        self.flow.source_sink_expectation(self.ui.doubleSpinBox_rain.value() / (1000 * 365))  # 一年的降雨量期望转化为一天的
        self.flow.source_sink_term(str(self.ui.doubleSpinBox_rain.value()))
        self.ui.textBrowser_rain_function.clear()  # 重置展示降雨量函数的框
        self.wave_list = []  # 重置降雨量波动列表
        self.rain_function = str(self.ui.doubleSpinBox_rain.value() / (1000 * 365))
        self.ui.textBrowser_rain_function.append(self.rain_function)

    def Set_new_wave(self):
        self.Set_new_wave_window.we = (self.ui.doubleSpinBox_rain.value() / 365)  # 一年的降雨量期望转化为一天的
        # 设置文本框支持
        self.Set_new_wave_window.ui.textBrowser.clear()
        self.Set_new_wave_window.ui.textBrowser.append(
            '振幅不得超过一天的降雨量平均值：' + str(
                self.Set_new_wave_window.we) + '\n计算时长为周期的整数倍，为解析解苛求条件。\n根据香农采样定理采样频率必须大于信号频率的两倍。\n所以降雨波动信号的周期确定生成必须大于采样周期的两倍，建议取三倍')
        self.Set_new_wave_window.ui.show()
        self.Set_new_wave_window.signal_wave_parameter.connect(self.get_wave_parameter)

    def get_wave_parameter(self, wave_parameter):  # 主窗口获得波动设置参数的槽函数
        self.wave_list.append('+' + str(wave_parameter[0]/1000) + '*sin(' + str(2 * 3.1415 / wave_parameter[1]) + '*t)')
        self.rain_function = str(self.ui.doubleSpinBox_rain.value() / (1000 * 365))  # 重置降雨量函数
        for i in self.wave_list:
            self.rain_function += i
        self.ui.textBrowser_rain_function.clear()
        self.ui.textBrowser_rain_function.append(self.rain_function)

    def delete_wave(self):
        del self.wave_list[-1]
        self.rain_function = str(self.ui.doubleSpinBox_rain.value() / (1000 * 365))  # 重置降雨量函数
        for i in self.wave_list:
            self.rain_function += i
        self.ui.textBrowser_rain_function.clear()
        self.ui.textBrowser_rain_function.append(self.rain_function)

    def random_new_wave(self):
        self.flow.t_length(self.ui.spinBox_t_length.value())
        amplitude, cycle, frequency = self.flow.random_w()
        self.wave_list.append('+' + str(amplitude) + '*sin(' + str(2 * 3.1415 / frequency) + '*t)')
        self.rain_function = str(self.ui.doubleSpinBox_rain.value() / (1000 * 365))
        for i in self.wave_list:
            self.rain_function += i
        self.ui.textBrowser_rain_function.clear()
        self.ui.textBrowser_rain_function.append(self.rain_function)

    def left_boundary(self):
        if self.ui.comboBox_left_boundary.currentText() == '一类边界（给定水头）':
            self.flow.l_boundary(self.ui.doubleSpinBox_left_boundary.value(), Dirichlet=True)
            self.ui.textBrowser.append(
                '左边界（给定水头边界）：\n水头值为：' + str(self.ui.doubleSpinBox_left_boundary.value()))
        elif self.ui.comboBox_left_boundary.currentText() == '二类边界（给定通量）':
            self.flow.l_boundary(self.ui.doubleSpinBox_left_boundary.value(), Neumann=True)
            self.ui.textBrowser.append(
                '左边界（给定通量边界）：\n通量值为：' + str(self.ui.doubleSpinBox_left_boundary.value()))

    def right_boundary(self):
        if self.ui.comboBox_right_boundary.currentText() == '一类边界（给定水头）':
            self.flow.r_boundary(self.ui.doubleSpinBox_right_boundary.value(), Dirichlet=True)
            self.ui.textBrowser.append(
                '右边界（给定水头边界）：\n水头值为：' + str(self.ui.doubleSpinBox_right_boundary.value()))
        elif self.ui.comboBox_right_boundary.currentText() == '二类边界（给定通量）':
            self.flow.r_boundary(self.ui.doubleSpinBox_right_boundary.value(), Neumann=True)
            self.ui.textBrowser.append(
                '右边界（给定通量边界）：\n通量值为：' + str(self.ui.doubleSpinBox_right_boundary.value()))

    def solve_FDM(self):
        self.left_boundary()
        self.right_boundary()
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        self.flow.x_length(self.ui.spinBox_x_length.value())
        self.flow.t_length(self.ui.spinBox_t_length.value())
        self.flow.source_sink_term(self.ui.textBrowser_rain_function.toPlainText())
        self.ui.textBrowser.append('正在进行数值解求解')
        start_time = time.perf_counter()
        self.solve_fdm = self.flow.solve()
        end_time = time.perf_counter()
        self.ui.textBrowser.append('计算完毕，用时' + str(end_time - start_time) + '秒')
        # 设置时刻选择条
        self.set_time_choose_box()

    def set_time_choose_box(self):  # 设置时刻选择条
        time_all = int(self.flow.tl / self.flow.st) + 1
        self.ui.spinBox_time.setMaximum(time_all - 1)
        self.ui.textBrowser_time.setPlainText(
            '计算时长为：' + str(self.flow.tl) + '天，时间分割步长为：' + str(self.flow.st) + '天，共计有时刻' + str(
                time_all) + '个。')

    def draw_rain_function(self):
        self.flow.t_length(self.ui.spinBox_t_length.value())
        # T轴单元格的数目
        m = int(self.flow.tl / self.flow.sl) + 1
        # T轴
        t = np.linspace(0, self.flow.tl, m)
        # 降雨量
        w = eval(self.rain_function)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_subplot()
        ax.plot(t, w, linewidth=1, antialiased=True)
        ax.set(ylabel='降雨量（m）', xlabel='时间轴（m）')
        plt.title('降雨量时域图像')
        plt.show()

    def draw_head_line(self):
        self.time_location = self.ui.spinBox_time.value()
        title = '数值解，空间差分步长为' + str(self.flow.sl) + '时间差分步长为' + str(
            self.flow.st) + '，绘图时刻为第' + str(self.time_location) + '时刻'
        self.flow.draw(self.solve_fdm, time=self.time_location, title=title)

    def draw_solve_line(self):
        # 判定能否绘图
        if self.solve_fdm is not None:
            pass
        else:
            QMessageBox.critical(self.ui, '错误', '请先进行数值解计算！')  # 未通过校验即报错
            return None  # 结束代码
        self.draw_head_line()
        self.ui.textBrowser.append('第' + str(self.time_location) + '时刻数值解绘图，当前时刻各点解值为：')
        self.ui.textBrowser.append(str(self.solve_fdm[self.time_location]))


# 按间距中的绿色按钮以运行脚本。
if __name__ == '__main__':
    # QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app = QApplication()
    window = Random_one_dimension_boussinesq_window()
    window.ui.show()
    sys.exit(app.exec_())
