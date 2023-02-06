from PySide2.QtWidgets import *
from PySide2.QtUiTools import QUiLoader
from PySide2.QtGui import *
import FDMundergroundwater.onedimensionflow as fo
import FDMundergroundwater.twodimensionsflow as ft
import time


class One_dimension_confined_aquifer_stable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/odcasf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Confined_aquifer_SF()

    def flow_draw(self):
        self.flow.transmissivity(self.ui.transmissivity.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.length(self.ui.length.toPlainText())
        self.flow.draw(self.flow.solve())

    def return_main(self):
        self.ui.close()


class One_dimension_confined_aquifer_unstable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/odcausf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 数值解存放
        self.solve_fdm = None
        # 解析解存放
        self.solve_as = None
        # 傅里叶级数存放
        self.fourier_series = None
        # 解析解分配CPU核心数
        self.cpu_cores = None
        # 监测按钮《计算数值解》
        self.ui.solve.clicked.connect(self.solve)
        # 监测按钮《计算解析解》
        self.ui.solve_analytic_solution.clicked.connect(self.solve_analytic_solution)
        # 监测按钮《绘制数值解》
        self.ui.draw_solve.clicked.connect(self.draw_solve)
        # 监测按钮《绘制解析解》
        self.ui.draw_solve_analytic_solution.clicked.connect(self.draw_solve_analytic_solution)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Confined_aquifer_USF()
        # 获取当前系统时间戳
        t = time.localtime()
        self.ui.textBrowser.append('程序运行时间（北京时间）：')
        self.ui.textBrowser.append(str(time.strftime("%Y-%m-%d %H:%M:%S", t)))

    def solve(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.step_time(self.ui.step_time.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        # 判断填写的是压力模量系数还是导水系数和贮水系数
        if self.ui.pressure_diffusion_coefficient.toPlainText() == '':
            self.flow.storativity(self.ui.storativity.toPlainText())
            self.flow.transmissivity(self.ui.transmissivity.toPlainText())
            self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        else:
            self.flow.pressure_diffusion_coefficient(self.ui.pressure_diffusion_coefficient.toPlainText())
            self.flow.transmissivity(1)
            self.flow.leakage_recharge("0")
        self.ui.textBrowser.append('正在进行数值解求解')
        start_time = time.perf_counter()
        self.solve_fdm = self.flow.solve()
        end_time = time.perf_counter()
        self.ui.textBrowser.append('计算完毕，用时' + str(end_time - start_time) + '秒')

    def draw_solve(self):
        self.flow.draw(self.solve_fdm)

    def solve_analytic_solution(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.step_time(self.ui.step_time.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        self.fourier_series = self.ui.spinBox.value()
        self.cpu_cores = self.ui.verticalSlider.value()
        # 判断填写的是压力模量系数还是导水系数和贮水系数
        if self.ui.pressure_diffusion_coefficient.toPlainText() == '':
            self.flow.storativity(self.ui.storativity.toPlainText())
            self.flow.transmissivity(self.ui.transmissivity.toPlainText())
            self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        else:
            self.flow.pressure_diffusion_coefficient(self.ui.pressure_diffusion_coefficient.toPlainText())
            self.flow.transmissivity(1)
            self.flow.leakage_recharge("0")
        self.ui.textBrowser.append('正在进行解析解求解')
        start_time = time.perf_counter()
        self.solve_as = self.flow.solve_multi(fourier_series=self.fourier_series, cpu_cores=self.cpu_cores)
        end_time = time.perf_counter()
        self.ui.textBrowser.append('计算完毕，用时' + str(end_time - start_time) + '秒')

    def draw_solve_analytic_solution(self):
        title = '傅里叶级数解（傅里叶级数取前' + str(self.fourier_series) + '项)，' + '分配CPU核心' + str(self.cpu_cores) + '个'
        self.flow.draw(self.solve_as, title=title)

    def return_main(self):
        self.ui.close()


class One_dimension_unconfined_aquifer_stable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/oduasf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Unconfined_aquifer_SF()

    def flow_draw(self):
        self.flow.hydraulic_conductivity(self.ui.hydraulic_conductivity.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.length(self.ui.length.toPlainText())
        self.flow.reference_thickness(self.ui.reference_thickness.toPlainText())
        self.flow.draw(self.flow.solve())

    def return_main(self):
        self.ui.close()


class One_dimension_unconfined_aquifer_unstable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/oduausf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Unconfined_aquifer_USF()

    def flow_draw(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        # 判断填写的是压力模量系数还是渗透系数，给水度和参考厚度
        if self.ui.pressure_diffusion_coefficient.toPlainText() == '':
            self.flow.storativity(self.ui.storativity.toPlainText())
            self.flow.hydraulic_conductivity(self.ui.hydraulic_conductivity.toPlainText())
            self.flow.reference_thickness(self.ui.reference_thickness.toPlainText())
            self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        else:
            self.flow.pressure_diffusion_coefficient(self.ui.pressure_diffusion_coefficient.toPlainText())
            self.flow.reference_thickness(self.ui.reference_thickness.toPlainText())
            self.flow.hydraulic_conductivity(1)
            self.flow.leakage_recharge("0")
        self.flow.draw(self.flow.solve())

    def return_main(self):
        self.ui.close()


class Two_dimension_confined_aquifer_stable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/tdcasf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = ft.Confined_aquifer_SF()

    def flow_draw(self):
        self.flow.transmissivity(self.ui.transmissivity.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.flow.t_boundary(self.ui.t_boundary.toPlainText())
        self.flow.b_boundary(self.ui.b_boundary.toPlainText())
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.y_length(self.ui.y_length.toPlainText())
        self.flow.draw(self.flow.solve())

    def return_main(self):
        self.ui.close()


class Two_dimension_unconfined_aquifer_stable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/tduasf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = ft.Unconfined_aquifer_SF()

    def flow_draw(self):
        self.flow.hydraulic_conductivity(self.ui.hydraulic_conductivity.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.flow.t_boundary(self.ui.t_boundary.toPlainText())
        self.flow.b_boundary(self.ui.b_boundary.toPlainText())
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.y_length(self.ui.y_length.toPlainText())
        self.flow.reference_thickness(self.ui.reference_thickness.toPlainText())
        self.flow.draw(self.flow.solve())

    def return_main(self):
        self.ui.close()


class Two_dimension_confined_aquifer_unstable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/tdcausf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 监测按钮《上一时刻》
        self.ui.previous_time.clicked.connect(self.previous_time)
        # 监测按钮《下一时刻》
        self.ui.next_time.clicked.connect(self.next_time)
        # 获取自编库中的类的用法
        self.flow = ft.Confined_aquifer_USF()
        # 存储水头解值的列表
        self.h_all_time = []
        self.time_location = 0  # 时刻位置
        self.time_all = 0  # 抛开初始时刻的所有时刻个数

    def flow_draw(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.t_boundary(self.ui.t_boundary.toPlainText())
        self.flow.b_boundary(self.ui.b_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.step_time(self.ui.step_time.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.y_length(self.ui.y_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        self.flow.storativity(self.ui.storativity.toPlainText())
        self.flow.transmissivity(self.ui.transmissivity.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.h_all_time = self.flow.solve()
        self.flow.draw(self.h_all_time[0])  # 绘制初始时刻的水头值
        self.time_all = len(self.h_all_time) - 1
        self.time_location = 0
        self.ui.progressBar.reset()
        self.ui.progressBar.setValue(0)  # 进度条置为0

    def next_time(self):
        self.time_location += 1
        self.flow.draw(self.h_all_time[self.time_location])
        self.ui.progressBar.reset()
        T = int((self.time_location / self.time_all) * 100)
        self.ui.progressBar.setValue(T)  # 设置进度条进度

    def previous_time(self):
        self.time_location -= 1
        self.flow.draw(self.h_all_time[self.time_location])
        self.ui.progressBar.reset()
        T = int((self.time_location / self.time_all) * 100)
        self.ui.progressBar.setValue(T)  # 设置进度条进度

    def return_main(self):
        self.ui.close()


class Two_dimension_unconfined_aquifer_unstable_flow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 从文件中加载ui格式
        self.ui = QUiLoader().load("ui/tduausf.ui")
        self.ui.setWindowIcon(QIcon("water.ico"))
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 监测按钮《上一时刻》
        self.ui.previous_time.clicked.connect(self.previous_time)
        # 监测按钮《下一时刻》
        self.ui.next_time.clicked.connect(self.next_time)
        # 获取自编库中的类的用法
        self.flow = ft.Unconfined_aquifer_USF()
        # 存储水头解值的列表
        self.h_all_time = []
        self.time_location = 0  # 时刻位置
        self.time_all = 0  # 抛开初始时刻的所有时刻个数

    def flow_draw(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.t_boundary(self.ui.t_boundary.toPlainText())
        self.flow.b_boundary(self.ui.b_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
        self.flow.step_time(self.ui.step_time.toPlainText())
        self.flow.x_length(self.ui.x_length.toPlainText())
        self.flow.y_length(self.ui.y_length.toPlainText())
        self.flow.t_length(self.ui.t_length.toPlainText())
        self.flow.initial_condition(self.ui.initial_condition.toPlainText())
        self.flow.storativity(self.ui.storativity.toPlainText())
        self.flow.hydraulic_conductivity(self.ui.hydraulic_conductivity.toPlainText())
        self.flow.reference_thickness(self.ui.reference_thickness.toPlainText())
        self.flow.leakage_recharge(self.ui.leakage_recharge.toPlainText())
        self.h_all_time = self.flow.solve()
        self.flow.draw(self.h_all_time[0])  # 绘制初始时刻的水头值
        self.time_all = len(self.h_all_time) - 1
        self.time_location = 0
        self.ui.progressBar.reset()
        self.ui.progressBar.setValue(0)  # 进度条置为0

    def next_time(self):
        self.time_location += 1
        self.flow.draw(self.h_all_time[self.time_location])
        self.ui.progressBar.reset()
        T = int((self.time_location / self.time_all) * 100)
        self.ui.progressBar.setValue(T)  # 设置进度条进度

    def previous_time(self):
        self.time_location -= 1
        self.flow.draw(self.h_all_time[self.time_location])
        self.ui.progressBar.reset()
        T = int((self.time_location / self.time_all) * 100)
        self.ui.progressBar.setValue(T)  # 设置进度条进度

    def return_main(self):
        self.ui.close()
