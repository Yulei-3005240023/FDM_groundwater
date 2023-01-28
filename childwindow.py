from PySide2.QtWidgets import *
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import *
from PySide2.QtGui import *
import FDMundergroundwater.onedimensionflow as fo
import FDMundergroundwater.twodimensionsflow as ft


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
        # 监测按钮《计算并绘图》
        self.ui.draw.clicked.connect(self.flow_draw)
        # 监测按钮《返回上一级》
        self.ui.back.clicked.connect(self.return_main)
        # 获取自编库中的类的用法
        self.flow = fo.Confined_aquifer_USF()

    def flow_draw(self):
        self.flow.l_boundary(self.ui.l_boundary.toPlainText())
        self.flow.r_boundary(self.ui.r_boundary.toPlainText())
        self.flow.step_length(self.ui.step_length.toPlainText())
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
        self.flow.draw(self.flow.solve())

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
