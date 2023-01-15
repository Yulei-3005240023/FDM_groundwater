import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import sympy as sy
from sympy import symbols


class Stableflow:
    def __init__(self):
        self.name_chinese = "稳定一维流"
        self.xl = None
        self.yl = None
        self.sl = None
        self.h_r = None
        self.h_l = None

    def l_boundary(self, h_l, Dirichlet=True, Neumann=False, Robin=False):  # 左边界
        if Dirichlet:
            self.h_l = float(h_l)

    def r_boundary(self, h_r, Dirichlet=True, Neumann=False, Robin=False):  # 右边界
        if Dirichlet:
            self.h_r = float(h_r)

    def step_length(self, sl):  # 差分步长,此处默认X,Y轴
        self.sl = float(sl)

    def x_length(self, xl):  # X轴轴长
        self.xl = float(xl)

    def y_length(self, yl):  # Y轴轴长
        self.yl = float(yl)
