import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import sympy as sy
from sympy import symbols


class Stableflow:
    def __init__(self):
        self.h_b = None
        self.h_t = None
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

    def t_boundary(self, h_t, Dirichlet=True, Neumann=False, Robin=False):  # 上边界
        if Dirichlet:
            self.h_t = float(h_t)

    def b_boundary(self, h_b, Dirichlet=True, Neumann=False, Robin=False):  # 下边界
        if Dirichlet:
            self.h_b = float(h_b)

    def step_length(self, sl):  # 差分步长,此处默认X,Y轴差分步长相同
        self.sl = float(sl)

    def x_length(self, xl):  # X轴轴长
        self.xl = float(xl)

    def y_length(self, yl):  # Y轴轴长
        self.yl = float(yl)


class Confined_aquifer_SF(Stableflow):
    def __init__(self):
        super().__init__()
        self.w = None
        self.a = None
        self.S = None
        self.T = None

    def transmissivity(self, T: float):  # 承压含水层导水系数的设定
        self.T = T

    def storativity(self, S: float):  # 承压含水层储水系数（弹性给水度）的设定
        self.S = S

    def pressure_diffusion_coefficient(self, a: float):  # 承压含水层压力扩散系数的设定，等于导水系数除以贮水系数T/S
        self.a = a

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定。
        self.w = w
    def solve(self):
        # X,Y轴单元格的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1
        # 创建一个全部值为0的矩阵
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n, m * n))
        # 系数a矩阵行数
        l_a = 0
        while l_a < m * n:
            for i in range(0, n):  # 对行进行扫描
                for j in range(0, m):  # 对列进行扫描
                    # 上下左右边界赋值
                    if (i - 1) < 0:
                        H_b[l_a] = H_b[l_a] - self.h_t
                    if (j - 1) < 0:
                        H_b[l_a] = H_b[l_a] - self.h_l
                    if (i + 1) == n:
                        H_b[l_a] = H_b[l_a] - self.h_b
                    if (j + 1) == m:
                        H_b[l_a] = H_b[l_a] - self.h_r
                    # 给位置为(i-1,j)处的水头赋上系数值
                    if (i - 1) >= 0:
                        H_a[l_a, (i - 1) * m + j] = 1
                    # 给位置为(i+1,j)处的水头赋上系数值
                    if (i + 1) < n:
                        H_a[l_a, (i + 1) * m + j] = 1
                    # 给位置为(i,j-1)处的水头赋上系数值
                    if (j - 1) >= 0:
                        H_a[l_a, i * m + j - 1] = 1
                    # 给位置为(i,j+1)处的水头赋上系数值
                    if (j + 1) < m:
                        H_a[l_a, i * m + j + 1] = 1
                    # 给位置为（i,j)处的水头赋上系数值
                    H_a[l_a, i * m + j] = -4
                    l_a += 1
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for i in range(0, n):  # 对行进行扫描
            for j in range(0, m):  # 对列进行扫描
                H_ALL[i, j] = H[i * n + j]

        return H_ALL
