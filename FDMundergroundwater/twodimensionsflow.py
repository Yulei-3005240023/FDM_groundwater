import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import sympy as sy
from sympy import symbols


class Stableflow:
    def __init__(self):
        self.h_b = None
        self.h_t = None
        self.name_chinese = "稳定二维流"
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

    def draw(self, H_ALL: np.ndarray):
        # X轴单元格的数目
        m = int(self.xl / self.sl) + 1
        # Y轴单元格的数目
        n = int(self.yl / self.sl) + 1
        # X轴
        X = np.linspace(0, self.xl, m)
        # Y轴
        Y = np.linspace(0, self.yl, n)
        # 定义初值
        X, Y = np.meshgrid(X, Y)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(projection='3d')

        def maxH_z(h_all):
            hz = 0
            for i in h_all:
                for j in i:
                    if j > hz:
                        hz = j
            return hz

        def minH_z(h_all):
            hz = 0
            for i in h_all:
                for j in i:
                    if j < hz:
                        hz = j
            return hz

        ax.set_zlim(minH_z(H_ALL), maxH_z(H_ALL))
        ax.plot_surface(X, Y, H_ALL, linewidth=0, antialiased=True, cmap=plt.get_cmap('rainbow'))
        plt.title("差分数值解(差分步长%s)" % self.sl)
        plt.show()


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

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定，可以设定为x,y的函数。
        self.w = w

    def solve(self):
        # 如果未设定压力扩散系数
        if self.a is None:
            self.a = self.T / self.S
        # 对于承压含水层二维稳定流，定义两个参数 x y
        x = symbols("x")
        y = symbols("y")
        # X,Y轴单元格的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1

        # 对函数W(x)为源汇项函数除以倒水系数
        def W(x, y):
            return eval(self.w) / self.T

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
                    # 源汇项赋值
                    H_b[l_a] = H_b[l_a] - W(j * self.sl, i * self.sl)
                    l_a += 1
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for i in range(0, n):  # 对行进行扫描
            for j in range(0, m):  # 对列进行扫描
                H_ALL[i, j] = H[i * n + j]

        return H_ALL


a = Confined_aquifer_SF()
a.x_length(10)
a.y_length(10)
a.step_length(1)
a.transmissivity(1)
a.storativity(1)
a.l_boundary(0)
a.r_boundary(0)
a.t_boundary(1)
a.b_boundary(0)
a.leakage_recharge("0")
b = a.solve()
a.draw(b)
