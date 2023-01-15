import numpy as np
import sympy as sy
from sympy import symbols
import numpy.linalg as nla
import matplotlib.pyplot as plt


class Stableflow:

    def __init__(self):
        self.name_chinese = "稳定一维流"
        self.l = None
        self.sl = None
        self.h_r = None
        self.h_l = None

    def l_boundary(self, h_l, Dirichlet=True, Neumann=False, Robin=False):  # 左边界
        if Dirichlet:
            self.h_l = float(h_l)

    def r_boundary(self, h_r, Dirichlet=True, Neumann=False, Robin=False):  # 右边界
        if Dirichlet:
            self.h_r = float(h_r)

    def step_length(self, sl):  # 差分步长
        self.sl = float(sl)

    def length(self, l):  # 轴长
        self.l = float(l)

    def draw(self, H_ALL):
        # X轴单元格的数目
        m = int(self.l // self.sl) + 1
        # X轴
        X = np.linspace(0, self.l, m)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot()
        ax.plot(X, H_ALL, linewidth=1, antialiased=True)

        def H_y(h_all):
            hy = 0
            for i in h_all:
                if i > hy:
                    hy = i

        ax.set_ylim(0, H_y(H_ALL))
        plt.title("差分数值解(差分步长%s)" % self.sl)
        plt.show()


class Confined_aquifer_SF(Stableflow):

    def __init__(self):
        super().__init__()
        self.name_chinese = "承压含水层稳定一维流"
        self.w = None
        self.t = None

    def transmissivity(self, t=1):  # 承压含水层导水系数的设定，可以为一个常数也可以为带有前缀为sy.的函数，如sy.sin(x)
        self.t = str(t)

    def leakage_recharge(self, w=0):  # 承压含水层越流补给源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
        self.w = str(w)

    def solve(self):
        # 对于承压含水层一维稳定流，定义参数 x
        x = symbols("x")

        # 对导水系数T(x)
        def T(x):
            return eval(self.t)

        # 对源汇项函数W(x)
        def W(x):
            return eval(self.w)

        # X轴差分点的数目
        m = int(self.l // self.sl) + 1
        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros(m)
        # 常数b矩阵w
        H_b = np.zeros((m, 1))
        # 系数a矩阵
        H_a = np.zeros((m, m))
        for i in range(0, m):  # 对列进行扫描
            # 源汇项赋值
            H_b[i] = H_b[i] - W(i * self.sl)
            # 左右边界赋值
            if (i - 1) < 0:
                H_b[i] = H_b[i] - (T(i * self.sl) - self.sl * (sy.diff(T(x), x).subs(x, i * self.sl))) * self.h_l
            if (i + 1) == m:
                H_b[i] = H_b[i] - T(self.l) * self.h_r
            # 给位置为(i-1)处的水头赋上系数值
            if (i - 1) >= 0:
                H_a[i, (i - 1)] = T(i * self.sl) - (self.sl * sy.diff(T(x), x).subs(x, i * self.sl))
            # 给位置为(i+1)处的水头赋上系数值
            if (i + 1) < m:
                H_a[i, (i + 1)] = T(self.l)
            # 给位置为(i)处的水头赋上系数值
            H_a[i, i] = self.sl * sy.diff(T(x), x).subs(x, i * self.sl) - 2 * T(self.l)
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for i in range(0, m):
            H_ALL[i] = H[i]
        print(H_a)
        print(H_b)
        print(H)
        return H_ALL


class Unstableflow:
    def __init__(self):
        self.tl = None
        self.st = None
        self.name_chinese = "非稳定一维流"
        self.xl = None
        self.sl = None
        self.h_r = None
        self.h_l = None

    def l_boundary(self, h_l, Dirichlet=True, Neumann=False, Robin=False):  # 左边界
        if Dirichlet:
            self.h_l = float(h_l)

    def r_boundary(self, h_r, Dirichlet=True, Neumann=False, Robin=False):  # 右边界
        if Dirichlet:
            self.h_r = float(h_r)

    def step_length(self, sl):  # X轴差分步长
        self.sl = float(sl)

    def step_time(self, st):  # 时间轴差分步长
        self.st = float(st)

    def x_length(self, xl):  # X轴轴长
        self.xl = float(xl)

    def t_length(self, tl):  # 时间轴轴长，原则上单位为天
        self.tl = float(tl)

    def draw(self, H_ALL):
        # X轴单元格的数目
        m = int(self.xl // self.sl) + 1
        # 时间轴单元格数目
        n = int(self.tl // self.st) + 1
        # X轴
        X = np.linspace(0, self.xl, m)
        # 时间轴
        T = np.linspace(0, self.tl, n)
        # 定义初值
        X, T = np.meshgrid(X, T)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot('3d')
        ax.plot_surface(X, T, H_ALL, linewidth=0, antialiased=True, cmap=plt.get_cmap('rainbow'))

        def H_z(h_all):
            hz = 0
            for i in h_all:
                if i > hz:
                    hz = i

        ax.set_zlim(0, H_z(H_ALL))
        plt.title("差分数值解(差分步长%s)" % self.sl)
        plt.show()


class Confined_aquifer_USF(Unstableflow):
    def __init__(self):
        super().__init__()
        self.name_chinese = "承压含水层一维非稳定流"
        self.w = None
        self.a = None

    def transmissivity(self, a=1):  # 承压含水层压力扩散系数的设定，等于
        self.a = str(a)

    def leakage_recharge(self, w=0):  # 承压含水层越流补给源汇项的设定。
        self.w = str(w)

    def solve(self):
        # 对于承压含水层一维非稳定流，定义两个参数 x t
        x = symbols("x")
        t = symbols("t")
        # X轴差分点的数目
        m = int(self.xl // self.sl) + 1
        # 时间轴差分点的数目
        n = int(self.tl // self.st) + 1

        # 对源汇项函数W(x)
        def W(x):
            return eval(self.w)

        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n, m * n))
        for k in range(0, n):  # 对行(时间轴)进行扫描
            for i in range(0, m):  # 对列(X轴)进行扫描
                # 源汇项赋值
                H_b[i] = H_b[i] - W(i * self.sl)
                # 左右边界赋值
                if (i - 1) < 0:
                    H_b[i] = H_b[i] - (T(i * self.sl) - self.sl * (sy.diff(T(x), x).subs(x, i * self.sl))) * self.h_l
                if (i + 1) == m:
                    H_b[i] = H_b[i] - T(self.l) * self.h_r
                # 给位置为(i-1)处的水头赋上系数值
                if (i - 1) >= 0:
                    H_a[i, (i - 1)] = T(i * self.sl) - (self.sl * sy.diff(T(x), x).subs(x, i * self.sl))
                # 给位置为(i+1)处的水头赋上系数值
                if (i + 1) < m:
                    H_a[i, (i + 1)] = T(self.l)
                # 给位置为(i)处的水头赋上系数值
                H_a[i, i] = self.sl * sy.diff(T(x), x).subs(x, i * self.sl) - 2 * T(self.l)
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for k in range(0, n):  # 对时间进行扫描
            for i in range(0, m):  # 对空间进行扫描
                H_ALL[k, i] = H[k * n + i]
        print(H_a)
        print(H_b)
        print(H)
        return H_ALL


a = Confined_aquifer_USF()
a.x_length(10)
a.t_length(10)
a.step_length(0.05)
a.step_time(1)
a.l_boundary(0)
a.r_boundary(3)
a.transmissivity(1)
a.leakage_recharge(3)
b = a.solve()
print(b)
a.draw(b)
