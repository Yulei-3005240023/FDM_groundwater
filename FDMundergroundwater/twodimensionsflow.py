import numpy as np
import numpy.linalg as nla
import matplotlib

matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
import sympy as sy


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
        self.name_chinese = "承压含水层稳定二维流"
        self.w = None
        self.S = None
        self.T = None

    def transmissivity(self, T):  # 承压含水层导水系数的设定
        self.T = float(T)

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定，可以设定为x,y的函数。
        self.w = w

    def solve(self):
        # 对于承压含水层二维稳定流，定义两个参数 x y
        x = sy.symbols("x")
        y = sy.symbols("y")
        # X,Y轴差分点的数目
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
                H_ALL[i, j] = H[i * m + j]

        return H_ALL


class Unconfined_aquifer_SF(Stableflow):
    def __init__(self):
        super().__init__()
        self.name_chinese = "潜水含水层稳定二维流"
        self.w = None
        self.K = None
        self.ha = None

    def reference_thickness(self, ha):  # 潜水含水层的参考厚度，使用参考厚度法来简化该方程
        self.ha = float(ha)

    def hydraulic_conductivity(self, K):  # 潜水含水层渗透系数的设定
        self.K = float(K)

    def leakage_recharge(self, w: str = "0"):  # 潜水含水层源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
        self.w = w

    def solve(self):
        # 对于潜水含水层二维稳定流，定义两个参数 x y
        x = sy.symbols("x")
        y = sy.symbols("y")
        # X,Y轴差分点的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1

        # 对函数W(x)为源汇项函数除以渗透系数和参考厚度
        def W(x, y):
            return eval(self.w) / (self.K * self.ha)

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
                H_ALL[i, j] = H[i * m + j] + self.ha

        return H_ALL


class Unstableflow:
    def __init__(self):
        self.yl = None
        self.h_t = None
        self.h_b = None
        self.ic = None
        self.tl = None
        self.st = None
        self.name_chinese = "非稳定二维流"
        self.xl = None
        self.sl = None
        self.h_r = None
        self.h_l = None

    def t_boundary(self, h_t, Dirichlet=True, Neumann=False, Robin=False):  # 上边界
        if Dirichlet:
            self.h_t = float(h_t)

    def b_boundary(self, h_b, Dirichlet=True, Neumann=False, Robin=False):  # 下边界
        if Dirichlet:
            self.h_b = float(h_b)

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

    def y_length(self, yl):  # Y轴轴长
        self.yl = float(yl)

    def t_length(self, tl):  # 时间轴轴长，原则上单位为天
        self.tl = float(tl)

    def initial_condition(self, ic):  # 初始条件的水头设定
        self.ic = float(ic)

    def draw(self, H_ALL: np.ndarray):
        # X轴单元格的数目
        m = int(self.xl / self.sl) + 1
        # y轴单元格数目
        n = int(self.yl / self.sl) + 1
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
        ax.plot_surface(X, T, H_ALL, linewidth=0, antialiased=True, cmap=plt.get_cmap('rainbow'))
        plt.title("差分数值解(差分步长%s)" % self.sl)
        plt.show()


class Confined_aquifer_USF(Unstableflow):
    def __init__(self):
        super().__init__()
        self.name_chinese = "承压含水层非稳定二维流"
        self.w = None
        self.T = None

    def transmissivity(self, T):  # 承压含水层导水系数的设定
        self.T = float(T)

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定，可以设定为x,y,t的函数。
        self.w = w

    def solve(self):
        # 对于承压含水层二维非稳定流，定义两个参数 x y
        x = sy.symbols("x")
        y = sy.symbols("y")
        t = sy.symbols("t")
        # X,Y,t轴差分点的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1
        p = int(self.tl / self.st) + 1

        # 对函数W(x)为源汇项函数除以倒水系数
        def W(x, y, t):
            return eval(self.w) / self.T

        # 创建一个全部值为0的矩阵
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n * p, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n * p, m * n * p))
        # 系数a矩阵行数
        l_a = 0
        while l_a < m * n:
            for k in range(0, p):  # 对时间进行扫描
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
                H_ALL[i, j] = H[i * m + j]

        return H_ALL


if __name__ == "__main__":
    a = Unconfined_aquifer_SF()
    a.x_length(20)
    a.y_length(10)
    a.step_length(1)
    a.hydraulic_conductivity(1)
    a.reference_thickness(10)
    a.l_boundary(1)
    a.r_boundary(1)
    a.t_boundary(3)
    a.b_boundary(-1)
    a.leakage_recharge("0")
    b = a.solve()
    print(b)
    a.draw(b)
