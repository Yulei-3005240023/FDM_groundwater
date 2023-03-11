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
        plt.suptitle(self.name_chinese)
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

    def draw(self, H_ALL: np.ndarray, title=''):
        # X轴单元格的数目
        m = int(self.xl / self.sl) + 1
        # y轴单元格数目
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
        ax.set(zlabel='水头（m）', ylabel='Y轴（m）', xlabel='X轴（m）')
        plt.suptitle(self.name_chinese)
        if title == '':
            plt.title("差分数值解(差分空间步长{0}，时间步长{1})".format(self.sl, self.st))
        else:
            plt.title(title)
        plt.show()


class Confined_aquifer_USF(Unstableflow):  # 承压含水层水流方程
    def __init__(self):
        super().__init__()
        self.S = None
        self.name_chinese = "承压含水层非稳定二维流"
        self.w = None
        self.T = None

    def transmissivity(self, T):  # 承压含水层导水系数的设定
        self.T = float(T)

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定，可以设定为x,y,t的函数。
        self.w = w

    def storativity(self, S):  # 承压含水层贮水系数（弹性给水度）的设定
        self.S = float(S)

    def solve(self):
        # 对于承压含水层二维非稳定流，定义两个参数 x y
        x = sy.symbols("x")
        y = sy.symbols("y")
        t = sy.symbols("t")
        # X,Y,t轴差分点的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1
        p = int(self.tl / self.st) + 1

        # 对函数W(x)为源汇项函数除以导水系数
        def W(x, y, t):
            return eval(self.w) / self.T

        # 创建一个列表来储存所有时刻的水头值
        H_all_time = []
        # 常数b矩阵
        H_b = np.zeros((m * n * p, 1))
        # 系数a矩阵5
        H_a = np.zeros((m * n * p, m * n * p))
        # 系数a矩阵行数
        l_a = 0
        while l_a < m * n * p:
            for k in range(0, p):  # 对时间进行扫描
                for i in range(0, n):  # 对行进行扫描
                    for j in range(0, m):  # 对列进行扫描
                        #  初值设定
                        if k == 0:
                            H_a[l_a, l_a] = 1
                            H_b[l_a] = self.ic
                        else:
                            # 上下左右边界赋值
                            if (i - 1) < 0:
                                H_b[l_a] = H_b[l_a] - self.h_t
                            if (j - 1) < 0:
                                H_b[l_a] = H_b[l_a] - self.h_l
                            if (i + 1) == n:
                                H_b[l_a] = H_b[l_a] - self.h_b
                            if (j + 1) == m:
                                H_b[l_a] = H_b[l_a] - self.h_r
                            # 给位置为(i-1,j,k)处的水头赋上系数值
                            if (i - 1) >= 0:
                                H_a[l_a, k * n * m + (i - 1) * m + j] = self.st
                            # 给位置为(i+1,j,k)处的水头赋上系数值
                            if (i + 1) < n:
                                H_a[l_a, k * n * m + (i + 1) * m + j] = self.st
                            # 给位置为(i,j-1,k)处的水头赋上系数值
                            if (j - 1) >= 0:
                                H_a[l_a, k * n * m + i * m + j - 1] = self.st
                            # 给位置为(i,j+1,k)处的水头赋上系数值
                            if (j + 1) < m:
                                H_a[l_a, k * n * m + i * m + j + 1] = self.st
                            # 给位置为(i,j,k)处的水头赋上系数值
                            H_a[l_a, k * n * m + i * m + j] = -4 * self.st - self.S * self.sl * self.sl / self.T
                            # 给位置为(i,j,k-1)处的水头赋上系数值
                            H_a[l_a, (k - 1) * n * m + i * m + j] = self.S * self.sl * self.sl / self.T
                            # 源汇项赋值
                            H_b[l_a] = H_b[l_a] - W(j * self.sl, i * self.sl, k * self.st) * self.sl * self.sl * self.st
                        l_a += 1
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for k in range(0, p):
            # 创建一个全部值为0的矩阵，用来存放每一个单独时刻的水头值
            H_ALL = np.zeros((n, m))
            for i in range(0, n):  # 对行进行扫描
                for j in range(0, m):  # 对列进行扫描
                    H_ALL[i, j] = H[k * m * n + i * m + j]
            H_all_time.append(H_ALL)

        return H_all_time


class Unconfined_aquifer_USF(Unstableflow):  # 潜水Boussinesq方程
    def __init__(self):
        super().__init__()
        self.K = None
        self.ha = None
        self.Sy = None
        self.name_chinese = "潜水含水层非稳定二维流"
        self.w = None

    def hydraulic_conductivity(self, K):  # 潜水含水层渗透系数的设定
        self.K = float(K)

    def reference_thickness(self, ha):  # 潜水含水层的参考厚度，使用参考厚度法来简化该方程
        self.ha = float(ha)

    def leakage_recharge(self, w: str = "0"):  # 潜水含水层越流补给源汇项的设定，可以设定为x,y,t的函数。
        self.w = w

    def specific_yield(self, Sy):  # 潜水含水层给水度（重力给水度）的设定
        self.Sy = float(Sy)

    def solve(self):
        # 对于潜水含水层二维非稳定流，定义两个参数 x y
        x = sy.symbols("x")
        y = sy.symbols("y")
        t = sy.symbols("t")
        # X,Y,t轴差分点的数目
        m = int(self.xl / self.sl) + 1
        n = int(self.yl / self.sl) + 1
        p = int(self.tl / self.st) + 1

        # 对函数W(x)为源汇项函数除以导水系数
        def W(x, y, t):
            if self.w == '0':
                return 0
            else:
                return eval(self.w) / self.K

        # 创建一个列表来储存所有时刻的水头值
        H_all_time = []
        l_a = 0
        while l_a < m * n * p:
            for k in range(0, p):  # 对时间进行扫描
                H_c = np.zeros((m * n, m * n))
                l_c = 0
                H_d = np.zeros((m * n, 1))
                for i in range(0, n):  # 对行进行扫描
                    for j in range(0, m):  # 对列进行扫描
                        #  初值设定
                        if k == 0:
                            H_c[l_c, l_c] = 1
                            H_d[l_c] = self.ic
                        else:
                            # 上下左右边界赋值
                            if (i - 1) < 0:
                                H_d[l_c] = H_d[l_c] - self.h_t * (H_all_time[k - 1][i, j] + self.h_t)
                                tt = self.h_t
                            else:
                                tt = H_all_time[k - 1][i - 1, j]
                            if (j - 1) < 0:
                                H_d[l_c] = H_d[l_c] - self.h_l * (H_all_time[k - 1][i, j] + self.h_l)
                                ll = self.h_l
                            else:
                                ll = H_all_time[k - 1][i, j - 1]
                            if (i + 1) == n:
                                H_d[l_c] = H_d[l_c] - self.h_b * (H_all_time[k - 1][i, j] + self.h_b)
                                bb = self.h_b
                            else:
                                bb = H_all_time[k - 1][i + 1, j]
                            if (j + 1) == m:
                                H_d[l_c] = H_d[l_c] - self.h_r * (H_all_time[k - 1][i, j] + self.h_r)
                                rr = self.h_r
                            else:
                                rr = H_all_time[k - 1][i, j + 1]
                            # 给位置为(i-1,j,k)处的水头赋上系数值
                            if (i - 1) >= 0:
                                H_c[l_c, (i - 1) * m + j] = H_all_time[k - 1][i, j] + H_all_time[k - 1][i - 1, j]
                            # 给位置为(i+1,j,k)处的水头赋上系数值
                            if (i + 1) < n:
                                H_c[l_c, (i + 1) * m + j] = H_all_time[k - 1][i, j] + H_all_time[k - 1][i + 1, j]
                            # 给位置为(i,j-1,k)处的水头赋上系数值
                            if (j - 1) >= 0:
                                H_c[l_c, i * m + j - 1] = H_all_time[k - 1][i, j] + H_all_time[k - 1][i, j - 1]
                            # 给位置为(i,j+1,k)处的水头赋上系数值
                            if (j + 1) < m:
                                H_c[l_c, i * m + j + 1] = H_all_time[k - 1][i, j] + H_all_time[k - 1][i, j + 1]
                            # 给位置为(i,j,k)处的水头赋上系数值
                            H_c[l_c, l_c] = -4 * H_all_time[k - 1][i, j] - rr - ll - bb - tt - 2 * self.sl * self.Sy / (self.K * self.st)
                            # 源汇项赋值
                            H_d[l_c] = H_d[l_c] - W(j * self.sl, i * self.sl, k * self.st) * self.sl * self.sl - H_all_time[k - 1][i, j] * (2 * self.sl * self.Sy / (self.K * self.st))
                        l_c += 1
                        l_a += 1
                # 解矩阵方程
                H = nla.solve(H_c, H_d)
                # 创建一个全部值为0的矩阵，用来存放每一个单独时刻的水头值
                H_ALL = np.zeros((n, m))
                for i in range(0, n):  # 对行进行扫描
                    for j in range(0, m):  # 对列进行扫描
                        H_ALL[i, j] = H[i * m + j]
                H_all_time.append(H_ALL)

        return H_all_time


class Toth_difficult_baisn:  # Tóth名字里面有非ASCII字符，故给他改个名字
    def __init__(self):
        self.H = None
        self.L = None
        self.h0 = None
        # 对于Tóth复杂盆地的多年平均水位方程，定义参数 x,z
        x = sy.symbols("x")
        z = sy.symbols("z")

    def basin_length(self, L):  # 设置盆地宽度
        self.L = float(L)

    def basin_high(self, H):  # 设置河谷高程
        self.H = float(H)

    def average_water_level_equation(self, h0: str):  # Tóth复杂盆地的多年平均水位方程
        # 在Tóth的假设下为x,z的函数，前缀需要带np.如：z0 + x * np.tan(sy.pi/3) +np.pi/3 * np.sin(2 * np.pi * x)
        self.h0 = h0

    def draw_water_level(self):  # 根据多年平均水位方程还绘制水位
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False

        # 对于多年平均水位方程
        def H0(x):
            return eval(self.h0)

        # X轴
        x = np.linspace(0, self.L, 200)
        # Z轴
        Z = H0(x)
        plt.plot(x, Z, linestyle='-', linewidth=1, antialiased=True)
        plt.title("Tóth复杂盆地的多年平均水位")
        plt.show()


if __name__ == "__main__":
    a = Toth_difficult_baisn()
    a.h0 = '5000 + x * np.tan(np.pi/200) +3 * np.sin((2 * np.pi * x) / 6 * np.cos(np.pi/200))/np.cos(np.pi/200)'
    a.L = 20000
    a.draw_water_level()
