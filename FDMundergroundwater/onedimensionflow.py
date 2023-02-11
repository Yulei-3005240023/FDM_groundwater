from multiprocessing import Process, Manager
import numpy as np
import sympy as sy
from sympy import symbols
import numpy.linalg as nla
import matplotlib

matplotlib.use('QtAgg')
import matplotlib.pyplot as plt


class Stableflow:

    def __init__(self):
        self.name_chinese = "稳定一维流"
        self.l = None
        self.sl = None
        self.h_r = None
        self.h_l = None

    def l_boundary(self, h_l, Dirichlet=True, Neumann=False, Robin=False):  # 左边界 或者潜水含水层相对于参考厚度的水头
        if Dirichlet:
            self.h_l = float(h_l)

    def r_boundary(self, h_r, Dirichlet=True, Neumann=False, Robin=False):  # 右边界 或者潜水含水层相对于参考厚度的水头
        if Dirichlet:
            self.h_r = float(h_r)

    def step_length(self, sl):  # 差分步长
        self.sl = float(sl)

    def length(self, l):  # 轴长
        self.l = float(l)

    def draw(self, H_ALL):
        # X轴单元格的数目
        m = int(self.l / self.sl) + 1
        # X轴
        X = np.linspace(0, self.l, m)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot()
        ax.plot(X, H_ALL, linewidth=1, antialiased=True)

        def maxH_y(h_all):
            hy = 0
            for i in h_all:
                if i > hy:
                    hy = i
            return hy

        def minH_y(h_all):
            hy = 0
            for i in h_all:
                if i < hy:
                    hy = i
            return hy

        ax.set_ylim(minH_y(H_ALL), maxH_y(H_ALL))
        plt.suptitle(self.name_chinese)
        plt.title("差分数值解(差分步长%s)" % self.sl)
        plt.show()


class Confined_aquifer_SF(Stableflow):

    def __init__(self):
        super().__init__()
        self.name_chinese = "承压含水层稳定一维流"
        self.w = None
        self.t = None

    def transmissivity(self, t: str = "1"):  # 承压含水层导水系数的设定，可以为一个常数也可以为带有前缀为sy.的函数，如sy.sin(x)
        self.t = str(t)

    def leakage_recharge(self, w: str = "0"):  # 承压含水层越流补给源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
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
        m = int(self.l / self.sl) + 1
        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros(m)
        # 常数b矩阵w
        H_b = np.zeros((m, 1))
        # 系数a矩阵
        H_a = np.zeros((m, m))
        for i in range(0, m):  # 对列进行扫描
            # 左右边界赋值
            if (i - 1) < 0:
                H_a[i, i] = 1
                H_b[i] = self.h_l
            elif (i + 1) == m:
                H_a[i, i] = 1
                H_b[i] = self.h_r
            else:
                # 源汇项赋值
                H_b[i] = H_b[i] - W(i * self.sl) * (self.sl * self.sl)
                # 给位置为(i-1)处的水头赋上系数值
                H_a[i, (i - 1)] = T(i * self.sl) - (self.sl * sy.diff(T(x), x).subs(x, i * self.sl))
                # 给位置为(i+1)处的水头赋上系数值
                H_a[i, (i + 1)] = T(i * self.sl)
                # 给位置为(i)处的水头赋上系数值
                H_a[i, i] = self.sl * sy.diff(T(x), x).subs(x, i * self.sl) - 2 * T(i * self.sl)
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for i in range(0, m):
            H_ALL[i] = H[i]
        return H_ALL


class Unconfined_aquifer_SF(Stableflow):
    def __init__(self):
        super().__init__()
        self.w = None
        self.K = None
        self.ha = None

    def reference_thickness(self, ha):  # 潜水含水层的参考厚度，使用参考厚度法来简化该方程
        self.ha = float(ha)

    def hydraulic_conductivity(self, K: str):  # 潜水含水层渗透系数的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
        self.K = K

    def leakage_recharge(self, w: str = "0"):  # 潜水含水层源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
        self.w = w

    def solve(self):
        # 对于潜水含水层一维稳定流，定义参数 x
        x = symbols("x")

        # 对于K(x)，定义为渗透系数K乘以参考厚度ha
        def K(x):
            return eval(self.K) * self.ha

        # 对源汇项函数W(x)
        def W(x):
            return eval(self.w)

        # X轴差分点的数目
        m = int(self.l / self.sl) + 1
        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros(m)
        # 常数b矩阵w
        H_b = np.zeros((m, 1))
        # 系数a矩阵
        H_a = np.zeros((m, m))
        for i in range(0, m):  # 对列进行扫描
            # 源汇项赋值
            H_b[i] = H_b[i] - W(i * self.sl) * (self.sl * self.sl)
            # 左右边界赋值
            if (i - 1) < 0:
                H_b[i] = H_b[i] - (K(i * self.sl) - self.sl * (sy.diff(K(x), x).subs(x, i * self.sl))) * self.h_l
            if (i + 1) == m:
                H_b[i] = H_b[i] - K(i * self.sl) * self.h_r
            # 给位置为(i-1)处的水头赋上系数值
            if (i - 1) >= 0:
                H_a[i, (i - 1)] = K(i * self.sl) - (self.sl * sy.diff(K(x), x).subs(x, i * self.sl))
            # 给位置为(i+1)处的水头赋上系数值
            if (i + 1) < m:
                H_a[i, (i + 1)] = K(i * self.sl)
            # 给位置为(i)处的水头赋上系数值
            H_a[i, i] = self.sl * sy.diff(K(x), x).subs(x, i * self.sl) - 2 * K(i * self.sl)
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for i in range(0, m):
            H_ALL[i] = H[i] + self.ha
        return H_ALL


class Unstableflow:
    def __init__(self):
        self.ic = None
        self.tl = None
        self.st = None
        self.name_chinese = "非稳定一维流"
        self.xl = None
        self.sl = None
        self.h_r = None
        self.h_l = None
        tl = self.tl

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

    def initial_condition(self, ic):  # 初始条件的水头设定
        self.ic = float(ic)

    def draw(self, H_ALL: np.ndarray, title=''):
        # X轴单元格的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴单元格数目
        n = int(self.tl / self.st) + 1
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
        plt.suptitle(self.name_chinese)
        if title == '':
            plt.title("差分数值解(差分空间步长{0}，时间步长{1})".format(self.sl, self.st))
        else:
            plt.title(title)
        plt.show()


def solve_as(a, xl, sl, n_, t_, h_l, h_r, ic, c, T, w, fourier_series, return_dict):  # 该函数用于多核心的傅里叶解析解计算
    # X轴差分点的数目
    m = int(xl / sl) + 1

    def M(x, t):
        h = (1 - x / xl)
        for i_ in range(1, fourier_series + 1):
            h -= 2 * np.sin(i_ * np.pi * x / xl) / (i_ * np.pi) * np.exp(
                -(a * i_ * i_ * np.pi * np.pi) / (xl * xl) * t)
        return h

    def N(x):
        if T == 0:
            return 0
        else:
            h = -1 * (float(w) / (2 * T)) * x * x + ((h_r - h_l) / xl + float(w) * xl / (2 * T)) * x + h_l
        return h

    # X轴
    X = np.linspace(0, xl, m)
    # 所有解值矩阵
    H = np.zeros((m * n_, 1))
    # 定义解值矩阵的行数
    l = 0
    while l < m * n_:
        for k in t_:
            for j in X:
                H[l] = (h_l - ic) * M(j, k) + (h_r - ic) * M(xl - j, k)  # + N(j)
                l += 1

    return_dict[c] = H
    return return_dict


class Confined_aquifer_USF(Unstableflow):
    def __init__(self):
        super().__init__()
        self.S = None
        self.T = None
        self.name_chinese = "承压含水层一维非稳定流"
        self.w = None
        self.a = None

    def transmissivity(self, T):  # 承压含水层导水系数的设定
        self.T = float(T)

    def storativity(self, S):  # 承压含水层贮水系数（弹性给水度）的设定
        self.S = float(S)

    def pressure_diffusion_coefficient(self, a_):  # 承压含水层压力扩散系数的设定，等于导水系数除以贮水系数T/S
        self.a = float(a_)

    def leakage_recharge(self, w: str = "0"):
        # 承压含水层越流补给源汇项的设定,承压含水层越流补给源汇项的设定，在解析解的条件下为一个常数，在数值解的条件下也可以为带有前缀为sy.的函数,如sy.sin(x) + sy.cos(t)
        self.w = w

    def solve(self):  # 使用有限差分法对该非稳定流方程进行求解
        # 如果未设定压力扩散系数
        if self.a is None or self.a == '':
            self.a = self.T / self.S
        # 对于承压含水层一维非稳定流，定义两个参数 x t
        x = symbols("x")
        t = symbols("t")
        # X轴差分点的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴差分点的数目
        n = int(self.tl / self.st) + 1

        # 对函数W(x)为源汇项函数除以导水系数
        def W(x, t):
            return eval(self.w) / self.T

        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n, m * n))
        # 定义系数a矩阵的行数
        l_a = 0
        while l_a < m * n:
            for k in range(0, n):  # 对行(时间轴)进行扫描
                for i in range(0, m):  # 对列(X轴)进行扫描
                    # 时间边界赋值(初始条件）
                    if k == 0:
                        H_a[l_a, l_a] = 1
                        H_b[l_a] = self.ic
                    else:
                        # 源汇项赋值
                        H_b[l_a] = H_b[l_a] - self.a * (self.sl * self.sl) * self.st * W(i * self.sl, k * self.st)
                        # 左右边界赋值
                        if (i - 1) < 0:
                            H_b[l_a] = H_b[l_a] - (self.a * self.st) * self.h_l
                        if (i + 1) == m:
                            H_b[l_a] = H_b[l_a] - self.a * self.st * self.h_r
                        # 给位置为(i, k)处的水头赋上系数值
                        H_a[l_a, l_a] = H_a[l_a, l_a] - 2 * self.a * self.st - self.sl * self.sl
                        # 给位置为(i-1，k)处的水头赋上系数值
                        if (i - 1) >= 0:
                            H_a[l_a, l_a - 1] = H_a[l_a, l_a - 1] + self.a * self.st
                        # 给位置为(i+1, k)处的水头赋上系数值
                        if (i + 1) < m:
                            H_a[l_a, l_a + 1] = H_a[l_a, l_a + 1] + self.a * self.st
                        # 给位置为(i, k-1)处的水头赋上系数值
                        if (k - 1) >= 0:
                            H_a[l_a, l_a - m] = H_a[l_a, l_a - m] + self.sl * self.sl
                    l_a += 1
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for k in range(0, n):  # 对时间进行扫描
            for i in range(0, m):  # 对空间进行扫描
                H_ALL[k, i] = H[k * m + i]
        return H_ALL

    def solve_analytic_solution(self, fourier_series=20):  # 默认取傅里叶级数前20项，单进程运行函数
        # 如果未设定压力扩散系数
        if self.a is None or self.a == '':
            self.a = self.T / self.S
        # X轴差分点的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴差分点的数目
        n = int(self.tl / self.st) + 1

        def M(x, t):
            h = (1 - x / self.xl)
            for i_ in range(1, fourier_series + 1):
                h -= 2 * np.sin(i_ * np.pi * x / self.xl) / (i_ * np.pi) * np.exp(
                    -(self.a * i_ * i_ * np.pi * np.pi) / (self.xl * self.xl) * t)
            return h

        def N(x):
            h = -1 * float(self.w) / (2 * self.T) * x * x + (
                        ((self.h_r - self.ic) - (self.h_l - self.ic)) / self.xl + float(self.w) * self.xl / (
                            2 * self.T)) * x + self.h_l - self.ic
            return h

        # X轴
        X = np.linspace(0, self.xl, m)
        # T轴
        T = np.linspace(0, self.tl, n)
        # 创建一个全部值为0的矩阵，用于存放各个与差分位置对等的解析解的水头值
        H_ALL = np.zeros((n, m))
        # 所有解值矩阵
        H = np.zeros((m * n, 1))
        # 定义解值矩阵的行数
        l = 0
        while l < m * n:
            for k in T:
                for j in X:
                    H[l] = (self.h_l - self.ic) * M(j, k) + (self.h_r - self.ic) * M(self.xl - j, k) #+ N(j)
                    l += 1

        for k in range(0, n):  # 对时间进行扫描
            for i in range(0, m):  # 对空间进行扫描
                H_ALL[k, i] = H[k * m + i] + self.ic
        return H_ALL

    def solve_multi(self, cpu_cores=2, fourier_series=20):  # 默认取傅里叶级数前20项，CPU运算核心为两个，多核运行选择
        # 如果未设定压力扩散系数
        if self.a is None or self.a == '':
            self.a = self.T / self.S
        # X轴差分点的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴差分点的数目
        n = int(self.tl / self.st) + 1
        # 分配到每一个cpu上的数目
        n_ = n // cpu_cores
        # 创建一个全部值为0的矩阵，用于存放各个与差分位置对等的解析解的水头值
        H_ALL = np.zeros((n, m))
        # 所有解值矩阵
        H = np.zeros((m * n, 1))
        # 创建类似字典的跨进程共享对象
        manager = Manager()
        return_dict = manager.dict()
        # 创建多进程工作列表
        plist = []
        # 多进程任务分配
        for c in range(0, cpu_cores):
            # 时间轴（以拆分）
            if c != (cpu_cores - 1):
                t_ = np.linspace(c * n_ * self.st, (c + 1) * (n_ - 1) * self.st, n_)
            else:
                t_ = np.linspace((c * n_) * self.st, self.tl, n - c * n_)
                n_ = n - c * n_
            p = Process(target=solve_as, args=(
                self.a, self.xl, self.sl, n_, t_, self.h_l, self.h_r, self.ic, c, self.T, self.w, fourier_series,
                return_dict))
            p.start()
            plist.append(p)

        # 多进程运算开始
        for p in plist:
            p.join()
        l = 0
        while l < m * n:
            for c in range(0, cpu_cores):
                for h in return_dict[c]:
                    H[l] = h
                    l += 1
        for k in range(0, n):  # 对时间进行扫描
            for i in range(0, m):  # 对空间进行扫描
                H_ALL[k, i] = H[k * m + i] + self.ic
        return H_ALL


class Unconfined_aquifer_USF(Unstableflow):
    def __init__(self):
        super().__init__()
        self.Sy = None
        self.K = None
        self.w = None
        self.a = None
        self.ha = None

    def reference_thickness(self, ha):  # 潜水含水层的参考厚度，使用参考厚度法来简化该方程
        self.ha = float(ha)

    def pressure_diffusion_coefficient(self, a):  # 潜水含水层压力扩散系数的设定。等于渗透系数乘初始水头常数除给水度Kh0/Sy
        self.a = float(a)

    def leakage_recharge(self, w: str):  # 潜水含水层源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x)
        self.w = w

    def hydraulic_conductivity(self, K):  # 潜水含水层渗透系数的设定
        self.K = float(K)

    def storativity(self, Sy):  # 潜水含水层储水系数（重力给水度）的设定
        self.Sy = float(Sy)

    def solve(self):
        # 如果未设定压力扩散系数
        if self.a is None or self.a == "":
            self.a = (self.K * self.ha) / self.Sy
        # 对于潜水含水层一维非稳定流，定义两个参数 x t
        x = symbols("x")
        t = symbols("t")
        # X轴差分点的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴差分步长
        self.st = self.sl * self.tl / self.xl
        # 时间轴差分点的数目
        n = int(self.tl / self.st) + 1

        # 对函数W(x)为源汇项函数除以K,ha
        def W(x):
            return eval(self.w) / (self.K * self.ha)

        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n, m * n))
        # 定义系数a矩阵的行数
        l_a = 0
        while l_a < m * n:
            for k in range(0, n):  # 对行(时间轴)进行扫描
                for i in range(0, m):  # 对列(X轴)进行扫描
                    # 时间边界赋值(初始条件）
                    if k == 0:
                        H_a[l_a, l_a] = 1
                        H_b[l_a] = self.ic
                    else:
                        # 源汇项赋值
                        H_b[l_a] = H_b[l_a] - self.a * (self.sl * self.sl) * self.st * W(i * self.sl)
                        # 左右边界赋值
                        if (i - 1) < 0:
                            H_b[l_a] = H_b[l_a] - (self.a * self.st) * self.h_l
                        if (i + 1) == m:
                            H_b[l_a] = H_b[l_a] - self.a * self.st * self.h_r
                        # 给位置为(i, k)处的水头赋上系数值
                        H_a[l_a, l_a] = H_a[l_a, l_a] - 2 * self.a * self.st - self.sl * self.sl
                        # 给位置为(i-1，k)处的水头赋上系数值
                        if (i - 1) >= 0:
                            H_a[l_a, l_a - 1] = H_a[l_a, l_a - 1] + self.a * self.st
                        # 给位置为(i+1, k)处的水头赋上系数值
                        if (i + 1) < m:
                            H_a[l_a, l_a + 1] = H_a[l_a, l_a + 1] + self.a * self.st
                        # 给位置为(i, k-1)处的水头赋上系数值
                        if (k - 1) >= 0:
                            H_a[l_a, l_a - m] = H_a[l_a, l_a - m] + self.sl * self.sl
                    l_a += 1
        # 解矩阵方程
        H = nla.solve(H_a, H_b)
        for k in range(0, n):  # 对时间进行扫描
            for i in range(0, m):  # 对空间进行扫描
                H_ALL[k, i] = H[k * m + i] + self.ha
        return H_ALL
