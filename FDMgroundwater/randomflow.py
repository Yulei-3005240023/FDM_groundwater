# 地下水流随机方程的数值模拟求解尝试
import numpy as np
import numpy.linalg as nla
from sympy import symbols
import sympy as sy
import matplotlib

matplotlib.use('QtAgg')
import matplotlib.pyplot as plt


class Random_flow:
    def __init__(self):
        self.ic = None
        self.tl = None
        self.st = None
        self.name_chinese = "非稳定随机一维流"
        self.xl = None
        self.sl = None
        self.h_r = None
        self.h_l = None
        self.B = 1  # 默认一维流的宽度为1个单位

    def l_boundary(self, h_l, Dirichlet=True, Neumann=False, Robin=False):  # 左边界
        if Dirichlet:
            self.h_l = [1, float(h_l)]
        elif Neumann:
            self.h_l = [2, float(h_l)]

    def r_boundary(self, h_r, Dirichlet=True, Neumann=False, Robin=False):  # 右边界
        if Dirichlet:
            self.h_r = [1, float(h_r)]
        elif Neumann:
            self.h_l = [2, float(h_r)]

    def step_length(self, sl):  # X轴差分步长
        self.sl = float(sl)

    def step_time(self, st):  # 时间轴差分步长
        self.st = float(st)

    def x_length(self, xl):  # X轴轴长
        self.xl = float(xl)

    def t_length(self, tl):  # 时间轴轴长，原则上单位为天
        self.tl = float(tl)

    def initial_condition(self, ic: str):  # 初始条件的水头设定
        self.ic = str(ic)

    def width(self, B):  # 含水层宽度的设定
        self.B = float(B)

    def draw(self, H_ALL: np.ndarray, time=0, title=''):  # 按给定的时刻绘制水头曲线
        # X轴单元格的数目
        m = int(self.xl / self.sl) + 1
        # X轴
        X = np.linspace(0, self.xl, m)
        # 可以plt绘图过程中中文无法显示的问题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        # 解决负号为方块的问题
        plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot()
        ax.plot(X, H_ALL[time], linewidth=1, antialiased=True)

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

        ax.set_ylim(minH_y(H_ALL[time]), maxH_y(H_ALL[time]))
        ax.set(ylabel='水头（m）', xlabel='X轴（m）')
        plt.suptitle(self.name_chinese)
        if title == '':
            plt.title("差分数值解，当前为第{0}时刻(差分空间步长{1}，时间步长{2})".format(time, self.sl, self.st))
        else:
            plt.title(title)
        plt.show()


class Random_one_dimension_boussinesq(Random_flow):
    def __init__(self):
        super().__init__()
        self.Sy = None
        self.K = None
        self.w = None
        self.a = None
        self.a_as = None
        self.ha = None
        self.name_chinese = '潜水含水层随机非稳定一维流'

    def reference_thickness(self, ha):  # 潜水含水层的参考厚度，解析解求解中使用参考厚度法线性化偏微分方程
        self.ha = float(ha)

    def pressure_diffusion_coefficient(self, a):  # 潜水含水层压力扩散系数的设定。等于渗透系数乘初始水头常数除给水度Kh0/Sy
        self.a = float(a)

    def leakage_recharge(self, w: str):  # 潜水含水层源汇项的设定，可以为一个常数也可以为带有前缀为sy.的函数,如sy.sin(x) + sy.cos(y)
        self.w = w

    def hydraulic_conductivity(self, K):  # 潜水含水层渗透系数的设定
        self.K = float(K)

    def specific_yield(self, Sy):  # 潜水含水层储水系数（重力给水度）的设定
        self.Sy = float(Sy)

    def solve(self):
        # 如果未设定压力扩散系数
        if self.a is None or self.a == "":
            self.a = self.K / self.Sy
        # 对于潜水含水层一维非稳定流，定义两个参数 x t
        x = symbols("x")
        t = symbols("t")
        # X轴差分点的数目
        m = int(self.xl / self.sl) + 1
        # 时间轴差分点的数目
        n = int(self.tl / self.st) + 1

        # 对函数W(x, t)定义为源汇项函数除以渗透系数K
        def W(x, t):
            return eval(self.w) / self.K

        # 函数IC定义为初始水头分布曲线
        def IC(x):
            return eval(self.ic)

        # 创建一个全部值为0的矩阵，用于存放各个差分位置的水头值
        H_ALL = np.zeros((n, m))
        # 常数b矩阵
        H_b = np.zeros((m * n, 1))
        # 系数a矩阵
        H_a = np.zeros((m * n, m * n))
        # 定义系数a矩阵的行数

        # 矩阵赋值
        for k in range(0, n):  # 对行(时间轴)进行扫描
            iteration_times = 0  # 迭代运算次数计数
            H_previous_iteration = np.zeros((1, m))
            # 迭代运算开始
            while True:
                H_a = np.zeros((m, m))
                l_a = 0
                H_b = np.zeros((m, 1))

                if iteration_times == 0 and k != 0:
                    H_previous_iteration = H_ALL[k - 1]  # 前次迭代的当前时刻水头数值,此处未开始计算，使用上一时刻的水头值进行近似

                for i in range(0, m):  # 对列(X轴)进行扫描
                    # 时间边界赋值(初始条件）
                    if k == 0:
                        H_a[l_a, l_a] = 1
                        H_b[l_a] = IC(i * self.sl)

                    # 左边界赋值
                    elif (i - 1) < 0 and self.h_l[0] == 1:  # 一类边界判断
                        H_a[l_a, l_a] = 1
                        H_b[l_a] = self.h_l[1]
                    elif (i - 1) < 0 and self.h_l[0] == 2:  # 二类边界判断
                        # 源汇项赋值
                        H_b[l_a] = H_b[l_a] - W(i * self.sl, k * self.st) - self.Sy / (self.K * self.st) * H_ALL[
                            k - 1, i] - 2 * self.sl * self.h_l[1] * H_previous_iteration[i] / (
                                           self.sl * self.sl)
                        # 给位置为(i, k)处的水头赋上系数值
                        H_a[l_a, l_a] = -(H_previous_iteration[i + 1] + H_previous_iteration[i]) / (
                                2 * self.sl * self.sl) - H_previous_iteration[i] / (
                                                self.sl * self.sl) - self.Sy / (self.K * self.st)
                        # 给位置为(i+1, k)处的水头赋上系数值
                        H_a[l_a, l_a + 1] = (H_previous_iteration[i + 1] + H_previous_iteration[i]) / (
                                2 * self.sl * self.sl) + H_previous_iteration[i] / (
                                                    self.sl * self.sl)

                    # 右边界赋值
                    elif (i + 1) == m and self.h_r[0] == 1:
                        H_a[l_a, l_a] = 1
                        H_b[l_a] = self.h_r[1]
                    elif (i + 1) == m and self.h_r[0] == 2:
                        # 源汇项赋值
                        H_b[l_a] = H_b[l_a] - W(i * self.sl, k * self.st) - self.Sy / (self.K * self.st) * H_ALL[
                            k - 1, i] + 2 * self.sl * self.h_r[1] * H_previous_iteration[i] / (
                                             self.sl * self.sl)
                        # 给位置为(i, k)处的水头赋上系数值
                        H_a[l_a, l_a] = -H_previous_iteration[i]/ (
                                2 * self.sl * self.sl) - (H_previous_iteration[i] + H_previous_iteration[i - 1]) / (
                                                2 * self.sl * self.sl) - self.Sy / (self.K * self.st)
                        # 给位置为(i-1, k)处的水头赋上系数值
                        H_a[l_a, l_a - 1] = (H_previous_iteration[i] + H_previous_iteration[i - 1]) / (
                                2 * self.sl * self.sl) + H_previous_iteration[i] / (
                                                    2 * self.sl * self.sl)
                    else:  # 非边界部分赋值
                        # 源汇项赋值
                        H_b[l_a] = H_b[l_a] - W(i * self.sl, k * self.st) - self.Sy / (self.K * self.st) * H_ALL[
                            k - 1, i]
                        # 给位置为(i, k)处的水头赋上系数值
                        H_a[l_a, l_a] = -(H_previous_iteration[i + 1] + H_previous_iteration[i]) / (
                                2 * self.sl * self.sl) - (H_previous_iteration[i] + H_previous_iteration[i - 1]) / (
                                                2 * self.sl * self.sl) - self.Sy / (self.K * self.st)
                        # 给位置为(i-1，k)处的水头赋上系数值
                        H_a[l_a, l_a - 1] = (H_previous_iteration[i] + H_previous_iteration[i - 1]) / (
                                2 * self.sl * self.sl)
                        # 给位置为(i+1, k)处的水头赋上系数值
                        H_a[l_a, l_a + 1] = (H_previous_iteration[i + 1] + H_previous_iteration[i]) / (
                                2 * self.sl * self.sl)
                    l_a += 1

                H = nla.solve(H_a, H_b)  # 进行当前时刻的水头计算结果
                if k == 0:  # 第零时刻不参与迭代计算
                    break

                # 判断是否满足精度需求
                precision = 0
                for u in range(0, m):
                    if abs(H_previous_iteration[u] - H[u]) > 0.001:
                        precision = 1
                if precision != 1:
                    break
                else:
                    iteration_times += 1
                    H_previous_iteration = H

                if iteration_times >= 100:
                    break
            for o in range(0, m):  # 对空间进行扫描，整合成所有适合的计算水头
                H_ALL[k, o] = H[o]
        return H_ALL


if __name__ == "__main__":
    flow = Random_one_dimension_boussinesq()
    flow.sl = 1
    flow.st = 1
    flow.ic = '10'
    flow.tl = 10
    flow.xl = 100
    flow.h_r = [1, 20]
    flow.h_l = [2, 0]
    flow.Sy = 0.08
    flow.K = 0.5
    flow.w = '0.095'
    h = flow.solve()
    flow.draw(H_ALL=h, time=6)
