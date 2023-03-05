#include<math.h>
//声明使用语言为C语言
extern"C"{

//傅里叶级数求解函数_参考厚度法
double M_rtm(double a, double x, double t, double xl, int fourier_series=1000){
    double h = (1-x/xl);
    for (int i_ = 1; i_<=fourier_series; i_++)
        h -= 2 * sin(i_ * 3.1415926 * x / xl) / (i_ * 3.1415926) * exp(-(a * i_ * i_ * 3.1415926 * 3.1415926) / (xl * xl) * t);
    return h;
}

//傅里叶级数求解函数_平方法
double M_sm(double a, double x, double t, double xl, double reference_thickness, double h_l, double h_r, int fourier_series=1000){
    double f;
    for (int i = 1; i<=fourier_series; i++){
        if((fourier_series % 2)==0){
            f += (2 / (i * 3.1415926)) * (reference_thickness * reference_thickness - 0.5 * (h_l * h_l + h_r * h_r)) * sin(i * 3.1415926 * x / xl) * exp(-(a * i * i * 3.1415926 * 3.1415926) / (xl * xl) * t);
        }
        else{
            f += (2 / (i * 3.1415926)) * (0.5 * (-h_l * h_l + h_r * h_r)) * sin(i * 3.1415926 * x / xl) * exp(-(a * i * i * 3.1415926 * 3.1415926) / (xl * xl) * t);
        }
    }
    return f;
}

//使用参考厚度法求解一维潜水含水层非稳定流解析解函数，无源汇项，注意在此条件下，初始时刻水头与参考厚度相等，所以仅有参考厚度输入
double Boussinesq_one_dimension_unstable_flow_reference_thickness_method(double a, double x, double t, double xl, double h_l, double h_r, double reference_thickness, int fourier_series=1000){
    double h;
    h = (h_l - reference_thickness) * M_rtm(a, x, t, xl, fourier_series) + (h_r - reference_thickness) * M_rtm(a, xl-x, t, xl, fourier_series);//《地下水运动方程》38页（王旭升）
    return h + reference_thickness;
}

//使用平方法求解一维潜水含水层非稳定流解析解函数，无源汇项，注意在此条件下，初始时刻水头与参考厚度相等，所以仅有参考厚度输入
double Boussinesq_one_dimension_unstable_flow_square_method(double a, double x, double t, double xl, double h_l, double h_r, double reference_thickness, int fourier_series=1000){
    double h;
    double g;
    g = M_sm(a, x, t, xl, reference_thickness, h_l, h_r, fourier_series) + 0.5 * h_l * h_l + (x / xl) * (h_r * h_r - h_l * h_l);
    h = sqrt(2 * g);
    return h;
}
}
