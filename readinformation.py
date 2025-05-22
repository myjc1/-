# -*- coding = utf-8 -*-
# @Time:  16:23
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：readinformation.py
# @Software: PyCharm
import re
import matrix
import numpy as np


def read_peer_motion(filename):  # peer地震动的数据读取，输出单位为m/s^2
    file = open(filename, 'r+')
    lines = file.read().splitlines()  # 读取AT2文件中所有的数据
    file.close()
    npts = int(re.findall(r'\d+', lines[3])[0])  # 加速度数据的个数
    step = float(re.findall(r'\.\d+', lines[3])[0])  # 加速度的时间间隔
    ag = []  # 初始化地面加速度
    for line in range(4, len(lines)):
        ag_str = lines[line].split('  ')
        b = list(filter(lambda x: x != '' and x != ' ', ag_str))
        for j in b:
            ag.append(float(j) * 9.807)  # 将加速度单位换算为m/s^2
    return step, npts, ag


def read_structure_information(information):
    """
    example:
    NDOF:	3
    ksai:	0.05
    Hyst:   UHYST02
    Mas1:	1000
    Mas2:	1000
    Mas3:	1000
    Spr1:	2.2e6	2.646e7	0.05
    Spr2:	2.2e6	1.205e7	0.05
    Spr3:	2.2e6	0.323e7	0.05
    :return:
    """
    inf = open(information, 'r+')
    lines = inf.read().splitlines()
    inf.close()
    m = []
    nds = int(re.findall(r'\d+', lines[0])[0])
    kesai = float(re.findall(r':(.+)', lines[1])[0])
    model = re.findall(r':(.*)', lines[2])[0]
    props = np.zeros([nds, 3])  # k,fy,eta
    for i in range(3, 3+nds):
        m.append(int(re.findall(r':(.*)', lines[i])[0]))
        prop = re.findall(r':(.*)', lines[i+nds])[0].split(',')
        props[i - 3, 0] = float(prop[0])
        props[i - 3, 1] = float(prop[1])
        props[i - 3, 2] = float(prop[2])
    m_matrix = np.diag(m)
    k_matrix = matrix.kmatrix(props[:, 0])
    c_matrix = matrix.cmatrix(m_matrix, k_matrix, kesai)  # 质量矩阵，刚度矩阵，阻尼矩阵
    return m_matrix, k_matrix, c_matrix, model, nds, props
