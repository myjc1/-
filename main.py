# -*- coding = utf-8 -*-
# @Time:  17:26
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：main.py
# @Software: PyCharm
from readinformation import *
from numericalmethod import *
dt, number, a_ground = read_peer_motion('ELC-NS.AT2')
M, K, C, hyst_model, nodes, p = read_structure_information('information1.txt')
NMDOF(M, K, C, dt, number, a_ground, hyst_model, nodes, p).newmark_beta()
# MDOF(M, K, C, dt, number, a_ground, nodes)
# NSDOF(a_ground, dt, number, 1, 1, 2.27, 0.01, 'bilinear model')
# SDOF(a_ground, dt, 1, 0.05, number).newmark_beta()
