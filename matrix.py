# -*- coding = utf-8 -*-
# @Time:  16:00
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：matrix.py
# @Software: PyCharm
import numpy as np


def cmatrix(m_matrix, k_matrix, kesai):
    lam, vet = np.linalg.eig(np.dot(np.linalg.inv(m_matrix), k_matrix))
    c_fre = np.sort(np.sqrt(lam))  # 结构圆频率
    a = 2 * kesai * c_fre[0] * c_fre[1]/(c_fre[0] + c_fre[1])
    b = 2 * kesai / (c_fre[0] + c_fre[1])
    c_matrix = a * m_matrix + b * k_matrix
    return c_matrix


def kmatrix(k):
    k1 = list(k.copy())
    k1.pop(0)
    k2 = k1.copy()
    k2.append(0)
    k_matrix = np.diag(k) - np.diag(k1, 1) - np.diag(k1, -1) + np.diag(k2)
    return k_matrix
