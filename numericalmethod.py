# -*- coding = utf-8 -*-
# @Time:  10:04
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：SDOF.py
# @Software: PyCharm
import pandas as pd
import numpy as np

import math
import resultsprocessing as rp
import matrix
from resiliencemodel import HysteresisModelMdof as Hm
from resiliencemodel import HysteresisModelSdof as Hs


class SDOF:
    def __init__(self, ag, dt, period, damping_ratio, npts):
        self.ag = ag  # 地面加速度
        self.dt = dt  # 步长
        self.period = period  # 结构周期
        self.damping_ratio = damping_ratio  # 结构阻尼比
        self.npts = npts  # 地震动记录数量

    def newmark_beta(self):
        dis = [0] * self.npts
        vel = [0] * self.npts
        acc = [0] * self.npts
        gama = 0.5
        beta = 0.25
        p1_m = [0] * (self.npts + 1)
        # ------------------------系数------------------------#
        a1_m = 1 / beta / (self.dt ** 2) + gama / beta / self.dt * 4 * np.pi * self.damping_ratio / self.period
        a2_m = 1 / beta / self.dt + (gama / beta - 1) * 4 * np.pi * self.damping_ratio / self.period
        a3_m = (1 / 2 / beta - 1) + self.dt * (gama / 2 / beta - 1) * 4 * np.pi * self.damping_ratio / self.period
        k1_m = 4 * np.pi ** 2 / self.period ** 2 + a1_m  # k1与质量m的比值
        for i in range(1, self.npts):
            p1_m[i] = -self.ag[i] + a1_m * dis[i - 1] + a2_m * vel[i - 1] + a3_m * acc[i - 1]
            # 求解位移，输出单位m
            dis[i] = p1_m[i] / k1_m
            # 求解速度，输出单位m/s
            vel[i] = 3 / self.dt * (dis[i] - dis[i - 1]) - 2 * vel[i - 1] - self.dt / 2 * acc[i - 1]
            # 求解加速度，输出单位m/s^2
            acc[i] = 6 / self.dt ** 2 * (dis[i] - dis[i - 1]) - 6 / self.dt * vel[i - 1] - 2 * acc[i - 1]
        th = pd.DataFrame({
            'time(s)': np.arange(0, round(self.dt * self.npts, 5), self.dt),
            'Acc(g)': list(map(lambda x: x / 9.807, acc)),
            'Vel(cm/s)': list(map(lambda x: x * 100, vel)),
            'Dis(mm)': list(map(lambda x: x * 1000, dis))
        })
        th.to_excel('Newmark.xlsx', index=False)
        rp.draw(th, 'Newmark')

    def duhamel(self):
        """
        由于定积分法采用梯形法求解积分，存在一定误差，
        相比而言，递推法的求解精度一般高于定积分法，因
        此本文选用递推法
        """

        w = 2 * np.pi / self.period  # 求解结构频率
        w_d = w * math.sqrt(1 - self.damping_ratio ** 2)
        # --------------------------------递推系数a1，a2，a3k，a4k，b1，b2，b3k，b4k--------------------------------
        a1 = (math.cos(w_d * self.dt) + self.damping_ratio / math.sqrt(1 - self.damping_ratio ** 2) *
              math.sin(w_d * self.dt)) * np.exp(-self.damping_ratio * w * self.dt)
        a2 = (math.sin(w_d * self.dt) / w_d) * math.exp(-self.damping_ratio * w * self.dt)
        a3k = (2 * self.damping_ratio / w / self.dt + (-(1 + 2 * self.damping_ratio / w / self.dt) *
               math.cos(w_d * self.dt) + ((1 - 2 * self.damping_ratio ** 2) / w_d / self.dt - self.damping_ratio /
               math.sqrt(1 - self.damping_ratio ** 2)) * math.sin(w_d * self.dt)) *
               math.exp(-self.damping_ratio * w * self.dt))
        a4k = (1 - 2 * self.damping_ratio / w / self.dt + (2 * self.damping_ratio / w / self.dt *
               math.cos(w_d * self.dt) + (2 * self.damping_ratio ** 2 - 1) / w_d / self.dt * math.sin(w_d * self.dt)) *
               math.exp(-self.damping_ratio * w * self.dt))
        b1 = math.exp(-self.damping_ratio * w * self.dt) * (-w / math.sqrt(1 - self.damping_ratio ** 2) *
                                                            math.sin(w_d * self.dt))
        b2 = math.exp(-self.damping_ratio * w * self.dt) * (math.cos(w_d * self.dt) - math.sin(w_d * self.dt) *
                                                            self.damping_ratio / math.sqrt(1 - self.damping_ratio ** 2))
        b3k = (-1 / self.dt + math.exp(-self.damping_ratio * w * self.dt) * (1 / self.dt * math.cos(w_d * self.dt) + (
               w / math.sqrt(1 - self.damping_ratio ** 2) + self.damping_ratio / self.dt /
               math.sqrt(1 - self.damping_ratio ** 2)) * math.sin(w_d * self.dt)))
        b4k = (1 - math.exp(-self.damping_ratio * w * self.dt) * (
                    math.cos(w_d * self.dt) + self.damping_ratio / math.sqrt(1 - self.damping_ratio ** 2)
                    * math.sin(w_d * self.dt))) / self.dt
        # ------------------------------------------------------------------------------------------------------------ #
        t = np.array([[a1 + a3k, a2 + 2 * self.damping_ratio / w * a3k, a3k / w ** 2],
                      [b1 + b3k, b2 + 2 * self.damping_ratio / w * b3k, b3k / w ** 2],
                      [-w ** 2 * (a1 + a3k) - 2 * self.damping_ratio * w * (b1 + b3k),
                       -w ** 2 * (a2 + 2 * self.damping_ratio / w * a3k) - 2 *
                       self.damping_ratio * w * (b2 + 2 * self.damping_ratio / w * b3k), -a3k - 2 *
                       self.damping_ratio * b3k / w]])
        hm = np.array([[a4k / w ** 2], [b4k / w ** 2], [1 - 2 * self.damping_ratio * b4k / w - a4k]])
        dis = [0]
        vel = [0]
        acc = [0]
        x0 = np.zeros((3, 1))
        for i in range(1, len(self.ag)):
            response = np.dot(t, x0) + hm * self.ag[i]
            x0 = response
            dis.append(-float(response[0]))  # 提取位移，输出单位m
            vel.append(-float(response[1]))  # 提取速度，输出单位m/s
            acc.append(-float(response[2]))  # 提取加速度，输出单位m/s^2
        th = pd.DataFrame({
            'time(s)': np.arange(0, round(self.dt * self.npts, 5), self.dt),
            'Acc(g)': list(map(lambda x: x / 9.807, acc)),
            'Vel(cm/s)': list(map(lambda x: x * 100, vel)),
            'Dis(mm)': list(map(lambda x: x * 1000, dis))
        })
        th.to_excel('Duhamel.xlsx', index=False)
        rp.draw(th, 'Duhamel')

    def central_difference_method(self):
        dis = [0] * self.npts
        vel = [0] * self.npts
        acc = [0] * self.npts
        for i in range(2, self.npts):
            # 求解位移，输出单位m
            dis[i] = (-self.ag[i] + (2 / self.dt ** 2 - 4 * math.pi ** 2) * dis[i - 1] + (
                    0.2 * math.pi / 2 / self.dt - 1 / self.dt ** 2) * dis[i - 2]) / (
                                      1 / self.dt ** 2 + 0.2 * math.pi / 2 / self.dt)
            # 求解速度，输出单位m/s
            vel[i-1] = (dis[i] - dis[i-2]) / self.dt / 2
            # 求解加速度，输出单位m/s^2
            acc[i-1] = (dis[i] - 2 * dis[i-1] + dis[i-2]) / self.dt ** 2
        th = pd.DataFrame({
            'time(s)': np.arange(0, round(self.dt * self.npts, 5), self.dt),
            'Acc(g)': list(map(lambda x: x / 9.807, acc)),
            'Vel(cm/s)': list(map(lambda x: x * 100, vel)),
            'Dis(mm)': list(map(lambda x: x * 1000, dis))
        })
        th.to_excel('Central Difference Method.xlsx', index=False)
        rp.draw(th, 'Central Difference Method')


class NSDOF:
    def __init__(self, ag, dt, npts, period, m, strength, eta, model):
        self.ag = ag
        self.dt = dt
        self.npts = npts
        self.period = period
        self.m = m
        self.strength = strength
        self.eta = eta
        self.model = model

    def newmark_beta(self):
        # ----------------系统参数---------------- #
        damping_ratio = 0.05  # 结构阻尼比
        k = self.m * 4 * np.pi ** 2 / self.period ** 2  # 结构刚度计算
        c = 2 * np.sqrt(k * self.m) * damping_ratio  # 结构阻尼计算
        # ------------------NewMark方法的参数选择------------------ #
        gama = 1 / 2
        beta = 1 / 4  # 平均加速度法
        # beta = 1 / 6  # 线性加速度法
        # --------------------材料状态变量--------------------- #
        props = [k, self.strength, self.eta]
        state = [0, 0, 0, 0, 0, 0, 0]
        # ----------------定义初始条件---------------- #
        force = [0] * self.npts  # 定义初始结构抗力
        dis = [0] * self.npts  # 定义初始结构位移
        vel = [0] * self.npts  # 定义初始结构速度
        acc = [0] * self.npts  # 定义初始结构加速度
        p_ = [0] * (self.npts + 1)
        # ---------------------常数计算--------------------- #
        a1 = self.m / beta / self.dt ** 2 + gama / beta / self.dt * c
        a2 = self.m / beta / self.dt + (gama / beta - 1) * c
        a3 = (1 / 2 / beta - 1) * self.m + self.dt * (gama / 2 / beta - 1) * c
        # ------------------迭代开始-----------------------#
        for i in range(0, self.npts - 1):
            p_[i + 1] = -self.ag[i + 1] * self.m + a1 * dis[i] + a2 * vel[i] + a3 * acc[i]
            force_j = force[i]
            dis_j = dis[i]
            if i == 0:
                kt = k
            else:
                hm = Hs(props, force_j, dis_j, 0, state, abs(dis[i]) - abs(dis[i - 1]))
                if self.model == 'bilinear model':
                    force_j, kt = hm.bilinear_model()
                elif self.model == 'clough model':
                    force_j, kt, state = hm.clough_model()
                elif self.model == 'slip model':
                    force_j, kt, state = hm.slip_model()
                elif self.model == 'origin oriented model':
                    force_j, kt = hm.origin_oriented_model()
                else:
                    force_j, kt = hm.bilinear_elastic_model()
            for j in range(1, 200):
                r_j = p_[i + 1] - force_j - a1 * dis_j
                if abs(r_j) > 0.0005:  # 收敛检查
                    ke = kt + a1
                    dt_dis = r_j / ke
                    dis_j = dis_j + dt_dis
                    hm = Hs(props, force_j, dis_j, dt_dis, state, dt_dis)
                    if self.model == 'bilinear model':
                        force_j, kt = hm.bilinear_model()
                    elif self.model == 'clough model':
                        force_j, kt, state = hm.clough_model()
                    elif self.model == 'slip model':
                        force_j, kt, state = hm.slip_model()
                    elif self.model == 'origin oriented model':
                        force_j, kt = hm.origin_oriented_model()
                    else:
                        force_j, kt = hm.bilinear_elastic_model()
                else:
                    dis[i + 1] = dis_j  # 最终求得的位移
                    break
            force[i + 1] = force_j  # 最终求得的结构抗力
            vel[i + 1] = (gama / beta / self.dt * (dis[i + 1] - dis[i]) + (1 - gama / beta) * vel[i] + self.dt *
                          (1 - gama / 2 / beta) * acc[i])  # 最终求得的结构速度
            acc.append((dis[i + 1] - dis[i]) / beta / self.dt ** 2 - vel[i] / beta / self.dt - (1 / 2 / beta - 1) *
                       acc[i])  # 最终求得的结构加速度
        # --------------------------------时程曲线的保存-------------------------------- #
        th = pd.DataFrame({
            'time(s)': np.arange(0, round(self.dt * self.npts, 5), self.dt),
            'Acc(g)': list(map(lambda x: x / 9.807, acc)),
            'Vel(cm/s)': list(map(lambda x: x * 100, vel)),
            'Dis(mm)': list(map(lambda x: x * 1000, dis)),
            'force(N)': force
            # 'force(N)': list(map(lambda x: x * 1000, force))
        })
        rp.draw(th, self.model)  # 调用绘图函数
        th.to_excel('TH {}.xlsx'.format(self.model), index=False)  # 将数据保存为Excel文件


class MDOF:
    def __init__(self, m_matrix, k_matrix, c_matrix, step, npts, ag, nds):
        self.M = m_matrix
        self.K = k_matrix
        self.C = c_matrix
        self.step = step
        self.npts = npts
        self.ag = ag
        self.nds = nds

    def central_difference_method(self):
        lam, vet = np.linalg.eig(np.dot(np.linalg.inv(self.M), self.K))
        idx = np.argsort(lam)
        fai = vet[:, idx][:, 0:2]
        M = np.dot(np.dot(fai.T, self.M), fai)
        K = np.dot(np.dot(fai.T, self.K), fai)
        C = np.dot(np.dot(fai.T, self.C), fai)
        q1 = np.zeros([2, self.npts+1])
        q2 = np.zeros([2, self.npts+1])
        q3 = np.zeros([2, self.npts+1])
        dis = np.zeros([self.nds, self.npts])
        vel = np.zeros([self.nds, self.npts])
        acc = np.zeros([self.nds, self.npts])
        K_ = 1 / self.step ** 2 * M + 1 / 2 / self.step * C
        a = 1 / self.step ** 2 * M - 1 / 2 / self.step * C
        b = K - 2 / self.step ** 2 * M
        for i in range(1, self.npts-1):
            P = (-np.dot(fai.T, np.dot(self.M * self.ag[i+1], np.ones(self.nds))) -
                 np.dot(a, q1[:, i-1]) - np.dot(b, q1[:, i]))
            q1[:, i+1] = np.dot(np.linalg.inv(K_), P)
            q2[:, i] = 1 / 2 / self.step * (q1[:, i+1] - q1[:, i-1])
            q3[:, i] = 1 / self.step ** 2 * (q1[:, i+1] - 2 * q1[:, i] + q1[:, i-1])
            dis[:, i+1] = np.dot(fai, q1[:, i+1])
            vel[:, i] = np.dot(fai, q2[:, i])
            acc[:, i] = np.dot(fai, q3[:, i])
        dm = rp.MdofResult(self.nds, np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis)
        dm.draw_history()
        rp.Write(np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, self.nds).mdof()

    def newmark_beta(self):
        lam, vet = np.linalg.eig(np.dot(np.linalg.inv(self.M), self.K))
        idx = np.argsort(lam)
        fai = vet[:, idx][:, 0:2]
        M = np.dot(np.dot(fai.T, self.M), fai)
        K = np.dot(np.dot(fai.T, self.K), fai)
        C = np.dot(np.dot(fai.T, self.C), fai)
        q1 = np.zeros([2, self.npts])
        q2 = np.zeros([2, self.npts])
        q3 = np.zeros([2, self.npts])
        dis = np.zeros([self.nds, self.npts])
        vel = np.zeros([self.nds, self.npts])
        acc = np.zeros([self.nds, self.npts])
        gama = 1/2
        beta = 1/4
        a1 = 1 / beta / self.step ** 2 * M + gama / beta / self.step * C
        a2 = 1 / beta / self.step * M + (gama / beta - 1) * C
        a3 = (1 / 2 / beta - 1) * M + (gama / 2 / beta - 1) * self.step * C
        K_ = K + a1
        for i in range(0, self.npts-1):
            P = (-np.dot(fai.T, np.dot(self.M * self.ag[i+1], np.ones(self.nds))) + np.dot(a1, q1[:, i]) +
                 np.dot(a2, q2[:, i]) + np.dot(a3, q3[:, i]))
            q1[:, i+1] = np.dot(np.linalg.inv(K_), P)
            q2[:, i+1] = (gama / beta / self.step * (q1[:, i+1] - q1[:, i]) + (1 - gama / beta) * q2[:, i] +
                          self.step * (1 - gama / 2 / beta) * q3[:, i])
            q3[:, i+1] = (1 / beta / self.step ** 2 * (q1[:, i+1] - q1[:, i]) - 1 / beta / self.step * q2[:, i]
                          - (1 / 2 / beta - 1) * q3[:, i])
            dis[:, i+1] = np.dot(fai, q1[:, i+1])
            vel[:, i+1] = np.dot(fai, q2[:, i+1])
            acc[:, i+1] = np.dot(fai, q3[:, i+1])
        dm = rp.MdofResult(self.nds, np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis)
        dm.draw_history()
        rp.Write(np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, self.nds).mdof()


class NMDOF:
    def __init__(self, m_matrix, k_matrix, c_matrix, step, npts, ag, model, nds, props):
        self.m_matrix = m_matrix
        self.k_matrix = k_matrix
        self.c_matrix = c_matrix
        self.step = step
        self.npts = npts
        self.ag = ag
        self.model = model
        self.nds = nds
        self.props = props

    def newmark_beta(self):
        force = np.zeros([self.nds, self.npts])
        dis = np.zeros([self.nds, self.npts])
        vel = np.zeros([self.nds, self.npts])
        acc = np.zeros([self.nds, self.npts])
        states = np.zeros([self.nds, 7])
        # -------------------------------------------NewMark-Raphson------------------------------------------- #
        # ---------------------常数计算--------------------- #
        a1 = 4 * self.m_matrix / self.step ** 2 + 2 * self.c_matrix / self.step
        a2 = self.m_matrix * 4 / self.step + self.c_matrix
        for i in range(0, self.npts - 1):
            p_ = (-np.dot(self.m_matrix * self.ag[i + 1], np.ones(self.nds)) + np.dot(a1, dis[:, i]) +
                  np.dot(a2, vel[:, i]))
            force_j = force[:, i]
            dis_j = dis[:, i]
            if i == 0:
                kt_j = self.k_matrix
            else:
                hm = Hm(self.props, force_j, dis_j, np.zeros(self.nds), states, self.nds,
                        abs(dis[:, i]) - abs(dis[:, i - 1]))
                if self.model == 'UHYST01':
                    force_j, kt_j = hm.bilinear_model()
                    kt_j = matrix.kmatrix(kt_j)
                else:
                    force_j, kt_j, states = hm.clough_model()
                    kt_j = matrix.kmatrix(kt_j)
            for j in range(0, 200):
                R_j = p_ - force_j - np.dot(a1, dis_j)
                if np.linalg.norm(R_j, 2) > 0.0001:
                    ke = kt_j + a1
                    dt_dis = np.dot(np.linalg.inv(ke), R_j)
                    hm = Hm(self.props, force_j, dis_j, dt_dis, states, self.nds, dt_dis)
                    dis_j = dis_j + dt_dis
                    if self.model == 'UHYST01':
                        force_j, kt_j = hm.bilinear_model()
                        kt_j = matrix.kmatrix(kt_j)
                    elif self.model == 'UHYST02':
                        force_j, kt_j, states = hm.clough_model()
                        kt_j = matrix.kmatrix(kt_j)
                    else:
                        force_j, kt_j = hm.origin_oriented_model()
                        kt_j = matrix.kmatrix(kt_j)
                else:
                    dis[:, i + 1] = dis_j  # 最终求得的位移
                    break
            force[:, i + 1] = force_j  # 最终求得的结构抗力
            vel[:, i + 1] = 2 / self.step * (dis[:, i + 1] - dis[:, i]) - vel[:, i]  # 最终求得的结构速度
            acc[:, i + 1] = (4 / self.step ** 2 * (dis[:, i + 1] - dis[:, i]) -
                             4 / self.step * vel[:, i] - acc[:, i])  # 最终求得的结构加速度
        dm = rp.MdofResult(self.nds, np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, force)
        dm.draw_hysteretic_curve()
        dm.draw_history()
        rp.Write(np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, self.nds, force).nmdof()

    def central_difference_method(self):
        force = np.zeros([self.nds, self.npts+1])
        dis = np.zeros([self.nds, self.npts+1])
        vel = np.zeros([self.nds, self.npts+1])
        acc = np.zeros([self.nds, self.npts+1])
        states = np.zeros([self.nds, 7])
        k_ = 1 / self.step ** 2 * self.m_matrix + 1 / 2 / self.step * self.c_matrix
        a = 1 / self.step ** 2 * self.m_matrix - 1 / 2 / self.step * self.c_matrix
        b = - 2 / self.step ** 2 * self.m_matrix
        for i in range(1, self.npts-1):
            p_ = (-np.dot(self.m_matrix * self.ag[i + 1], np.ones(self.nds)) - np.dot(a, dis[:, i-1]) -
                  np.dot(b, dis[:, i]) - force[:, i])
            dis[:, i+1] = np.dot(np.linalg.inv(k_), p_)
            hm = Hm(self.props, force[:, i], dis[:, i], np.zeros(self.nds), states, self.nds,
                    abs(dis[:, i]) - abs(dis[:, i - 1]))
            if self.model == 'UHYST01':
                force[:, i+1] = hm.bilinear_model()[0]
            else:
                force[:, i+1] = hm.clough_model()[0]
            vel[:, i] = 1 / 2 / self.step * (dis[:, i+1] - dis[:, i-1])
            acc[:, i] = 1 / self.step ** 2 * (dis[:, i+1] - 2 * dis[:, i] + dis[:, i-1])
        dm = rp.MdofResult(self.nds, np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, force)
        dm.draw_hysteretic_curve()
        dm.draw_history()
        rp.Write(np.arange(0, round(self.step * self.npts, 5), self.step), acc, vel, dis, self.nds, force).nmdof()
