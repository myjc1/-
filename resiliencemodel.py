# -*- coding = utf-8 -*-
# @Time:  10:22
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：resiliencemodel.py
# @Software: PyCharm
import numpy as np


class HysteresisModelMdof:
    """
    该类中包含了五种骨架曲线模型，适用于多自由度模型
    props：初始刚度，屈服强度，应力变化率
    s：结构抗力
    e：结构位移
    de：第j次迭代的位移变化量
    dv：结构最终位移值dis(i+1)和dis(i)的绝对值的差值
    """
    def __init__(self, prop, s, e, de, state, nds, dv):
        self.props = prop  # n*3
        self.s = s  # n*1
        self.e = e  # n*1
        self.de = de  # n*n
        self.state = state
        self.nodes = nds
        self.dv = dv

    def bilinear_model(self):
        """
        force为结构抗力
        kt为结构现阶段下的刚度
        evs为进入屈服状态后不同位移下的作用力
        eve为屈服后的刚度
        """
        F = np.zeros(self.nodes)
        force = np.zeros(self.nodes)
        kt = np.zeros(self.nodes)
        drift = [self.e[0]]
        d_drift = [self.de[0]]
        for i in range(1, self.nodes):
            drift.append(self.e[i]-self.e[i-1])
            d_drift.append(self.de[i] - self.de[i - 1])
        for i in range(0, self.nodes):
            for j in range(i, self.nodes):
                F[i] = F[i] + self.s[j]
        for i in range(0, self.nodes):
            force[i] = F[i] + self.props[i, 0] * d_drift[i]
            kt[i] = self.props[i, 0]
            eve = self.props[i, 0] * self.props[i, 2]
            if d_drift[i] >= 0:
                evs = self.props[i, 1] + self.props[i, 2] * self.props[i, 0] * (drift[i] + d_drift[i] - self.props[i, 1]
                                                                                / self.props[i, 0])
                if force[i] >= evs:
                    force[i] = evs
                    kt[i] = eve
            else:
                evs = -self.props[i, 1] + self.props[i, 2] * self.props[i, 0] * (drift[i] + d_drift[i] +
                                                                                 self.props[i, 1] / self.props[i, 0])
                if force[i] <= evs:
                    force[i] = evs
                    kt[i] = eve
        for i in range(0, self.nodes):
            if i == (self.nodes - 1):
                F[i] = force[i]
            else:
                F[i] = force[i] - force[i+1]
        return F, kt

    def clough_model(self):
        F = np.zeros(self.nodes)
        force = np.zeros(self.nodes)
        kt = np.zeros(self.nodes)
        drift = [self.e[0]]
        d_drift = [self.de[0]]
        state = np.zeros([self.nodes, 7])
        for i in range(1, self.nodes):
            drift.append(self.e[i] - self.e[i - 1])
            d_drift.append(self.de[i] - self.de[i - 1])
        for i in range(0, self.nodes):
            for j in range(i, self.nodes):
                F[i] = F[i] + self.s[j]
        for story in range(0, self.nodes):
            e_max = self.state[story, 0]
            e_min = self.state[story, 1]
            ert = self.state[story, 2]
            srt = self.state[story, 3]
            erc = self.state[story, 4]
            src = self.state[story, 5]
            kon = self.state[story, 6]
            if kon == 0:
                e_max = self.props[story, 1] / self.props[story, 0]
                e_min = -self.props[story, 1] / self.props[story, 0]
                if d_drift[story] >= 0:
                    kon = 1
                else:
                    kon = 2
            elif kon == 1 and d_drift[story] < 0:
                kon = 2
                if F[story] > 0:
                    erc = drift[story]
                    src = F[story]
                if drift[story] > e_max:
                    e_max = drift[story]
            elif kon == 2 and d_drift[story] >= 0:
                kon = 1
                if F[story] < 0:
                    ert = drift[story]
                    srt = F[story]
                if drift[story] < e_min:
                    e_min = drift[story]
            force[story] = F[story] + self.props[story, 0] * d_drift[story]
            kt[story] = self.props[story, 0]
            eve = self.props[story, 0] * self.props[story, 2]
            sres = 0
            if d_drift[story] >= 0:
                evs = self.props[story, 1] + (drift[story] + d_drift[story] - self.props[story, 1] /
                                              self.props[story, 0]) * self.props[story, 2] * self.props[story, 0]
                if force[story] >= evs:
                    force[story] = evs
                    kt[story] = eve
                s_max = max(self.props[story, 1], self.props[story, 1] +
                            (e_max - self.props[story, 1] / self.props[story, 0]) *
                            self.props[story, 2] * self.props[story, 0])
                e_res = ert - (srt - sres) / self.props[story, 0]
                if e_res <= e_max - s_max / self.props[story, 0]:
                    s_rel = (drift[story] + d_drift[story] - e_res) / (e_max - e_res) * (s_max - sres) + sres
                    if force[story] > s_rel:
                        force[story] = s_rel
                        kt[story] = (s_max - sres) / (e_max - e_res)
            else:
                evs = -self.props[story, 1] + (drift[story] + d_drift[story] + self.props[story, 1] /
                                               self.props[story, 0]) * self.props[story, 2] * self.props[story, 0]
                if force[story] <= evs:
                    force[story] = evs
                    kt[story] = eve
                s_min = min(-self.props[story, 1], -self.props[story, 1] +
                            (e_min + self.props[story, 1] / self.props[story, 0]) *
                            self.props[story, 2] * self.props[story, 0])
                e_res = erc - (src - sres) / self.props[story, 0]
                if e_res >= e_min - s_min / self.props[story, 0]:
                    s_rel = (drift[story] + d_drift[story] - e_res) / (e_min - e_res) * (s_min - sres) + sres
                    if force[story] < s_rel:
                        force[story] = s_rel
                        kt[story] = (s_min - sres) / (e_min - e_res)
            state[story, :] = [e_max, e_min, ert, srt, erc, src, kon]
        for i in range(0, self.nodes):
            if i == (self.nodes - 1):
                F[i] = force[i]
            else:
                F[i] = force[i] - force[i+1]
        return F, kt, state

    def origin_oriented_model(self):
        F = np.zeros(self.nodes)
        force = np.zeros(self.nodes)
        kt = np.zeros(self.nodes)
        drift = [self.e[0]]
        d_drift = [self.de[0]]
        d_dis = [self.dv[0]]
        for i in range(1, self.nodes):
            drift.append(self.e[i]-self.e[i-1])
            d_drift.append(self.de[i] - self.de[i - 1])
            d_dis.append(self.dv[i] - self.dv[i - 1])
        for i in range(0, self.nodes):
            for j in range(i, self.nodes):
                F[i] = F[i] + self.s[j]
        for i in range(0, self.nodes):
            if drift[i] * d_drift[i] > 0:
                force[i] = self.props[i, 0] * (drift[i] + d_drift[i])
                kt[i] = self.props[i, 0]
                evs = (self.props[i, 1] + (abs(drift[i] + d_drift[i]) - self.props[i, 1] /
                                           self.props[i, 0]) * self.props[i, 2] * self.props[i, 0])
                eve = self.props[i, 2] * self.props[i, 0]
                if abs(force[i]) >= evs:
                    if d_drift[i] == 0:
                        force[i] = evs
                        kt[i] = eve
                    else:
                        force[i] = np.sign(d_drift[i]) * evs
                        kt[i] = np.sign(d_drift[i]) * eve
            elif d_drift[i] == 0:
                if d_dis[i] < 0:
                    kt[i] = F[i] / drift[i]
                    force[i] = kt[i] * drift[i]
                else:
                    force[i] = self.props[i, 0] * drift[i]
                    kt[i] = self.props[i, 0]
                    evs = (self.props[i, 1] + (abs(drift[i] + d_drift[i]) - self.props[i, 1] /
                                               self.props[i, 0]) * self.props[i, 2] * self.props[i, 0])
                    eve = self.props[i, 2] * self.props[i, 0]
                    if abs(force[i]) >= evs:
                        force[i] = evs * np.sign(drift[i])
                        kt[i] = eve
            elif F[i] * d_drift[i] < 0 and drift[i] != 0:  # 卸载过程
                kt[i] = F[i] / drift[i]
                force[i] = F[i] + kt[i] * d_drift[i]
            else:
                kt[i] = self.props[i, 0]
                force[i] = 0
        for i in range(0, self.nodes):
            if i == (self.nodes - 1):
                F[i] = force[i]
            else:
                F[i] = force[i] - force[i+1]
        return F, kt

    def bilinear_elastic_model(self):
        pass


class HysteresisModelSdof:
    """
    该类中包含了五种骨架曲线模型，应用于单自由度
    props：初始刚度，屈服强度，应力变化率
    s：结构抗力
    e：结构位移
    de：第j次迭代的位移变化量
    dv：结构最终位移值dis(i+1)和dis(i)的绝对值的差值
    """
    def __init__(self, props, s, e, de, state, dv):
        self.props = props
        self.s = s
        self.e = e
        self.de = de
        self.state = state
        self.dv = dv

    def bilinear_model(self):
        """
        双线性模型
        force为结构抗力
        kt为结构现阶段下的刚度
        evs为进入屈服状态后不同位移下的作用力
        eve为屈服后的刚度
        """
        force = self.s + self.props[0] * self.de
        kt = self.props[0]
        if self.de >= 0:
            evs = self.props[1] + self.props[2] * self.props[0] * (self.e + self.de - self.props[1] / self.props[0])
            eve = self.props[0] * self.props[2]
            if force >= evs:
                force = evs
                kt = eve
        else:
            evs = -self.props[1] + self.props[2] * self.props[0] * (self.e + self.de + self.props[1] / self.props[0])
            eve = self.props[0] * self.props[2]
            if force <= evs:
                force = evs
                kt = eve
        return force, kt

    def clough_model(self):
        """
        克拉夫模型
        :return: 新的恢复力以及刚度
        """
        e_max = self.state[0]
        e_min = self.state[1]
        ert = self.state[2]
        srt = self.state[3]
        erc = self.state[4]
        src = self.state[5]
        kon = self.state[6]
        if kon == 0:
            e_max = self.props[1] / self.props[0]
            e_min = -self.props[1] / self.props[0]
            if self.de >= 0:
                kon = 1
            else:
                kon = 2
        elif kon == 1 and self.de < 0:
            kon = 2
            if self.s > 0:
                erc = self.e
                src = self.s
            if self.e > e_max:
                e_max = self.e
        elif kon == 2 and self.de >= 0:
            kon = 1
            if self.s < 0:
                ert = self.e
                srt = self.s
            if self.e < e_min:
                e_min = self.e
        force = self.s + self.props[0] * self.de
        kt = self.props[0]
        eve = self.props[2] * self.props[0]
        sres = 0
        if self.de >= 0:
            evs = self.props[1] + (self.e + self.de - self.props[1] / self.props[0]) * self.props[2] * self.props[0]
            if force >= evs:
                force = evs
                kt = eve
            s_max = max(self.props[1], self.props[1] + (e_max - self.props[1] / self.props[0]) * self.props[2] *
                        self.props[0])
            e_res = ert - (srt - sres) / self.props[0]
            if e_res <= e_max - s_max / self.props[0]:
                s_rel = (self.e + self.de - e_res) / (e_max - e_res) * (s_max - sres) + sres
                if force > s_rel:
                    force = s_rel
                    kt = (s_max - sres) / (e_max - e_res)
        else:
            evs = -self.props[1] + (self.e + self.de + self.props[1] / self.props[0]) * self.props[2] * self.props[0]
            if force <= evs:
                force = evs
                kt = eve
            s_min = min(-self.props[1], -self.props[1] + (e_min + self.props[1] / self.props[0]) * self.props[2] *
                        self.props[0])
            e_res = erc - (src - sres) / self.props[0]
            if e_res >= e_min - s_min / self.props[0]:
                s_rel = (self.e + self.de - e_res) / (e_min - e_res) * (s_min - sres) + sres
                if force < s_rel:
                    force = s_rel
                    kt = (s_min - sres) / (e_min - e_res)
        state = [e_max, e_min, ert, srt, erc, src, kon]
        return force, kt, state

    def slip_model(self):
        """
        滑移型模型
        :return:
        """
        e_max = 0.01 * self.props[1] / self.props[0]
        e_min = -0.01 * self.props[1] / self.props[0]
        ert = self.state[0]
        srt = self.state[1]
        erc = self.state[2]
        src = self.state[3]
        kon = round(self.state[4])
        if kon == 0:
            if self.de >= 0:
                kon = 1
            else:
                kon = 2
        elif kon == 1 and self.de < 0:
            kon = 2
            if self.s > 0:
                erc = self.e
                src = self.s
        elif kon == 2 and self.de >= 0:
            kon = 1
            if self.s < 0:
                ert = self.e
                srt = self.s
        force = self.s + self.props[0] * self.de
        kt = self.props[0]
        sres = 0
        if self.de >= 0:
            evs = self.props[1] + (self.e + self.de - self.props[1] / self.props[0]) * self.props[2] * self.props[0]
            eve = self.props[2] * self.props[0]
            if force >= evs:
                force = evs
                kt = eve
            s_max = 0.01 * self.props[1]
            e_res = ert - (srt - sres) / self.props[0]
            if e_res <= e_max - s_max/self.props[0]:
                s_rel = (self.e + self.de - e_res) / (e_max - e_res) * (s_max - sres) + sres
                if force > s_rel and (self.e + self.de) < e_max:
                    force = s_rel
                    kt = (s_max - sres) / (e_max - e_res)
        else:
            evs = -self.props[1] + (self.e + self.de + self.props[1] / self.props[0]) * self.props[2] * self.props[0]
            eve = self.props[2] * self.props[0]
            if force <= evs:
                force = evs
                kt = eve
            s_min = -0.01 * self.props[1]
            e_res = erc - (src - sres) / self.props[0]
            if e_res >= e_min - s_min / self.props[0]:
                s_rel = (self.e + self.de - e_res) / (e_min - e_res) * (s_min - sres) + sres
                if force < s_rel and (self.e + self.de) > e_max:
                    force = s_rel
                    kt = (s_min - sres) / (e_min - e_res)
        return force, kt, [ert, srt, erc, src, kon]

    def origin_oriented_model(self):
        """
        原点指向性模型
        :return:
        """
        if self.s * self.de > 0:  # 加载过程
            force = self.props[0] * (self.e + self.de)
            kt = self.props[0]
            evs = (self.props[1] + (abs(self.e + self.de) - self.props[1] / self.props[0]) * self.props[2] *
                   self.props[0])
            eve = self.props[2] * self.props[0]
            if abs(force) >= evs:
                if self.de == 0:
                    force = evs
                    kt = eve
                else:
                    force = np.sign(self.de) * evs
                    kt = np.sign(self.de) * eve
        elif self.de == 0:
            if self.dv < 0:
                kt = self.s / self.e
                force = kt * self.e
            else:
                force = self.props[0] * self.e
                kt = self.props[0]
                evs = (self.props[1] + (abs(self.e + self.de) - self.props[1] / self.props[0]) * self.props[2] *
                       self.props[0])
                eve = self.props[2] * self.props[0]
                if abs(force) >= evs:
                    force = evs * np.sign(self.e)
                    kt = eve
        elif self.s * self.de < 0 and self.e != 0:  # 卸载过程
            kt = self.s / self.e
            force = self.s + kt * self.de
        else:
            kt = self.props[0]
            force = 0
        return force, kt

    def bilinear_elastic_model(self):
        if abs(self.e + self.de) <= self.props[1] / self.props[0]:  # 线性阶段
            force = self.props[0] * (self.e + self.de)
            kt = self.props[0]
        else:  # 屈服阶段
            force = np.sign(self.e + self.de) * (self.props[1] + (abs(self.e + self.de) - self.props[1] / self.props[0])
                                                 * self.props[0] * self.props[2])
            kt = self.props[2] * self.props[0]
        return force, kt
