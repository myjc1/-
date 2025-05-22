# -*- coding = utf-8 -*-
# @Time:  16:12
# @Author:Wang Maocen
# @E-mail:wangmaocen_1999@163.com
# @File：processing.py
# @Software: PyCharm
import matplotlib.pyplot as plt
import numpy as np
import xlsxwriter
# --------------------------定义全局图像字体----------------------------------- #
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 步骤一(替换sans-serif字体)
plt.rcParams['axes.unicode_minus'] = False  # 步骤二(解决坐标轴负数的负号显示问题)
plt.rcParams['lines.linewidth'] = 2.5  # 定义线宽
plt.rcParams['xtick.labelsize'] = 30  # 定义x轴坐标字体大小
plt.rcParams['ytick.labelsize'] = 30  # 定义y轴坐标字体大小
plt.rcParams['figure.figsize'] = (15, 10)  # 定义图像大小
plt.rcParams['axes.labelsize'] = 35  # 定义坐标轴标题字体大小
plt.rcParams['legend.fontsize'] = 25  # 定义标签字体大小


class MdofResult:
    def __init__(self, nds, time, acc, vel, dis, f=np.zeros([1, 1])):
        self.time = time
        self.acc = acc
        self.vel = vel
        self.dis = dis
        self.f = f
        self.nds = nds

    def draw_hysteretic_curve(self):
        F = np.zeros(self.f.shape)
        story_drift = np.zeros(self.f.shape)
        for story in range(0, self.nds):
            for j in range(story, self.nds):
                F[story, :] = F[story, :] + self.f[j, :]
            if story == 0:
                story_drift[story, :] = self.dis[story, :]
                plt.figure(story)
                plt.plot(story_drift[story, :] * 1000, F[story, :] / 1000)
                plt.xlabel('Dis(mm)')
                plt.ylabel('F(kN)')
                plt.title('hysteretic curve of {} story'.format(story), fontsize=40)
                plt.show()
            else:
                story_drift[story, :] = self.dis[story, :] - self.dis[story - 1, :]
                plt.figure(story)
                plt.plot(story_drift[story, :] * 1000, F[story, :] / 1000)
                plt.xlabel('Dis(mm)')
                plt.ylabel('F(kN)')
                plt.title('hysteretic curve of {} story'.format(story + 1), fontsize=40)
                plt.show()

    def draw_history(self):
        for story in range(0, self.nds):
            plt.figure(0)
            plt.plot(self.time, self.acc[story, :] / 9.807)
            plt.xlabel('time(s)')
            plt.ylabel('acc(g)')
            plt.title('Acc of {} story'.format(story + 1), fontsize=40)
            plt.figure(1)
            plt.plot(self.time, self.vel[story, :] * 100)
            plt.xlabel('time(s)')
            plt.ylabel('Vel(cm/s)')
            plt.title('Vel of {} story'.format(story + 1), fontsize=40)
            plt.figure(2)
            plt.plot(self.time, self.dis[story, :] * 1000)
            plt.xlabel('time(s)')
            plt.ylabel('Dis(mm)')
            plt.title('Dis of {} story'.format(story + 1), fontsize=40)
            plt.show()


def draw(history_data, method_name):
    time = history_data.iloc[:, 0]
    ylabel = ['Acc(g)', 'Vel(cm/s)', 'Dis(mm)']
    for i in range(0, 3):
        plt.figure(i)
        plt.plot(time, history_data.iloc[:, i+1])
        plt.xlim(min(time), max(time))
        plt.xlabel('time(s)')
        plt.ylabel(ylabel[i])
        plt.title(method_name, fontsize=40)
        plt.show()


class Write:
    def __init__(self, t, a, v, d, ndf, f=0, method_name=None):
        self.t = t
        self.f = f
        self.a = a
        self.v = v
        self.d = d
        self.ndf = ndf
        self.mn = method_name

    def sdof(self):
        pass

    def mdof(self):
        th = xlsxwriter.Workbook('th.xlsx')
        for sheet in ['acc', 'vel', 'dis']:
            th_sheet = th.add_worksheet(sheet)
            th_sheet.write(0, 0, 'time(s)')
            th_sheet.write_column(1, 0, self.t)
            if sheet == 'acc':
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.a[story, :])
            elif sheet == 'vel':
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.v[story, :])
            else:
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.d[story, :])
        th.close()

    def nmdof(self):
        th = xlsxwriter.Workbook('th.xlsx')
        for sheet in ['acc', 'vel', 'dis', 'force']:
            th_sheet = th.add_worksheet(sheet)
            th_sheet.write(0, 0, 'time(s)')
            th_sheet.write_column(1, 0, self.t)
            if sheet == 'acc':
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.a[story, :])
            elif sheet == 'vel':
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.v[story, :])
            elif sheet == 'dis':
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.d[story, :])
            else:
                for story in range(0, self.ndf):
                    th_sheet.write(0, story + 1, 'floor{}'.format(story + 1))
                    th_sheet.write_column(1, story + 1, self.f[story, :])
        th.close()
