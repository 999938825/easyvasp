import numpy as np
import math
from scipy import constants as con
import matplotlib.pyplot as plt
import os

class Toolkit(object):
    def __init__(self):
        pass

    @staticmethod
    def kcal_to_eV(value):
        '''
        :param value: 电子伏特
        :return: 千卡每摩尔
        '''
        return value/23.060386

    @staticmethod
    def eV_to_kcal(value):
        '''
        :param value: 千卡每摩尔
        :return: 电子伏特
        '''
        return value*23.060386

    @staticmethod
    def eV_to_kj(value):
        '''
        :param value: 电子伏特
        :return: 千焦每摩尔
        '''
        return value*96.484587

    @staticmethod
    def kj_to_eV(value):
        '''
        :param value: 千焦每摩尔
        :return: 电子伏特
        '''
        return value/96.484587

    @staticmethod
    def kcal_to_kj(value):
        '''
        :param value: 千卡每摩尔
        :return: 千焦每摩尔
        '''
        return value*4.183997

    @staticmethod
    def kj_to_kcal(value):
        '''
        :param value: 千焦每摩尔
        :return: 千卡每摩尔
        '''
        return value/4.183997

    @staticmethod
    def eV_to_v(eV):
        return eV*8049

    @staticmethod
    def v_to_eV(v):
        return v/8049

    @staticmethod
    def calculate_reaction_rate_constant_by_free_energy_barrier(G, T,sigma=1,te=1):
        '''
        计算单分子反应的化学反应常数，参考 http://bbs.keinsci.com/thread-14313-1-6.html
        :param 活化能G,(单位电子伏特eV), 反应温度T(单位开尔文K):
        :param sigma 反应路径简并度，te隧穿效应矫正系数
        :return 反应常数k:
        '''
        k = sigma * te * 2.1 * 10 ** 10 * T * math.exp(-1000 * 23.060386 * G / (T * 1.9859))
        return k

    @staticmethod
    def calculate_half_life_period_by_k(k,p):
        '''
        :param k: 化学反应速率常数
        :param p: 化学反应进度,如完成百分之95则取0.95
        :return: 所需时间，单位(小时)
        '''
        return -math.log(1-p)/k

    @staticmethod
    def calculate_transmission_coefficient_1(ifreq,T):
        '''
        :param ifreq: 体系的虚频 单位(cm-1)
        :param T: 温度(开尔文)
        :return: 透射系数
        使用Wigner方法计算透射系数，低温时精度差
        本部分内容参考Tian Lu, TSTcalculator, http://sobereva.com/310 (accessed 月 日, 年)
        '''
        return 1+(con.h/con.k/T*ifreq*30000000000)**2/24

    @staticmethod
    def calculate_transmission_coefficient_2(ifreq,T):
        '''
        :param ifreq: 体系的虚频 单位(cm-1)
        :param T: 温度(开尔文)
        :return: 透射系数
        使用近似的Skodje-Truhlar方法计算透射系数，低温时精度差
        本部分内容参考Tian Lu, TSTcalculator, http://sobereva.com/310 (accessed 月 日, 年)
        '''
        alpha=2*math.pi/con.h/(ifreq*30000000000)
        beta=1/(con.k*T)
        if alpha>=beta:
            return beta*math.pi/alpha/math.sin(beta*math.pi/alpha)
        else:
            print('The approximate Skodje-Truhlar method is not applicable, and the Wigner method is used')
            return Toolkit.calculate_transmission_coefficient_1(ifreq,T)

    @staticmethod
    def calculate_transmission_coefficient_3(ifreq, T,dG,dU):
        '''
        :param ifreq: 体系的虚频 单位(cm-1)
        :param T: 温度(开尔文)
        :param dG: 正向势垒(开尔文)
        :param dU: 放热反应能量为0，吸热反应为吸热数值，(开尔文)
        :return: 透射系数
        本部分内容参考Tian Lu, TSTcalculator, http://sobereva.com/310 (accessed 月 日, 年)
        '''
        dU=max(0,dU)
        alpha=2*math.pi/con.h/(ifreq*30000000000)
        beta=1/(con.k*T)
        if alpha<=beta:
            return beta/(beta-alpha)*(math.exp((beta-alpha)*(dG-dU)*1000/con.N_A)-1)
        else:
            return Toolkit.calculate_transmission_coefficient_2(ifreq,T)-beta/(alpha-beta)*(math.exp((beta-alpha)*(dG-dU)*1000/con.N_A))

    @staticmethod
    def find(file='CONTCAR',path=None):
        if path is not None:
            rootpath=path
        else:
            rootpath=os.getcwd()
        list=[]
        for top, dirs, files in os.walk(rootpath):
            topsz = top.split('\\')
            for i in files:
                if i == 'CONTCAR'and file=='CONTCAR':
                    list.append((top,topsz[-1],os.path.join(top,i)))
                elif i.split('.')[-1]=='cif' and file=='cif':
                    list.append((top,i.split('.')[0],os.path.join(top,i)))
                elif i == 'POSCAR' and file == 'POSCAR':
                    list.append((top,topsz[-1],os.path.join(top,i)))
                elif i.split('.')[-1] == 'xyz' and file == 'xyz':
                    list.append((top,i.split('.')[0],os.path.join(top, i)))
        return list

    @staticmethod
    def draw_reaction_path(namelist,valuelist,width=1,space=1.5,offsetup=0.1,offsetdown=0.2,color='#1111DD',draw_value=False):
        newx=[]
        newy=[]
        aa=space
        for i in valuelist:
            newx.append(aa-width/2)
            newx.append(aa+width/2)
            newy.append(i)
            newy.append(i)
            aa+=space
        plt.plot(newx, newy ,color, linestyle=':',linewidth = 2)
        for k in range(0,len(newx),2):
            plt.plot([newx[k],newx[k+1]],[newy[k],newy[k]],color,linewidth = 3)
        plt.xticks([],fontsize=20)  # 嗯调调字体
        plt.yticks(fontsize=20)
        plt.xlabel('Reaction coordinate ', fontsize=25)  # x轴名称
        plt.ylabel('Relative energy (eV)', fontsize=25)  # y轴名称
        for index in range(0,len(newx),2):
            plt.text(newx[index] + width / 2, newy[index] + offsetup, namelist[int(index / 2)], ha='center', va='bottom', fontsize=12)
            if draw_value:
                plt.text(newx[index] + width / 2, newy[index] - offsetdown, newy[index], ha='center', va='bottom', fontsize=12)
        ymax=max(newy)
        ymin=min(newy)
        plt.ylim(ymin-0.5,ymax+0.5)
        plt.gcf().subplots_adjust(left=0.2)
        plt.gcf().subplots_adjust(bottom=0.12)

    @staticmethod
    def save_fig(ylim=None,xlim=None):
        if ylim is not None:
            print(ylim)
            plt.ylim(ylim[0],ylim[1])
        if xlim is not None:
            plt.xlim(xlim[0],xlim[1])
        plt.savefig('rc.png')

    @staticmethod
    def show_fig(ylim=None,xlim=None):
        if ylim is not None:
            print(ylim)
            plt.ylim(ylim[0],ylim[1])
        if xlim is not None:
            plt.xlim(xlim[0],xlim[1])
        plt.show()

if __name__ == '__main__':
    s=Toolkit.kj_to_eV(70)
    Toolkit.draw_reaction_path(['*CH2O','*CH3O','*HCO2','*C2O','*'],[0,2.11,3.21,1.11,2.12],0.5,1)