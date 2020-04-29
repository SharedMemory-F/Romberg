# -*- coding:utf-8 -*-
"""
@author      : zhang yifan
@data        : 2020-04-01 10:27
@brief       : 利用Romberg算法求解数值积分
"""
import numpy as np 
from sympy import *

class Romberg:
    def __init__(self, deviation, interval):
        self.deviation = deviation #误差上限
        self.interval = interval
        #self.XY = np.array([[]])
        self.table = np.array([[]])
        self.ratio = np.array([[4.0/3, 1.0/3], [16.0/15, 1.0/15], [64.0/63, 1.0/63]])


    def function(self, x):
        #return np.sin(x)/x
        return 10/pow(x,2)*np.sin(10/x)

    def divide(self, n):#分割区间
        a = self.interval[0]
        b = self.interval[1]
        delta = (b-a)/n 
        # 区间左端点a特殊处理[0,1]
        #np_XY0 = np.array([[0,1]])
        np_X = np.array([a+i*delta for i in range(0, n+1)]).reshape(n+1,1)
        np_Y = self.function(np_X)
        return np.hstack((np_X,np_Y))

    
    def get_int_value(self, n):#复化梯形求积公式
        XY = self.divide(n)
        a = self.interval[0]
        b = self.interval[1]
        return 0.5*((b-a)/n)*(XY[0][1]+XY[-1][1]+2*np.sum([XY[i][1] for i in range(1,n)]))

    def save_table(self):
        # 更新第0列的T
        #self.table[:][0] = [self.get_int_value(pow(2,i)) for i in range(5)]
        self.table = np.array([self.get_int_value(pow(2,i)) for i in range(5)]).reshape(5,1) #T0
        # print("2"+table_T0)
        # 根据计算出的T0-T0(4)计算对应的T1、T2、T3
        for i in range(3):
            table_T = np.zeros(shape=(5,1))
            for j in range(i,4):
                table_T[j+1]= self.table[j+1][i]*self.ratio[i][0] - self.table[j][i]*self.ratio[i][1]
            self.table = np.hstack((self.table,table_T))
       # print(self.table)

    def update_table(self):
        self.save_table()
        i = 5
        while(np.abs(self.table[-1][-1]-self.table[-2][-1]) > 255*self.deviation):
            #从第6行继续更新table[i][0], 初始化第6行
            table_C = np.zeros(shape = (1, 4))
            self.table = np.append(self.table, table_C, axis=0)
            #print(self.table)
            self.table[i][0] = self.get_int_value(pow(2,i))
            for j in range(3):
                self.table[-1][j+1] = self.table[i][j]*self.ratio[j][0] - self.table[i-1][j]*self.ratio[j][1]
            i = i+1
        # if debug == 'True':
        #     print(self.table)
    
    def print_table(self):
        print("\n----------------Romberg算法的T表如下-------------------------------------------------\n")
        for i in self.table:
            print([j for j in i if j!=0])
        print("\n----------------Romberg算法结果如下--------------------------------------------------\n")
        loss = abs((self.table[-1][-1]-self.table[-2][-1])/255)
        print("Romberg算法理论误差约为 %.6e \nRomberg算法结果为： %.8f"%(loss, self.table[-1][-1]))
        x = symbols('x')
        actrue = integrate(10/pow(x,2)*sin(10/x),(x,2,5))
        print("真实积分为%s\n真实积分约为%.8f"%(actrue,actrue))
        print("与真实积分的误差约为：%.6e"%((actrue)-self.table[-1][-1]))
        print("\n--------------------------------------------------------------------------------------")

if __name__ == "__main__":
    romberg = Romberg(0.0001, [2, 5])
    romberg.update_table()
    romberg.print_table()
    # XY = np.array([[1,2],[1,2,3]])
    # print(XY)
    #romberg.update_table()
    #romberg.print_table()
    # x = symbols('x')
    # print(integrate(10/pow(x,2)*sin(10/x),(x,2,5)))
    