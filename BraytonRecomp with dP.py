import numpy as np
import scipy.optimize
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb,heat_source,heat_exch,split,mix
from scipy.optimize import root_scalar,root,minimize


global data
data = open("res.txt", "a")

def Sensitivity(inp):
    P1 = inp[0]
    x = inp[1]
    P7 = inp[2]
    dP = 50e3
    G1 = 1
    T7 = 32+273.15
    T2 = 505+273.15
    KPDcomp = 0.9
    KPDturb = 0.9
    dTh = 10
    dTl = 10
    def Calc(input):
        T1 = input[0]
        Node('1',G=G1,P=P1-2*dP,T=T1,fluid='CO2')
        heat_source('HEAT', '1', '2', T=T2,dP=dP)
        turb('TURB','2','3',P7+3*dP,KPDturb)
        T12 = input[1]
        Node('12',G=G1,P=P1-dP,T=T12,fluid='CO2')
        heat_exch('HTHE', '3', '4', '12', '1', T22=T1,dP1=dP,dP2=dP)
        T5 = input[2]
        Node('5',G=G1,P=P7+dP,T=T5,fluid='CO2')
        T8 = input[3]
        Node('8',G=G1*x,P=P1,T=T8,fluid='CO2')
        heat_exch('LTHE', '4', '5', '8', '9', T12=T5,dP1=dP,dP2=dP)
        split('SPLIT','5','6','10',x)
        heat_source('COOL', '6', '7', T=T7,dP=dP)
        comp('MCOMP','7','8',P1,KPDcomp)
        comp('ACOMP','10','11',P1-dP,KPDcomp)
        mix('MIX','9','11','12')
        Equation1 = nodes.loc['8','T'] - T8
        Equation2 = nodes.loc['12','T'] - T12
        Equation3 = blocks.loc['HTHE', 'dT'] - dTh
        Equation4 = blocks.loc['LTHE','dT'] - dTl
        return Equation3,Equation2,Equation1,Equation4
    root(Calc,x0=[800, 500, 400, 350])
    N_TURB = blocks.loc['TURB']['N']
    N_MCOMP = blocks.loc['MCOMP']['N']
    N_RCOMP = blocks.loc['ACOMP']['N']
    Q_COND = blocks.loc['COOL','Q']
    Q_HEAT = blocks.loc['HEAT','Q']
    KPD1 = (N_TURB - N_MCOMP - N_RCOMP)/Q_HEAT*100
    # print(Q_HEAT + N_MCOMP + N_RCOMP - N_TURB - Q_COND)
    print(P1/1e6,x,P7,KPD1)
    data.write( '\n'+ str(P1/1e6)
               +'\t'+ str(x)
               +'\t'+ str(P7)
               +'\t'+ str(KPD1)
               +'\t'+ str(Q_HEAT + N_MCOMP + N_RCOMP - N_TURB - Q_COND))
    # print(nodes)

    return -KPD1

from multiprocessing import Process, active_children
from time import sleep

if __name__ == '__main__':
    for P1 in np.linspace(20e6,35e6,16):
        for x in np.linspace(0.6,0.8,21):
            for P7 in np.linspace(7.6e6,8.2e6,11):
                calculating = Process(target=Sensitivity, args=([P1,x,P7]))
                calculating.start()
                while len(active_children()) > 59:
                    sleep(2)
            data.close()



