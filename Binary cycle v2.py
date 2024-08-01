import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb, heat_exch_2streams,comb_stoic,mix,heat_exch
from scipy.optimize import root_scalar
from scipy.optimize import root
def Optimize(PS):
    P3g = 3e6
    dP_KU = 0#0.025
    P10g = 0.1013e6
    P6g = P10g*(1+dP_KU)
    T6g = 1060+273.15
    KPDcomp = 0.85
    KPDturb = 0.9

    X1g = 'REFPROP::O2[0.21]&N2[0.79]'
    T1g = 15+273.15
    P1g = 1e5
    G1g = 490
    H1g = prop("H", "P", P1g, "T", T1g, X1g)
    S1g = prop("S", "P", P1g, "T", T1g, X1g)
    Q1g = prop("Q", "P", P1g, "T", T1g, X1g)
    nodes.loc['1g'] = [T1g, P1g, H1g, S1g, Q1g, G1g, X1g]
    X2g = 'REFPROP::Methane[1]'
    T2g = 50+273.15
    P2g = 1.2e6
    H2g = prop("H", "P", P2g, "T", T2g, X2g)
    S2g = prop("S", "P", P2g, "T", T2g, X2g)
    Q2g = prop("Q", "P", P2g, "T", T2g, X2g)
    nodes.loc['2g'] = [T2g, P2g, H2g, S2g, Q2g, '', X2g]
    def Troot(G2g):
        nodes.loc['2g','G'] = float(G2g)
        comp('AirCOMP', '1g', '3g', P3g, KPDcomp)
        comp('FuelCOMP', '2g', '4g', P3g, KPDcomp)
        comb_stoic('COMB', '3g', '4g','5g',dP=0)
        return nodes.loc['5g']['T'] - T6g
    root_scalar(Troot,x0=1, xtol=10**-5)
    turb('GTURB', '5g', '6g', P6g, KPDturb)

    Pk = 0.004246971e6
    P1_ = PS  # 4.138998301e6
    Pd = 0.12e6
    dP_deair = 0.05
    # dT_deair = 5
    T8_ = 100 + 273.15
    dT_econ = 10
    dP_econ = 0  # 0.25e6
    T7_ = 60 + 273.15
    T1_ = nodes.loc['6g']['T'] - 30#515 + 273.15
    dT_pinch = 10
    dP_PE = 0  # 0.25e6
    dP_GPK = 0  # 0.025
    def Root(input):
        D0_ = input[0]
        Gotb_ = input[1]
        Grec_ = input[2]
        T8g_ = input[3]
        T7g_ = input[4]
        X1 = 'H2O'#1
        T1 = T1_
        P1 = P1_
        G1 = D0_
        H1 = prop("H", "P", P1, "T", T1, X1)
        S1 = prop("S", "P", P1, "T", T1, X1)
        Q1 = prop("Q", "P", P1, "T", T1, X1)
        nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
        turb('CVD','1','2',Pd*(1+dP_deair),0.922722651)#2
        nodes.loc['3'] = nodes.loc['2']#3
        nodes.loc['3','G'] = nodes.loc['2']['G'] - Gotb_
        turb('CND','3','4',Pk,0.790994202)#4
        heat_exch('COND','4','5',x=0)#5
        comp('CONPUMP','5','6',Pd*(1+dP_GPK),KPDcomp)#6
        X8 = 'H2O'#8
        T8 = T8_
        P8 = Pd

        G8 = nodes.loc['6']['G'] + Grec_
        H8 = prop("H", "P", P8, "T", T8, X8)
        S8 = prop("S", "P", P8, "T", T8, X8)
        Q8 = prop("Q", "P", P8, "T", T8, X8)
        nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
        nodes.loc['9'] = nodes.loc['8']#9
        nodes.loc['9','G'] = Grec_
        comp('RECPUMP', '9', '10', Pd * (1 + dP_GPK), 0.85)#10
        mix('REC', '6', '10', '7')#7

        nodes.loc['11'] = nodes.loc['8'] #11
        nodes.loc['11','G'] = nodes.loc['6']['G']
        nodes.loc['12'] = nodes.loc['2']  #12
        nodes.loc['12','G'] = Gotb_
        X13 = 'H2O'#13
        H13 = nodes.loc['2','H']
        P13 = Pd
        G13 = nodes.loc['12','G']
        T13 = prop("T", "P", P13, "H", H13, X13)
        S13 = prop("S", "P", P13, "H", H13, X13)
        Q13 = prop("Q", "P", P13, "H", H13, X13)
        nodes.loc['13'] = [T13, P13, H13, S13, Q13, G13, X13]
        mix('DEAIR','13','11','14')#14
        Hd = prop('H','P',Pd,'Q',0,'H2O')
        comp('PPUMP','14','15',P1+dP_econ+dP_PE,0.85)#15
        X8g = nodes.loc['6g','fluid']#8g
        T8g = T8g_
        P8g = P10g*(1+dP_KU*2/4)
        G8g = nodes.loc['6g','G']
        H8g = prop("H", "P", P8g, "T", T8g, X8g)
        S8g = prop("S", "P", P8g, "T", T8g, X8g)
        Q8g = prop("Q", "P", P8g, "T", T8g, X8g)
        nodes.loc['8g'] = [T8g, P8g, H8g, S8g, Q8g, G8g, X8g]
        X17 = 'H2O' # 17
        Q17 = 0
        P17 = nodes.loc['15','P'] - dP_econ
        G17 = D0_
        H17 = prop("H", "P", P17, "Q", Q17, X17)
        S17 = prop("S", "P", P17, "Q", Q17, X17)
        T17 = prop("T", "P", P17, "Q", Q17, X17)
        nodes.loc['17'] = [T17, P17, H17, S17, Q17, G17, X17]
        heat_exch_2streams('ECON','8g','9g','15','16',T22=T17-dT_econ)#16
        X7g = nodes.loc['6g','fluid']#7g
        T7g = T7g_
        P7g = P10g*(1+dP_KU*3/4)
        G7g = nodes.loc['6g','G']
        H7g = prop("H", "P", P7g, "T", T7g, X7g)
        S7g = prop("S", "P", P7g, "T", T7g, X7g)
        Q7g = prop("Q", "P", P7g, "T", T7g, X7g)
        nodes.loc['7g'] = [T7g, P7g, H7g, S7g, Q7g, G7g, X7g]

        heat_exch_2streams('EVAP','7g','8g','17','18',Q22=1)#18
        nodes.loc['19'] = nodes.loc['18']#19
        heat_exch_2streams('PP','6g','7g','19','1',T22=T1)#19
        heat_exch_2streams('GPK', '9g', '10g', '7', '8', T22=T8_)#8


        eq1 = T7g_ - nodes.loc['7g']['T']
        eq2 = nodes.loc['14']['H'] - Hd
        eq3 = T7_ - nodes.loc['7']['T']
        eq4 = T8g_ - nodes.loc['8g']['T']
        eq5 = T8g_ - nodes.loc['17']['T'] - dT_pinch
        print(eq1,eq2,eq3,eq4,eq5)
        return eq1,eq2,eq3,eq4,eq5
    root(Root, x0=([10,0.1,3,500, 400]), method='hybr',tol=10**-15)

    #print(nodes.iloc[:, 0:6])
    # print(blocks)
    Q_PP = blocks.loc['PP']['Q']
    Q_EVAP = blocks.loc['EVAP']['Q']
    Q_ECON = blocks.loc['ECON']['Q']
    Q_GPK = blocks.loc['GPK']['Q']
    Q_COND = -blocks.loc['COND']['Q']
    N_CVD = blocks.loc['CVD']['N']
    N_CND = blocks.loc['CND']['N']
    N_PN = blocks.loc['PPUMP']['N']
    N_CN = blocks.loc['CONPUMP']['N']
    N_RN = blocks.loc['RECPUMP']['N']
    Ngasturb = blocks.loc['GTURB','N']
    Naircomp = blocks.loc['AirCOMP','N']
    Nfuelcomp = blocks.loc['FuelCOMP','N']
    Qcomb = blocks.loc['COMB','Q']
    # print(Q_PP + Q_EVAP + Q_ECON + Q_GPK + N_PN + N_CN + N_RN - N_CVD - N_CND - Q_COND+(nodes.loc['17']['H']-nodes.loc['16']['H'])*nodes.loc['17']['G'])
    KPD_GTU = (Ngasturb-Naircomp-Nfuelcomp)/Qcomb
    KPD_KU = (nodes.loc['6g','H'] - nodes.loc['10g','H'])/(nodes.loc['6g','H'] - prop('H','P',P10g,'T',T1g,nodes.loc['6g']['fluid']))
    KPD_PTU = (N_CVD+N_CND-N_PN-N_CN-N_RN)/(nodes.loc['6g','G']*(nodes.loc['6g','H'] - nodes.loc['10g','H']))
    KPD_PGU1 = (Ngasturb-Naircomp-Nfuelcomp+N_CVD+N_CND-N_PN-N_CN-N_RN)/Qcomb
    KPD_PGU2 = KPD_GTU + (1-KPD_GTU)*KPD_KU*KPD_PTU
    print(PS,KPD_GTU,KPD_KU,KPD_PTU,KPD_PGU1,KPD_PGU2,nodes.loc['4']['Q'],nodes.loc['6g']['T']-273.15,N_CVD,N_CND,Ngasturb,Naircomp,Nfuelcomp,N_PN,N_CN,N_RN)

Optimize(4.5e6)
# for PS in np.linspace(3e6,7e6,10):
#     Optimize(PS)