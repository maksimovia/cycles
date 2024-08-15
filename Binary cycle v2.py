import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb, heat_exch_2streams,comb_stoic,mix,heat_exch
from scipy.optimize import root_scalar
from scipy.optimize import root
def Sensitivity(PS):
    P2g = 2.5e6
    P6g = 1e5
    T5g = 1200 + 273.15
    KPDcomp = 0.8
    KPDturb = 0.9

    X1g = 'REFPROP::O2[0.2]&N2[0.77]&CO2[0.01]&Ar[0.01]&H2O[0.01]'
    T1g = 15 + 273.15
    P1g = 1e5
    G1g = 120
    H1g = prop("H", "P", P1g, "T", T1g, X1g)
    S1g = prop("S", "P", P1g, "T", T1g, X1g)
    Q1g = prop("Q", "P", P1g, "T", T1g, X1g)
    nodes.loc['1g'] = [T1g, P1g, H1g, S1g, Q1g, G1g, X1g]
    X3g = 'REFPROP::Methane[1]&H2[0]&CO[0]'
    T3g = 50 + 273.15
    P3g = 1.2e5
    G3g = 1
    H3g = prop("H", "P", P3g, "T", T3g, X3g)
    S3g = prop("S", "P", P3g, "T", T3g, X3g)
    Q3g = prop("Q", "P", P3g, "T", T3g, X3g)
    nodes.loc['3g'] = [T3g, P3g, H3g, S3g, Q3g, G3g, X3g]

    def Calc0(G3):
        nodes.loc['3g', 'G'] = float(G3)
        comp('AirCOMP', '1g', '2g', P2g, KPDcomp)
        comp('FuelCOMP', '3g', '4g', P2g, KPDcomp)
        comb_stoic('COMB', '2g', '4g', '5g', dP=0)
        return nodes.loc['5g']['T'] - T5g

    root_scalar(Calc0, bracket=[0.2, 5], xtol=10 ** -9, method='bisect')
    turb('GTURB', '5g', '6g', P6g, KPDturb)

    P4 = 0.005e6
    P1 = PS
    dT_pp = 30
    P14 = 0.12e6
    dP_deair = 0.05
    T8 = 100 + 273.15
    dT_econ = 10
    dP_GPK = 0.025
    dP_econ = 0.025
    P10g = 1e5
    T7 = 60 + 273.15

    dT_pinch = 10
    dP_KU = 0
    dP_PE = 0  # 0.25e6
    dP_GPK = 0.025
    def Calc(input):
        D = input[0]
        X1 = 'Water'#1
        T1 = nodes.loc['6g']['T'] - dT_pp
        G1 = D
        H1 = prop("H", "P", P1, "T", T1, X1)
        S1 = prop("S", "P", P1, "T", T1, X1)
        Q1 = prop("Q", "P", P1, "T", T1, X1)
        nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
        turb('CVD','1','2',P14*(1+dP_deair),KPDturb)#2

        Gotb = input[1]
        nodes.loc['3'] = nodes.loc['2']  # 3
        nodes.loc['3','G'] = nodes.loc['2']['G'] - Gotb
        turb('CND','3','4',P4,KPDturb)#4
        heat_exch('COND','4','5',x=0)#5
        comp('CONPUMP','5','6',P14*(1+dP_GPK),KPDcomp)#6
        Grec = input[2]
        X8 = 'Water'#8
        P8 = P14
        G8 = nodes.loc['6']['G'] + Grec
        H8 = prop("H", "P", P8, "T", T8, X8)
        S8 = prop("S", "P", P8, "T", T8, X8)
        Q8 = prop("Q", "P", P8, "T", T8, X8)
        nodes.loc['8'] = [T8, P8, H8, S8, Q8, G8, X8]
        nodes.loc['9'] = nodes.loc['8']#9
        nodes.loc['9','G'] = Grec
        comp('RECPUMP', '9', '10', P14 * (1 + dP_GPK), KPDcomp)#10
        mix('REC', '6', '10', '7')#7
        nodes.loc['11'] = nodes.loc['8'] #11
        nodes.loc['11','G'] = nodes.loc['6']['G']
        nodes.loc['12'] = nodes.loc['2']  #12
        nodes.loc['12','G'] = Gotb
        X13 = 'Water'#13
        H13 = nodes.loc['2','H']
        P13 = P14
        G13 = nodes.loc['12','G']
        T13 = prop("T", "P", P13, "H", H13, X13)
        S13 = prop("S", "P", P13, "H", H13, X13)
        Q13 = prop("Q", "P", P13, "H", H13, X13)
        nodes.loc['13'] = [T13, P13, H13, S13, Q13, G13, X13]
        mix('DEAIR','13','11','14')#14
        H14 = prop('H','P',P14,'Q',0,'Water')
        comp('PPUMP','14','15',P1*(1+dP_econ),KPDcomp)#15
        T8g = input[3]
        X8g = nodes.loc['6g','fluid']#8g
        P8g = P10g
        G8g = nodes.loc['6g','G']
        H8g = prop("H", "P", P8g, "T", T8g, X8g)
        S8g = prop("S", "P", P8g, "T", T8g, X8g)
        Q8g = prop("Q", "P", P8g, "T", T8g, X8g)
        nodes.loc['8g'] = [T8g, P8g, H8g, S8g, Q8g, G8g, X8g]
        X17 = 'Water' # 17
        Q17 = 0
        P17 = nodes.loc['15','P']/(1+dP_econ)
        G17 = D
        H17 = prop("H", "P", P17, "Q", Q17, X17)
        S17 = prop("S", "P", P17, "Q", Q17, X17)
        T17 = prop("T", "P", P17, "Q", Q17, X17)
        nodes.loc['17'] = [T17, P17, H17, S17, Q17, G17, X17]
        heat_exch_2streams('ECON','8g','9g','15','16',dP2=P1*(1+dP_econ)-P1,T22=T17-dT_econ)#16
        X7g = nodes.loc['6g','fluid']#7g
        T7g = input[4]
        P7g = P10g
        G7g = nodes.loc['6g','G']
        H7g = prop("H", "P", P7g, "T", T7g, X7g)
        S7g = prop("S", "P", P7g, "T", T7g, X7g)
        Q7g = prop("Q", "P", P7g, "T", T7g, X7g)
        nodes.loc['7g'] = [T7g, P7g, H7g, S7g, Q7g, G7g, X7g]

        heat_exch_2streams('EVAP','7g','8g','17','18',Q22=1)#18
        nodes.loc['19'] = nodes.loc['18']#19
        heat_exch_2streams('PP','6g','7g','19','1',T22=T1)#19
        heat_exch_2streams('GPK', '9g', '10g', '7', '8',dP2 = P8*(1+dP_GPK)-P8, T22=T8)#8


        eq1 = T7g - nodes.loc['7g']['T']
        eq2 = nodes.loc['14']['H'] - H14
        eq3 = T7 - nodes.loc['7']['T']
        eq4 = T8g - nodes.loc['8g']['T']
        eq5 = blocks.loc['EVAP','dT'] - dT_pinch
        #print(eq1,eq2,eq3,eq4,eq5)
        #print(input)
        return eq1,eq2,eq3,eq4,eq5
    root(Calc, x0=([5,0.1,5,500, 600]), method='hybr',tol=10**-5)

    #print(nodes.iloc[:, 0:6])
    # print(blocks)
    Q_KU = nodes.loc['6g', 'G']*(nodes.loc['6g', 'H'] - nodes.loc['10g', 'H'])
    Q_BAR = nodes.loc['17', 'G'] * (nodes.loc['17', 'H'] - nodes.loc['16', 'H'])
    Q_PP = blocks.loc['PP']['Q']
    Q_EVAP = blocks.loc['EVAP']['Q']
    Q_ECON = blocks.loc['ECON']['Q']
    Q_GPK = blocks.loc['GPK']['Q']
    Q_COND = blocks.loc['COND']['Q']
    N_CVD = blocks.loc['CVD']['N']
    N_CND = blocks.loc['CND']['N']
    N_PN = blocks.loc['PPUMP']['N']
    N_CN = blocks.loc['CONPUMP']['N']
    N_RN = blocks.loc['RECPUMP']['N']
    Ngasturb = blocks.loc['GTURB','N']
    Naircomp = blocks.loc['AirCOMP','N']
    Nfuelcomp = blocks.loc['FuelCOMP','N']
    Qcomb = blocks.loc['COMB','Q']
    # print(Q_KU,Q_PP + Q_EVAP + Q_ECON + Q_GPK + (nodes.loc['17']['H']-nodes.loc['16']['H'])*nodes.loc['17']['G'])
    # print(Q_PP, Q_EVAP, Q_ECON, Q_GPK, N_PN, N_CN, N_RN, N_CVD, N_CND, Q_COND,(nodes.loc['17']['H']-nodes.loc['16']['H'])*nodes.loc['17']['G'])
    # print(Q_PP + Q_EVAP + Q_ECON + Q_GPK + N_PN + N_CN + N_RN - N_CVD - N_CND - Q_COND + (nodes.loc['17']['H']-nodes.loc['16']['H'])*nodes.loc['17']['G'])
    # print(Q_KU + N_PN + N_CN + N_RN - N_CVD - N_CND - Q_COND + (nodes.loc['17']['H']-nodes.loc['16']['H'])*nodes.loc['17']['G'])
    #print('Тепловой баланс:', Q_KU + Q_BAR + N_PN + N_CN + N_RN - N_CVD - N_CND - Q_COND)
    KPD_GTU = (Ngasturb-Naircomp-Nfuelcomp)/Qcomb
    KPD_KU = (nodes.loc['6g','H'] - nodes.loc['10g','H'])/(nodes.loc['6g','H'] - prop('H','P',P1g,'T',T1g,nodes.loc['6g']['fluid']))
    KPD_PTU = (N_CVD+N_CND-N_PN-N_CN-N_RN)/(Q_KU)
    KPD_PGU1 = (Ngasturb-Naircomp-Nfuelcomp+N_CVD+N_CND-N_PN-N_CN-N_RN)/Qcomb
    KPD_PGU2 = KPD_GTU + (1-KPD_GTU)*KPD_KU*KPD_PTU
    # print(*[blocks.loc['PP', 'Q'] /1e6/ 20 * i for i in range(21)])
    # print(*[x - 273.15 for x in blocks.loc['PP','T1']])
    # print(*[x - 273.15 for x in blocks.loc['PP', 'T2']])
    # print(*[blocks.loc['EVAP', 'Q'] / 1e6 / 20 * i + blocks.loc['PP', 'Q']/1e6  for i in range(21)])
    # print(*[x - 273.15 for x in blocks.loc['EVAP', 'T1']])
    # print(*[x - 273.15 for x in blocks.loc['EVAP', 'T2']])
    # print(*[blocks.loc['ECON', 'Q'] / 1e6 / 20 * i + (blocks.loc['PP', 'Q'] + blocks.loc['EVAP', 'Q'])/1e6  for i in range(21)])
    # print(*[x - 273.15 for x in blocks.loc['ECON', 'T1']])
    # print(*[x - 273.15 for x in blocks.loc['ECON', 'T2']])
    # print(*[blocks.loc['GPK', 'Q'] / 1e6 / 20 * i + (blocks.loc['PP', 'Q'] + blocks.loc['EVAP', 'Q'] + blocks.loc['ECON', 'Q']) / 1e6 for i in
    #         range(21)])
    # print(*[x - 273.15 for x in blocks.loc['GPK', 'T1']])
    # print(*[x - 273.15 for x in blocks.loc['GPK', 'T2']])
    print(PS,KPD_GTU,KPD_KU,KPD_PTU,KPD_PGU1,KPD_PGU2,nodes.loc['4']['Q'],nodes.loc['6g']['T']-273.15,N_CVD,N_CND,Ngasturb,Naircomp,Nfuelcomp,N_PN,N_CN,N_RN,nodes.loc['4','Q'])

Sensitivity(3.5e6)
print(nodes.iloc[:,5])
print(blocks)
# for PS in np.linspace(3e6,7e6,10):
#     Sensitivity(PS)