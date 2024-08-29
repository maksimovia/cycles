import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import Node, comp, turb, heat_exch,comb_stoic2,mix,heat_source,throttle
from scipy.optimize import root_scalar, root
def Sensitivity(PS):
    P2 = 2e6
    T4 = 1200 + 273.15
    P5 = 1e5
    KPDcomp = 0.8
    KPDturb = 0.9

    T1 = 15 + 273.15
    P1 = 1e5
    G1 = 100
    Node('1',G=G1,P=P1,T=T1,fluid='REFPROP::O2[0.21]&N2[0.79]')
    comp('COMP', '1', '2', P2, KPDcomp)
    T3 = 100 + 273.15
    P3 = P2
    Node('3', G=None, P=P3, T=T3, fluid='REFPROP::Methane[1]&H2[0]&CO[0]')
    def Calc0(G3):
        nodes.loc['3', 'G'] = float(G3)
        comb_stoic2('COMB', '2', '3', '4',w_N2 = 0.79)
        return nodes.loc['4','T'] - T4
    root_scalar(Calc0, bracket=[0.2, 3], xtol=10 ** -5, method='bisect')
    turb('TURB', '4', '5', P5, KPDturb)
    P13 = 0.005e6
    P10 = PS
    dT_pp = 30
    P23 = 0.12e6
    dP_deair = 0.05
    T17 = 100 + 273.15
    #dT_econ = 0
    dP_GPK = 0.025
    dP_econ = 0.025
    P9 = 1e5
    T16 = 60 + 273.15
    dT_pinch = 10
    dP_GPK = 0.025
    def Calc(input):
        D = input[0]
        T10 = nodes.loc['5','T'] - dT_pp
        Node('10',G=D,P=P10,T=T10,fluid='Water')
        turb('CVD','10','11',P23*(1+dP_deair),KPDturb)#11
        Gotb = input[1]
        nodes.loc['12'] = nodes.loc['11']  # 12
        nodes.loc['12','G'] = nodes.loc['11','G'] - Gotb
        turb('CND','12','13',P13,KPDturb)#13
        heat_source('COND', '13', '14', x=0)#14
        comp('CONPUMP','14','15',P23*(1+dP_GPK),KPDcomp)#15
        Grec = input[2]
        P17 = P23
        G17 = nodes.loc['15','G'] + Grec
        Node('17',G=G17,P=P17,T=T17,fluid='Water')
        nodes.loc['18'] = nodes.loc['17']#18
        nodes.loc['18','G'] = Grec
        comp('RECPUMP', '18', '19', P23 * (1 + dP_GPK), KPDcomp)#19
        mix('REC', '15', '19', '16')#16
        nodes.loc['20'] = nodes.loc['17'] #20
        nodes.loc['20','G'] = nodes.loc['15','G']
        nodes.loc['21'] = nodes.loc['11']  #21
        nodes.loc['21','G'] = Gotb
        throttle('DROSS','21','22',P23)
        mix('DEAIR','22','20','23')#23
        H23 = prop('H','P',P23,'Q',0,'Water')
        comp('PPUMP','23','24',P10*(1+dP_econ),KPDcomp)#24
        T7 = input[3]
        P7 = P9
        G7 = nodes.loc['5','G']
        Node('7',G=G7,P=P7,T=T7,fluid=nodes.loc['5','fluid'])
        Q26 = 0
        P26 = nodes.loc['24','P']/(1+dP_econ)
        G26 = D
        Node('26',G=G26,P=P26,Q=Q26,fluid='Water')
        heat_exch('ECON', '7', '8', '24', '25', dP2=P10 * (1 + dP_econ) - P10, Q22=0)#25
        T6 = input[4]
        P6 = P9
        G6 = nodes.loc['5','G']
        Node('6',G=G6,P=P6,T=T6,fluid=nodes.loc['5','fluid'])
        heat_exch('EVAP', '6', '7', '26', '27', Q22=1)#27
        nodes.loc['28'] = nodes.loc['27']#28
        heat_exch('PP', '5', '6', '28', '10', T22=T10)#28
        heat_exch('GPK', '8', '9', '16', '17', dP2 =P17 * (1 + dP_GPK) - P17, T22=T17)#8
        Equation1 = T6 - nodes.loc['6','T']
        Equation2 = nodes.loc['23','H'] - H23
        Equation3 = T7 - nodes.loc['7','T']
        Equation4 = T16 - nodes.loc['16','T']
        Equation5 = blocks.loc['EVAP','dT'] - dT_pinch
        return Equation1,Equation2,Equation3,Equation4,Equation5
    root(Calc, x0=([5,0.1,5,500, 600]))
    Q_KU = nodes.loc['5', 'G']*(nodes.loc['5', 'H'] - nodes.loc['9', 'H'])
    Q_COND = blocks.loc['COND','Q']
    N_CVD = blocks.loc['CVD','N']
    N_CND = blocks.loc['CND','N']
    N_PN = blocks.loc['PPUMP','N']
    N_CN = blocks.loc['CONPUMP','N']
    N_RN = blocks.loc['RECPUMP','N']
    N_GASTURB = blocks.loc['TURB','N']
    N_COMP = blocks.loc['COMP','N']
    Q_COMB = blocks.loc['COMB','Q']
    #print('Тепловой баланс:', Q_KU + N_PN + N_CN + N_RN - N_CVD - N_CND - Q_COND)
    KPD_GTU = (N_GASTURB-N_COMP)/Q_COMB * 100
    KPD_KU = (nodes.loc['5','H'] - nodes.loc['9','H'])/(nodes.loc['5','H'] - prop('H','P',P9,'T',T1,nodes.loc['5','fluid'])) * 100
    KPD_PTU = (N_CVD+N_CND-N_PN-N_CN-N_RN)/(Q_KU) * 100
    KPD_PGU = (N_GASTURB-N_COMP+N_CVD+N_CND-N_PN-N_CN-N_RN)/Q_COMB * 100
    # print('КПД ГТУ:', KPD_GTU,'КПД КУ:', KPD_KU,'КПД ПТУ:', KPD_PTU,'КПД ПГУ:', KPD_PGU)
    print(*[blocks.loc['PP', 'Q'] /1e6/ 20 * i for i in range(21)])
    print(*[x - 273.15 for x in blocks.loc['PP','T1']])
    print(*[x - 273.15 for x in blocks.loc['PP', 'T2']])
    print(*[blocks.loc['EVAP', 'Q'] / 1e6 / 20 * i + blocks.loc['PP', 'Q']/1e6  for i in range(21)])
    print(*[x - 273.15 for x in blocks.loc['EVAP', 'T1']])
    print(*[x - 273.15 for x in blocks.loc['EVAP', 'T2']])
    print(*[blocks.loc['ECON', 'Q'] / 1e6 / 20 * i + (blocks.loc['PP', 'Q'] + blocks.loc['EVAP', 'Q'])/1e6  for i in range(21)])
    print(*[x - 273.15 for x in blocks.loc['ECON', 'T1']])
    print(*[x - 273.15 for x in blocks.loc['ECON', 'T2']])
    print(*[blocks.loc['GPK', 'Q'] / 1e6 / 20 * i + (blocks.loc['PP', 'Q'] + blocks.loc['EVAP', 'Q'] + blocks.loc['ECON', 'Q']) / 1e6 for i in range(21)])
    print(*[x - 273.15 for x in blocks.loc['GPK', 'T1']])
    print(*[x - 273.15 for x in blocks.loc['GPK', 'T2']])

    print(PS/1e6,KPD_GTU*100,KPD_KU*100,KPD_PTU*100,KPD_PGU*100,nodes.loc['13','Q']*100,nodes.loc['5','T']-273.15,N_CVD,N_CND,N_GASTURB,N_COMP,N_PN,N_CN,N_RN,nodes.loc['13','Q'])

Sensitivity(3.5e6)
#print(nodes.iloc[:,0:6])
# print(nodes.iloc[:,0:4])
# print(blocks)
# for PS in np.linspace(1.5e6,7e6,20):
#     Sensitivity(PS)