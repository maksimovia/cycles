import numpy as np
from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, turb, heat_exch_2streams,comb_stoic,mix,heat_exch
from scipy.optimize import root_scalar
P3 = 3e6
dP_KU = 0.025
P25 = 1e5
P6 = P25*(1+dP_KU)
T6 = 1060+273.15
KPDcomp = 0.8
KPDturb = 0.9

X1 = 'REFPROP::O2[0.20]&N2[0.77]&CO2[0.01]&H2O[0.01]&Ar[0.01]'
T1 = 15+273.15
P1 = 1e5
G1 = 120
H1 = prop("H", "P", P1, "T", T1, X1)
S1 = prop("S", "P", P1, "T", T1, X1)
Q1 = prop("Q", "P", P1, "T", T1, X1)
nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
def T_CC(G2):
    X2 = 'REFPROP::Methane[1]&H2[0]&CO[0]'
    T2 = 50+273.15
    P2 = 1.2e6
    G2 = float(G2)
    H2 = prop("H", "P", P2, "T", T2, X2)
    S2 = prop("S", "P", P2, "T", T2, X2)
    Q2 = prop("Q", "P", P2, "T", T2, X2)
    nodes.loc['2'] = [T2, P2, H2, S2, Q2, G2, X2]
    comp('AirCOMP', '1', '3', P3, KPDcomp)
    comp('FuelCOMP', '2', '4', P3, KPDcomp)
    comb_stoic('COMB', '3', '4','5',dP=0)
    return nodes.loc['5']['T'] - T6
root_scalar(T_CC,x0=1, xtol=10**-5)
turb('GTURB', '5', '6', P6, KPDturb)
print(nodes.loc['6','fluid'])
X2 = 'REFPROP::N2[.756199238198073]&O2[.139736653897233]&WATER[6.31134721631234E-02]&CO2[3.19053756565723E-02]&AR[9.04526008499821E-03]'
T2 = 543+273.15
P2 = 0.1013e6
G2 = 502
H2 = prop("H", "P", P2, "T", T2, X2)
S2 = prop("S", "P", P2, "T", T2, X2)
Q2 = prop("Q", "P", P2, "T", T2, X2)
nodes.loc['6'] = [T2, P2, H2, S2, Q2, G2, X2]

Pk = 0.004246971e6
P10 = 4.138998301e6
Pd = 0.12e6
dP_deair = 0.3
dT_deair = 5
dT_econ = 10
dP_econ = 0#0.25e6
T15 = 60 + 273.15
T10 = 515 + 273.15
dT_pinch = 15
dP_PE2 = 0#0.25e6
dP_GPK = 0#0.025

X22 = 'Water'
P22 = P10 + dP_PE2
Q22 = 0
T22 = prop('T','Q',Q22,'P',P22,X22)
H22 = prop('H','Q',Q22,'P',P22,X22)
S22 = prop('S','Q',Q22,'P',P22,X22)
nodes.loc['22'] = [T22, P22, H22, S22, Q22, 0, X22]
P23 = P22
Q23 = 1
T23 = prop('T','Q',Q23,'P',P23,X22)
H23 = prop('H','Q',Q23,'P',P23,X22)
S23 = prop('S','Q',Q23,'P',P23,X22)
nodes.loc['23'] = [T23, P23, H23, S23, Q23, 0, X22]
def Groot(G):
    nodes.loc['22','G'] = G
    nodes.loc['24'] = nodes.loc['23']
    heat_exch_2streams('PP','6', '7', '24', '10',dP1 = P25*dP_KU/4,dP2 = dP_PE2, T22 = T10)
    heat_exch_2streams('EVAP','7', '8', '22', '23',dP1 = P25*dP_KU/4, Q22 = 1)
    return nodes.loc['8']['T'] - nodes.loc['22']['T'] - dT_pinch
root_scalar(Groot,x0=4, xtol=10**-5)
#18
X18 = nodes.loc['22']['fluid']
P18 = Pd
Q18 = 0
G18 = nodes.loc['22']['G']
T18 = prop('T','Q',Q18,'P',P18,X18)
H18 = prop('H','Q',Q18,'P',P18,X18)
S18 = prop('S','Q',Q18,'P',P18,X18)
nodes.loc['18'] = [T18, P18, H18, S18, Q18, G18, X18]
#19
comp('PUMP','18','19',P22+dP_econ,KPDcomp)
#20
nodes.loc['20'] = nodes.loc['19']
#21
heat_exch_2streams('ECON','8','9','20','21',dP1 = P25*dP_KU/4,dP2=dP_econ,T22 = nodes.loc['22']['T'] - dT_econ)
#17
T17 = nodes.loc['18']['T'] - dT_deair
P17 = Pd
X17 = nodes.loc['18']['fluid']
Q17 = prop('Q','T',T17,'P',P17,X17)
H17 = prop('H','T',T17,'P',P17,X17)
S17 = prop('S','T',T17,'P',P17,X17)
nodes.loc['17'] = [T17, P17, H17, S17, Q17, 1, X17]
#11
turb('CVD','10','11',Pd*(1+dP_deair),KPDturb)
nodes.loc['11','G'] = nodes.loc['10']['G']
#11.5
X115 = nodes.loc['11']['fluid']
H115 = nodes.loc['11']['H']
P115 = Pd
T115 = prop('T','H',H115,'P',P115,X115)
H115 = prop('H','H',H115,'P',P115,X115)
S115 = prop('S','H',H115,'P',P115,X115)
Q115 = prop('Q','H',H115,'P',P115,X115)
nodes.loc['11.5'] = [T115, P115, H115, S115, Q115, 1, X115]
def Groot2(Gd):
    nodes.loc['17','G'] = nodes.loc['10']['G'] - Gd
    nodes.loc['11.5','G'] = Gd
    #print(nodes.loc['17']['G'],nodes.loc['11.5']['G'])
    H1 = nodes.loc['18']['H']
    mix('DEAR','17','11.5','18')
    return nodes.loc['18']['H'] - H1
root_scalar(Groot2,x0 = 1,xtol=10**-5)
#116
nodes.loc['11.6'] = nodes.loc['11']
nodes.loc['11.6','G'] = nodes.loc['11']['G'] - nodes.loc['11.5']['G']
#12
turb('CND','11.6','12',Pk,KPDturb)
#13
heat_exch('COND','12','13',x=0)
#14
comp('CONPUMP','13','14',Pd*(1+dP_GPK),KPDcomp)
def Groot3(Grec):
    nodes.loc['16'] = nodes.loc['17']
    nodes.loc['16', 'G'] = Grec + nodes.loc['14']['G']
    # 16.1
    nodes.loc['16.1'] = nodes.loc['17']
    nodes.loc['16.1', 'G'] = Grec
    # 16.2
    comp('RECPUMP', '16.1', '16.2', Pd * (1 + dP_GPK), KPDcomp)
    # 15
    mix('REC', '14', '16.2', '15')
    # 16
    heat_exch_2streams('GPK', '9', '10.5', '15', '16',dP1=P25*dP_KU/4, T22=nodes.loc['17']['T'])
    return nodes.loc['15']['T'] - T15
root_scalar(Groot3,x0=1,xtol=10**-5)


NGturb = blocks.loc['GTURB','N']
NGcomp = (blocks.loc['AirCOMP','N']+blocks.loc['FuelCOMP','N'])
Q = blocks.loc['COMB','Q']
GKU = nodes.loc['6']['G']
HinKU = nodes.loc['6']['H']
HoutKU = nodes.loc['10.5']['H']

Q_PP = blocks.loc['PP']['Q']
Q_EVAP = blocks.loc['EVAP']['Q']
Q_ECON= blocks.loc['ECON']['Q']
Q_GPK = blocks.loc['GPK']['Q']
Q_COND = -blocks.loc['COND']['Q']
N_CVD = blocks.loc['CVD']['N']
N_CND = blocks.loc['CND']['N']
N_PN = blocks.loc['PUMP']['N']
N_CN = blocks.loc['CONPUMP']['N']
N_RN = blocks.loc['RECPUMP']['N']
print(Q_PP+Q_EVAP+Q_ECON+Q_GPK+N_PN+N_CN+N_RN-N_CVD-N_CND-Q_COND)
KPD_GTU = (NGturb-NGcomp)/Q*100
KPD_PTU =(N_CVD+N_CND+N_PN+N_CN+N_RN)/(GKU*(HinKU-HoutKU))*100
print(KPD_GTU,KPD_PTU)

#TQ
T_gor = np.array([nodes.loc['6','T'],nodes.loc['7','T'],nodes.loc['8','T'],nodes.loc['9','T'],nodes.loc['10.5','T']])
T_hol = np.array([nodes.loc['10','T'],nodes.loc['24','T'],nodes.loc['22','T'],nodes.loc['21','T'],nodes.loc['20','T'],nodes.loc['16','T'],nodes.loc['15','T']])
Q_gorhol = np.array([0,Q_PP,Q_EVAP,Q_ECON,Q_GPK])
# print(T_gor)
# print(T_hol)
# print(Q_gorhol)

print(nodes.iloc[:,0:6])
print(blocks)