from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
import numpy as np
from scipy.optimize import root_scalar
import re
def comp(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1]['fluid']
    S1 = nodes.loc[node1]['S']
    H1 = nodes.loc[node1]['H']
    G = nodes.loc[node1]['G']

    H2t = prop("H", "P", P2, "S", S1, fluid)
    H2 = H1 + (H2t - H1) / eff
    T2 = prop('T','H', H2, 'P', P2, fluid)
    S2 = prop('S', 'H', H2, 'P', P2, fluid)
    Q2 = prop('Q', 'H', H2, 'P', P2, fluid)
    nodes.loc[node2] = [T2, P2, H2, S2, Q2, G, fluid]
    blocks.loc[name, 'N'] = G * (H2 - H1)
    pass


def heat(name, node11, node12, node21, node22, T12):
    n = 20
    fluid1 = nodes.loc[node11]['fluid']
    fluid2 = nodes.loc[node21]['fluid']
    H11 = nodes.loc[node11]['H']
    T11 = nodes.loc[node11]['T']
    P11 = nodes.loc[node11]['P']
    P21 = nodes.loc[node21]['P']
    T21 = nodes.loc[node21]['T']
    G1 = nodes.loc[node11]['G']
    H12 = prop("H", "T", T12, "P", P11, fluid1)
    G2 = nodes.loc[node21]['G']
    H21 = nodes.loc[node21]['H']
    step = (H11 - H12) / (n-1)
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    Q = np.zeros(n)
    h11 = H11
    h21 = H21
    t1[0] = T11
    t2[-1] = T21
    for i in range(n-1):
        h12 = h11 - step
        t1[i+1] = prop('T', 'H', h12, 'P', P11, fluid1)
        Q[i+1] = Q[i] + G1*(h11-h12)
        h11 = h12
    for i in range(n-1):
        h22 = h21 + (Q[-1-i] - Q[-2-i])/G2
        t2[-2-i] = prop('T', 'H', h22, 'P', P21, fluid2)
        h21 = h22
    T22 = t2[0]
    H22 = h22
    S12 = prop('S', 'H', H12, 'P', P11, fluid1)
    Q12 = prop('Q', 'H', H12, 'P', P11, fluid1)
    S22 = prop('S', 'H', H22, 'P', P21, fluid2)
    Q22 = prop('Q', 'H', H22, 'P', P21, fluid2)
    nodes.loc[node12] = [T12, P11, H12, S12, Q12, G1, fluid1]
    nodes.loc[node22] = [T22, P21, H22, S22, Q22, G2, fluid2]
    blocks.loc[name, 'Q'] = Q
    blocks.loc[name, 'T1'] = t1
    blocks.loc[name, 'T2'] = t2
    pass

def heat(name, node11, node12, node21, node22, T12):
    n = 20
    fluid1 = nodes.loc[node11]['fluid']
    fluid2 = nodes.loc[node21]['fluid']
    H11 = nodes.loc[node11]['H']
    T11 = nodes.loc[node11]['T']
    P11 = nodes.loc[node11]['P']
    P21 = nodes.loc[node21]['P']
    T21 = nodes.loc[node21]['T']
    G1 = nodes.loc[node11]['G']
    H12 = prop("H", "T", T12, "P", P11, fluid1)
    G2 = nodes.loc[node21]['G']
    H21 = nodes.loc[node21]['H']
    step = (H11 - H12) / (n-1)
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    Q = np.zeros(n)
    h11 = H11
    h21 = H21
    t1[0] = T11
    t2[-1] = T21
    for i in range(n-1):
        h12 = h11 - step
        t1[i+1] = prop('T', 'H', h12, 'P', P11, fluid1)
        Q[i+1] = Q[i] + G1*(h11-h12)
        h11 = h12
    for i in range(n-1):
        h22 = h21 + (Q[-1-i] - Q[-2-i])/G2
        t2[-2-i] = prop('T', 'H', h22, 'P', P21, fluid2)
        h21 = h22
    T22 = t2[0]
    H22 = h22
    S12 = prop('S', 'H', H12, 'P', P11, fluid1)
    Q12 = prop('Q', 'H', H12, 'P', P11, fluid1)
    S22 = prop('S', 'H', H22, 'P', P21, fluid2)
    Q22 = prop('Q', 'H', H22, 'P', P21, fluid2)
    nodes.loc[node12] = [T12, P11, H12, S12, Q12, G1, fluid1]
    nodes.loc[node22] = [T22, P21, H22, S22, Q22, G2, fluid2]
    blocks.loc[name, 'Q'] = Q
    blocks.loc[name, 'T1'] = t1
    blocks.loc[name, 'T2'] = t2
    pass
def heat_exch_2streams(name, node11, node12, node21, node22, **out):
    fluid1 = nodes.loc[node11]['fluid']
    fluid2 = nodes.loc[node21]['fluid']
    G1 = nodes.loc[node11]['G']
    G2 = nodes.loc[node21]['G']
    H11 = nodes.loc[node11]['H']
    H21 = nodes.loc[node21]['H']
    P11 = nodes.loc[node11]['P']
    P21 = nodes.loc[node21]['P']
    P12 = P11 - out['dP1'] if 'dP1' in out else P11
    P22 = P21 - out['dP2'] if 'dP2' in out else P21
    if 'T22' in out:
        T22 = out['T22']
        H22 = prop("H", "T", T22, "P", P22, fluid2)
        S22 = prop('S', 'T', T22, 'P', P22, fluid2)
        Q22 = prop('Q', 'T', T22, 'P', P22, fluid2)
        Q = G2 * (H22 - H21)
        H12 = H11 - Q/G1
        T12 = prop('T', 'H', H12, 'P', P12, fluid1)
        S12 = prop('S', 'H', H12, 'P', P12, fluid1)
        Q12 = prop('Q', 'H', H12, 'P', P12, fluid1)
    if 'T12' in out:
        T12 = out['T12']
        H12 = prop('H', 'T', T12, 'P', P12, fluid1)
        S12 = prop('S', 'T', T12, 'P', P12, fluid1)
        Q12 = prop('Q', 'T', T12, 'P', P12, fluid1)
        Q = G1 * (H11 - H12)
        H22 = H21 + Q/G2
        T22 = prop('T', 'H', H22, 'P', P22, fluid2)
        S22 = prop('S', 'H', H22, 'P', P22, fluid2)
        Q22 = prop('Q', 'H', H22, 'P', P22, fluid2)
    if 'Q22' in out:
        Q22 = out['Q22']
        H22 = prop("H", "Q", Q22, "P", P22, fluid2)
        S22 = prop('S', 'Q', Q22, 'P', P22, fluid2)
        T22 = prop('T', 'Q', Q22, 'P', P22, fluid2)
        Q = G2 * (H22 - H21)
        H12 = H11 - Q/G1
        T12 = prop('T', 'H', H12, 'P', P12, fluid1)
        S12 = prop('S', 'H', H12, 'P', P12, fluid1)
        Q12 = prop('Q', 'H', H12, 'P', P12, fluid1)
    if 'Q12' in out:
        Q12 = out['Q12']
        H12 = prop('H', 'Q', Q12, 'P', P12, fluid1)
        S12 = prop('S', 'Q', Q12, 'P', P12, fluid1)
        T12 = prop('T', 'Q', Q12, 'P', P12, fluid1)
        Q = G1 * (H11 - H12)
        H22 = H21 + Q/G2
        T22 = prop('T', 'H', H22, 'P', P22, fluid2)
        S22 = prop('S', 'H', H22, 'P', P22, fluid2)
        Q22 = prop('Q', 'H', H22, 'P', P22, fluid2)
    nodes.loc[node12] = [T12, P12, H12, S12, Q12, G1, fluid1]
    nodes.loc[node22] = [T22, P22, H22, S22, Q22, G2, fluid2]
    step = 20

    T1 = [prop('T','P',P11-(P11-P12)/step*i,'H',H11-Q/step*i/G1,fluid1) for i in range(step+1)]
    T2 = [prop('T','P',P21-(P21-P22)/step*i,'H',H22-Q/step*i/G2,fluid2) for i in range(step+1)]
    dT = [T1[i] - T2[i] for i in range(step+1)]
    mitta = min(dT)
    blocks.loc[name]['dT'] = mitta
    blocks.loc[name]['T1'] = T1
    blocks.loc[name]['T2'] = T2
    blocks.loc[name]['Q'] = Q
    pass

def mix(name, node11, node12,node2):
    H11 = nodes.loc[node11]['H']
    H12 = nodes.loc[node12]['H']
    G11 = nodes.loc[node11]['G']
    G12 = nodes.loc[node12]['G']
    P11 = nodes.loc[node11]['P']
    P12 = nodes.loc[node12]['P']
    X11 = nodes.loc[node11]['fluid']
    X12 = nodes.loc[node12]['fluid']
    if P11 != P12: print("Давления входящих потоков в ",name," не равны!")
    if X11 != X12: print("Среды входящих потоков в ",name," отличаются!")
    H2 = (G11*H11 + G12*H12)/(G11+G12)
    T2 = prop('T', 'H', H2, 'P', P12, X11)
    S2 = prop('S', 'H', H2, 'P', P12, X11)
    Q2 = prop('Q', 'H', H2, 'P', P12, X11)
    nodes.loc[node2] = [T2, P12, H2, S2, Q2, G11+G12, X11]
    pass
def turb(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1]['fluid']
    S1 = nodes.loc[node1]['S']
    H1 = nodes.loc[node1]['H']
    G = nodes.loc[node1]['G']
    H2t = prop("H", "P", P2, "S", S1, fluid)
    #print(H2t)
    H2 = H1 - (H1 - H2t)*eff
    T2 = prop('T', 'H', H2, 'P', P2, fluid)
    S2 = prop('S', 'H', H2, 'P', P2, fluid)
    Q2 = prop('Q', 'H', H2, 'P', P2, fluid)
    nodes.loc[node2] = [T2, P2, H2, S2, Q2, G, fluid]
    blocks.loc[name, 'N'] = G * (H1 - H2)
    pass

def cond(name, node1, node2):
    P = nodes.loc[node1]['P']
    H1 = nodes.loc[node1]['H']
    fluid = nodes.loc[node1]['fluid']
    G = nodes.loc[node1]['G']
    T2 = prop('T', 'Q', 0, 'P', P, fluid)
    H2 = prop('H', 'Q', 0, 'P', P, fluid)
    S2 = prop('S', 'Q', 0, 'P', P, fluid)
    nodes.loc[node2] = [T2, P, H2, S2, 0, G, fluid]
    blocks.loc[name, 'Q'] = G*(H1 - H2)
    pass

def heat_exch(name, node1, node2,**out):
    P1 = nodes.loc[node1]['P']
    H1 = nodes.loc[node1]['H']
    fluid = nodes.loc[node1]['fluid']
    G = nodes.loc[node1]['G']

    P2 = P1 - out['dP'] if 'dP' in out else P1
    if 'Q' in out:
        Q = out['Q']
        H2 = H1 + Q/G
        T2 = prop('T', 'H', H2, 'P', P2, fluid)
        Q2 = prop('Q', 'H', H2, 'P', P2, fluid)
        S2 = prop('S', 'H', H2, 'P', P2, fluid)
    if 'T' in out:
        T2 = out['T']
        H2 = prop('H', 'T', T2, 'P', P2, fluid)
        Q2 = prop('Q', 'T', T2, 'P', P2, fluid)
        S2 = prop('S', 'T', T2, 'P', P2, fluid)
    if 'x' in out:
        Q2 = out['x']
        H2 = prop('H', 'Q', Q2, 'P', P2, fluid)
        T2 = prop('T', 'Q', Q2, 'P', P2, fluid)
        S2 = prop('S', 'Q', Q2, 'P', P2, fluid)
    nodes.loc[node2] = [T2, P2, H2, S2, Q2, G, fluid]
    blocks.loc[name, 'Q'] = abs(G * (H2 - H1))
    pass
def comb_stoic(name, node11, node12,node2,dP):
    H11 = nodes.loc[node11]['H']
    P11 = nodes.loc[node11]['P']
    F11 = nodes.loc[node11]['fluid']
    F12 = nodes.loc[node12]['fluid']
    H12 = nodes.loc[node12]['H']
    Gox = nodes.loc[node11]['G']
    Gf = nodes.loc[node12]['G']
    P2 = P11 + dP

    Q_CH4 = 55515100
    Q_H2 = 141783257
    Q_CO = 10103390

    Q_CH4L = 50030044
    Q_H2L = 119957537
    Q_COL = 10103390

    M_CH4 = prop('M', 'CH4') * 1000
    M_H2 = prop('M', 'H2') * 1000
    M_CO = prop('M', 'CO') * 1000
    M_O2 = prop('M', 'O2') * 1000
    M_CO2 = prop('M', 'CO2') * 1000
    M_H2O = prop('M', 'H2O') * 1000
    M_N2 = prop('M', 'N2') * 1000
    M_Ar = prop('M', 'Ar') * 1000
    M11 = prop('M', F11) * 1000

    f11num = re.sub('<[^>]+>', ' ', '<'+ F11.replace(']','<').replace('[','>')+'a>').split(' ')[1:-1]
    f11name = re.sub("\[[^]]*\]", '', F11.replace('REFPROP::','')).split('&')
    M12 = prop('M', F12) * 1000
    f12num = re.sub('<[^>]+>', ' ', '<' + F12.replace(']', '<').replace('[', '>') + 'a>').split(' ')[1:-1]
    f12name = re.sub("\[[^]]*\]", '', F12.replace('REFPROP::', '')).split('&')
    Mmol12 = np.zeros(len(f12name))
    m12 = np.zeros(len(f12name))
    G_O2need = np.zeros(len(f12name))
    G_CO2frFuel = np.zeros(len(f12name))
    G_H2OfrFuel = np.zeros(len(f12name))
    for i in range(0,len(f12name)):
        Mmol12[i] = prop('M', f12name[i]) * 1000
        m12[i] = float(Mmol12[i])*float(f12num[i])/M12
    #methane - h2 - co
    if 'Methane' in f12name:
        G_CH4 = Gf*m12[0]
        G_O2need[0] = G_CH4 * (2 * M_O2 / M_CH4)
        G_CO2frFuel[0] = G_CH4 * (M_CO2 / M_CH4)
        G_H2OfrFuel[0] = G_CH4 * (2 * M_H2O / M_CH4)
    if 'H2' in f12name:
        G_H2 = Gf*m12[1]
        G_O2need[1] = G_H2 * (0.5 * M_O2 / M_H2)
        G_H2OfrFuel[1] = G_H2 * (M_H2O / M_H2)
    if 'CO' in f12name:
        G_CO = Gf*m12[2]
        G_O2need[2] = G_CO * (0.5 * M_O2 / M_CO)
        G_CO2frFuel[2] = G_CO * (M_CO2 / M_CO)
    #O2 - 1st
    Mmol11 = np.zeros(len(f11name))
    m11 = np.zeros(len(f11name))
    for i in range(0,len(f11name)):
        Mmol11[i] = prop('M', f11name[i]) * 1000
        m11[i] = float(Mmol11[i])*float(f11num[i])/M11


    G_O2in = m11[f11name.index('O2')]*Gox
    G_O2 = G_O2in - sum(G_O2need)
    G_CO2 = m11[f11name.index('CO2')]*Gox + sum(G_CO2frFuel) if 'CO2' in f11name else sum(G_CO2frFuel)
    G_H2O = m11[f11name.index('H2O')]*Gox + sum(G_H2OfrFuel) if 'H2O' in f11name else sum(G_H2OfrFuel)
    G_N2 = m11[f11name.index('N2')]*Gox if 'N2' in f11name else 0
    G_Ar = m11[f11name.index('Ar')]*Gox if 'Ar' in f11name else 0

    m_N2 = G_N2 / (Gox + Gf)
    m_O2 = G_O2 / (Gox + Gf)
    m_CO2 = G_CO2 / (Gox + Gf)
    m_H2O = G_H2O / (Gox + Gf)
    m_Ar = G_Ar / (Gox + Gf)

    mole_mix = m_N2 / M_N2 + m_CO2 / M_CO2 + m_H2O / M_H2O + m_O2 / M_O2 + m_Ar / M_Ar
    w_N2 = m_N2 / M_N2 / mole_mix
    w_CO2 = m_CO2 / M_CO2 / mole_mix
    w_H2O = m_H2O / M_H2O / mole_mix
    w_O2 = m_O2 / M_O2 / mole_mix
    w_Ar = m_Ar / M_Ar / mole_mix
    fluid = "REFPROP::N2[" + str(float(w_N2)) + "]&CO2[" + str(float(w_CO2)) + "]&H2O[" + str(
        float(w_H2O)) + "]&O2[" + str(float(w_O2)) + "]&Ar[" + str(float(w_Ar))+ "]"
    qCH4 = m12[f12name.index('Methane')]* Q_CH4 if 'Methane' in f12name else 0
    qH2 = m12[f12name.index('H2')] * Q_H2 if 'H2' in f12name else 0
    qCO = m12[f12name.index('CO')] * Q_CO if 'CO' in f12name else 0

    qCH4L = m12[f12name.index('Methane')] * Q_CH4L if 'Methane' in f12name else 0
    qH2L = m12[f12name.index('H2')] * Q_H2L if 'H2' in f12name else 0
    qCOL = m12[f12name.index('CO')] * Q_COL if 'CO' in f12name else 0

    H2 = (Gox * H11 + Gf * (H12 + qCH4 + qH2 + qCO)) / (Gox + Gf)
    T2 = prop('T', 'H', H2, 'P', P2, fluid)
    S2 = prop('S', 'H', H2, 'P', P2, fluid)
    Q2 = prop('Q', 'H', H2, 'P', P2, fluid)
    G2 = Gox + Gf
    nodes.loc[node2] = [T2, P2, H2, S2, Q2, G2, fluid]
    blocks.loc[name, 'Q'] = Gf * (qCH4L + qH2L + qCOL)
    pass




