from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
import numpy as np
from scipy.optimize import root_scalar

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

def turb(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1]['fluid']
    S1 = nodes.loc[node1]['S']
    H1 = nodes.loc[node1]['H']
    G = nodes.loc[node1]['G']
    H2t = prop("H", "P", P2, "S", S1, fluid)
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

def comb_stoic_CH4_O2N2(name, node11, node12,node2,dP,O2share):
    Q_comb = 55515100
    P11 = nodes.loc[node11]['P']
    T11 = nodes.loc[node11]['T']
    P12 = nodes.loc[node12]['P']
    T12 = nodes.loc[node12]['T']
    P2 = P11 + dP
    Ga =  nodes.loc[node11]['G']
    Bt = nodes.loc[node12]['G']

    N2share = 1-O2share
    Air = nodes.loc[node11]['fluid']

    M_CH4 = prop('M', 'CH4') * 1000
    M_O2 = prop('M', 'O2') * 1000
    M_CO2 = prop('M', 'CO2') * 1000
    M_H2O = prop('M', 'H2O') * 1000
    M_N2 = prop('M', 'N2') * 1000

    L_O2 = M_O2 / M_CH4 * 2
    L_CO2 = M_CO2 / M_CH4 * 1
    L_H2O = M_H2O / M_CH4 * 2
    L_N2 = M_N2 / M_CH4 * N2share/O2share*2
    L_air = L_N2 + L_O2

    def Comb(T2):
        C_CH4 = prop('Cpmass', 'T', T12, 'P', P12, 'CH4')
        C_CO2 = prop('Cpmass', 'T', T2, 'P', P2, 'CO2')
        C_N2 = prop('Cpmass', 'T', T2, 'P', P2, 'N2')
        C_H2O = prop('Cpmass', 'T', T2, 'P', P2, 'H2O')
        C_airin = prop('Cpmass', 'T', T11, 'P', P11, Air)
        C_airout = prop('Cpmass', 'T', T2, 'P', P2, Air)

        global L_pg,alpha
        alpha = (T2 * (L_CO2 * C_CO2 + L_H2O * C_H2O + L_N2 * C_N2 - L_air * C_airout) - Q_comb - C_CH4 * T12) / (
                T11 * L_air * C_airin - T2 * L_air * C_airout)
        L_pg = (L_CO2 + L_N2 + L_H2O + L_air * (alpha - 1))
        L_pg1 = (Ga+Bt)/Bt
        return L_pg-L_pg1
    T2 = root_scalar(Comb, x0=1473, xtol=10**-3).root
    m_N2 = L_N2 / L_pg[0] + L_air * (alpha[0] - 1) * 0.767 / L_pg[0]
    m_CO2 = L_CO2 / L_pg[0]
    m_H2O = L_H2O / L_pg[0]
    m_O2 = 1 - (m_N2 + m_CO2 + m_H2O)

    mole_mix = m_N2 / M_N2 + m_CO2 / M_CO2 + m_H2O / M_H2O + m_O2 / M_O2
    w_N2 = m_N2 / M_N2 / mole_mix
    w_CO2 = m_CO2 / M_CO2 / mole_mix
    w_H2O = m_H2O / M_H2O / mole_mix
    w_O2 = m_O2 / M_O2 / mole_mix
    fluid="REFPROP::N2[" + str(w_N2) + "]&CO2[" + str(w_CO2) + "]&H2O[" + str(w_H2O) + "]&O2[" + str(w_O2) + "]"

    H2 = prop('H', 'T', T2, 'P', P2, fluid)
    S2 = prop('S', 'T', T2, 'P', P2, fluid)
    Q2 = prop('Q', 'T', T2, 'P', P2, fluid)
    G2 = Ga+Bt
    nodes.loc[node2] = [T2, P2, H2, S2, Q2, G2, fluid]
    blocks.loc[name, 'alpha'] = alpha[0]
    pass

def comb_stoic_CH4_O2N22(name, node11, node12,node2,dP):
    Q_comb = 55515100
    H11 = nodes.loc[node11]['H']
    P11 = nodes.loc[node11]['P']
    H12 = nodes.loc[node12]['H']
    P2 = P11 + dP
    Ga =  nodes.loc[node11]['G']
    Bt = nodes.loc[node12]['G']

    H2 = (Ga*H11+Bt*(H12+Q_comb))/(Ga+Bt)

    M_CH4 = prop('M', 'CH4') * 1000
    M_O2 = prop('M', 'O2') * 1000
    M_CO2 = prop('M', 'CO2') * 1000
    M_H2O = prop('M', 'H2O') * 1000
    M_N2 = prop('M', 'N2') * 1000

    m_N2_air = 0.79*M_N2/(0.79*M_N2+0.21*M_O2)
    m_O2_air = 0.21 * M_O2 / (0.79 * M_N2 + 0.21 * M_O2)
    G_N2 = Ga*m_N2_air
    G_O2need = Bt*(2*M_O2/M_CH4)
    G_O2in = Ga*m_O2_air
    G_O2out = G_O2in - G_O2need
    G_CO2 = Bt*M_CO2/M_CH4
    G_H2O = Bt*2*M_H2O/M_CH4

    m_N2 = G_N2 / (Ga+Bt)
    m_O2 = G_O2out / (Ga + Bt)
    m_CO2 = G_CO2 / (Ga + Bt)
    m_H2O = G_H2O / (Ga + Bt)

    mole_mix = m_N2 / M_N2 + m_CO2 / M_CO2 + m_H2O / M_H2O + m_O2 / M_O2
    w_N2 = m_N2 / M_N2 / mole_mix
    w_CO2 = m_CO2 / M_CO2 / mole_mix
    w_H2O = m_H2O / M_H2O / mole_mix
    w_O2 = m_O2 / M_O2 / mole_mix

    fluid="REFPROP::N2[" + str(float(w_N2)) + "]&CO2[" + str(float(w_CO2)) + "]&H2O[" + str(float(w_H2O)) + "]&O2[" + str(float(w_O2)) + "]"

    T2 = prop('T', 'H', H2, 'P', P2, fluid)
    S2 = prop('S', 'H', H2, 'P', P2, fluid)
    Q2 = prop('Q', 'H', H2, 'P', P2, fluid)
    G2 = Ga+Bt
    nodes.loc[node2] = [T2, P2,H2, S2, Q2, G2, fluid]
    pass

