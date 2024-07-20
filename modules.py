from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
import numpy as np

# Функция сжатия
def comp(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1]['fluid']
    S1 = nodes.loc[node1]['S']
    H1 = nodes.loc[node1]['H']
    G = nodes.loc[node1]['G']

    H2t = prop("H", "P", P2, "S", S1, fluid)
    H2 = H1 + (H2t - H1) / eff
    T2 = prop('T', 'H', H2, 'P', P2, fluid)
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
