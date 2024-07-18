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


def heat(name, node11, node12, node21, node22, T12, dT):
    n = 51

    fluid1 = nodes.loc[node11]['fluid']
    fluid2 = nodes.loc[node21]['fluid']

    H11 = nodes.loc[node11]['H']
    T11 = nodes.loc[node11]['T']
    P11 = nodes.loc[node11]['P']
    P21 = nodes.loc[node21]['P']

    G1 = nodes.loc[node11]['G']
    H12 = prop("H", "T", T12, "P", P11, fluid1)
    G2 = nodes.loc[node21]['G']
    H21 = nodes.loc[node21]['H']

    step = (H11 - H12) / (n-1)
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    Q = np.zeros(n)
    print(step)
    h11 = H11
    h21 = H21
    t1[0] = T11
    for i in range(n-1):
        h12 = h11 - step
        t1[i+1] = prop('T', 'H', h12, 'P', P11, fluid1)
        Q[i+1] = Q[i] + G1*(h11-h12)
        print(i)
    for i in range(n-1):
        h22 = h21 + (Q[n-1 - i] - Q[n-1 - i - 1])/G2
        t2[n-1-i] = prop('T', 'H', h22, 'P', P21, fluid2)




#         step = (self.H11-self.H12)/self.hsteps
#         t1 = np.zeros(self.hsteps + 1)
#         t2 = np.zeros(self.hsteps + 1)
#         Q = np.zeros(self.hsteps + 1)
#         h11 = self.H11
#         h21 = self.H21
#         for i in range(self.hsteps + 1):
#             t1[i] = prop.h_p(h11, self.P11, self.fluid1)["T"]
#             if i < self.hsteps:
#                 h12 = h11 - step
#                 dQ = self.G1 * (h11 - h12)
#                 h11 = h12
#                 Q[i + 1] = Q[i] + dQ
#         for i in range(self.hsteps + 1):
#             t2[self.hsteps - i] = prop.h_p(h21, p2, self.fluid2)["T"]
#             if i < self.hsteps:
#                 h22 = h21 + (Q[self.hsteps - i] - Q[self.hsteps - i - 1]) * self.dQ / G2
#                 h21 = h22
#                 p2 = p2 - dP2 / (self.hsteps)
#         DT = t1 - t2
#         min_dt = min(DT[:-1])
#         T12 = t1[-1]
#         S12 = prop.h_p(h12, self.P11, self.fluid1)["S"]
#         Q12 = prop.h_p(h12, self.P11, self.fluid1)["Q"]
#         S22 = prop.h_p(h22, p2, self.fluid2)["S"]
#         Q22 = prop.h_p(h22, p2, self.fluid2)["Q"]
#         self.streams.loc[self.stream12] = [T12, self.P11, h12, S12, Q12, self.G1, self.fluid1]
#         self.streams.loc[self.stream22] = [t2[0], p2, h22, S22, Q22, G2, self.fluid2]
#         self.blocks.loc['HEAT'] = [0, Q[-1], min_dt, t1, t2, Q]
#         pass