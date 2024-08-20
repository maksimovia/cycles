from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
import numpy as np
from scipy.optimize import root_scalar
import re
import CoolProp

def Node(node,G,fluid,**out):
    if 'P' in out and 'T' in out:
        nodes.loc[node, 'T'] = out['T']
        nodes.loc[node, 'P'] = out['P']
        nodes.loc[node, 'H'] = prop('H', 'T', out['T'], 'P', out['P'], fluid)
        nodes.loc[node, 'S'] = prop('S', 'T', out['T'], 'P', out['P'], fluid)
        nodes.loc[node, 'Q'] = prop('Q', 'T', out['T'], 'P', out['P'], fluid)
    elif 'P' in out and 'H' in out:
        nodes.loc[node, 'H'] = out['H']
        nodes.loc[node, 'P'] = out['P']
        nodes.loc[node, 'T'] = prop('T', 'H', out['H'], 'P', out['P'], fluid)
        nodes.loc[node, 'S'] = prop('S', 'H', out['H'], 'P', out['P'], fluid)
        nodes.loc[node, 'Q'] = prop('Q', 'H', out['H'], 'P', out['P'], fluid)
    elif 'P' in out and 'Q' in out:
        nodes.loc[node, 'Q'] = out['Q']
        nodes.loc[node, 'P'] = out['P']
        nodes.loc[node, 'T'] = prop('T', 'Q', out['Q'], 'P', out['P'], fluid)
        nodes.loc[node, 'S'] = prop('S', 'Q', out['Q'], 'P', out['P'], fluid)
        nodes.loc[node, 'H'] = prop('H', 'Q', out['Q'], 'P', out['P'], fluid)
    elif 'P' in out and 'S' in out:
        nodes.loc[node, 'S'] = out['S']
        nodes.loc[node, 'P'] = out['P']
        nodes.loc[node, 'T'] = prop('T', 'S', out['S'], 'P', out['P'], fluid)
        nodes.loc[node, 'Q'] = prop('Q', 'S', out['S'], 'P', out['P'], fluid)
        nodes.loc[node, 'H'] = prop('H', 'S', out['S'], 'P', out['P'], fluid)
    elif 'T' in out and 'S' in out:
        nodes.loc[node, 'S'] = out['S']
        nodes.loc[node, 'T'] = out['T']
        nodes.loc[node, 'P'] = prop('P', 'S', out['S'], 'T', out['T'], fluid)
        nodes.loc[node, 'Q'] = prop('Q', 'S', out['S'], 'T', out['T'], fluid)
        nodes.loc[node, 'H'] = prop('H', 'S', out['S'], 'T', out['T'], fluid)
    elif 'T' in out and 'Q' in out:
        nodes.loc[node, 'Q'] = out['Q']
        nodes.loc[node, 'T'] = out['T']
        nodes.loc[node, 'P'] = prop('P', 'Q', out['Q'], 'T', out['T'], fluid)
        nodes.loc[node, 'S'] = prop('Q', 'Q', out['Q'], 'T', out['T'], fluid)
        nodes.loc[node, 'H'] = prop('H', 'Q', out['Q'], 'T', out['T'], fluid)
    nodes.loc[node, 'fluid'] = fluid
    nodes.loc[node, 'G'] = G
    pass
def comp(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1, 'fluid']
    G = nodes.loc[node1, 'G']
    H2t = prop("H", "P", P2, "S", nodes.loc[node1, 'S'], fluid)
    Node(node2,G=G,P=P2,H=nodes.loc[node1, 'H']+(H2t-nodes.loc[node1, 'H'])/eff,fluid=fluid)
    blocks.loc[name, 'N'] = G * (nodes.loc[node2, 'H']-nodes.loc[node1, 'H'])
    pass
def heat_exch_2streams(name, node11, node12, node21, node22, **out):
    fluid1 = nodes.loc[node11]['fluid']
    fluid2 = nodes.loc[node21]['fluid']
    G1 = nodes.loc[node11]['G']
    G2 = nodes.loc[node21]['G']
    P11 = nodes.loc[node11]['P']
    P21 = nodes.loc[node21]['P']
    P12 = P11 - out['dP1'] if 'dP1' in out else P11
    P22 = P21 - out['dP2'] if 'dP2' in out else P21
    if 'T22' in out:
        Node(node22, G=G2, P=P22, T=out['T22'], fluid=fluid2)
        Q = G2 * (nodes.loc[node22,'H'] - nodes.loc[node21, 'H'])
        Node(node12, G=G1, P=P12, H=nodes.loc[node11, 'H']-Q/G1, fluid=fluid1)
    elif 'T12' in out:
        Node(node12, G=G1, P=P12, T=out['T12'], fluid=fluid1)
        Q = G1 * (nodes.loc[node11, 'H'] - nodes.loc[node12,'H'])
        Node(node22, G=G2, P=P22, H=nodes.loc[node21, 'H']+Q/G2, fluid=fluid2)
    elif 'Q22' in out:
        Node(node22, G=G2, P=P22, Q=out['Q22'], fluid=fluid2)
        Q = G2 * (nodes.loc[node22,'H'] - nodes.loc[node21, 'H'])
        Node(node12, G=G1, P=P12, H=nodes.loc[node11, 'H']-Q/G1, fluid=fluid1)
    elif 'Q12' in out:
        Node(node12, G=G1, P=P12, Q=out['Q12'], fluid=fluid1)
        Q = G1 * (nodes.loc[node11, 'H'] - nodes.loc[node12,'H'])
        Node(node22, G=G2, P=P22, H=nodes.loc[node21, 'H']+Q/G2, fluid=fluid2)
    step = 20
    T1 = [prop('T','P',P11-(P11-P12)/step*i,'H',nodes.loc[node11,'H']-Q/step*i/G1,fluid1) for i in range(step+1)]
    T2 = [prop('T','P',P21-(P21-P22)/step*i,'H',nodes.loc[node22,'H']-Q/step*i/G2,fluid2) for i in range(step+1)]
    dTmin = min([T1[i] - T2[i] for i in range(step+1)])
    blocks.loc[name] = [None,Q,dTmin,T1,T2]
    pass

def mix(name, node11, node12,node2):
    G11 = nodes.loc[node11]['G']
    G12 = nodes.loc[node12]['G']
    P11 = nodes.loc[node11]['P']
    P12 = nodes.loc[node12]['P']
    fluid11 = nodes.loc[node11]['fluid']
    fluid12 = nodes.loc[node12]['fluid']
    if P11 != P12: print("Давления входящих потоков в ",name," не равны!")
    if fluid11 != fluid12: print("Среды входящих потоков в ",name," отличаются!")
    H2 = (G11*nodes.loc[node11, 'H'] + G12*nodes.loc[node12, 'H'])/(G11+G12)
    Node(node2,G=G11+G12,P=P11,H=H2,fluid=fluid11)
    pass
def split(name, node1, node21,node22,massfrac):
    G1 = nodes.loc[node1]['G']
    nodes.loc[node21] = nodes.loc[node1]
    nodes.loc[node22] = nodes.loc[node1]
    nodes.loc[node21, 'G'] = G1*massfrac
    nodes.loc[node22, 'G'] = G1*(1-massfrac)
    pass

def turb(name, node1, node2, P2, eff):
    fluid = nodes.loc[node1, 'fluid']
    G = nodes.loc[node1, 'G']
    H2t = prop("H", "P", P2, "S", nodes.loc[node1, 'S'], fluid)
    Node(node2, G=G, P=P2, H=nodes.loc[node1, 'H']-(nodes.loc[node1, 'H']-H2t)*eff, fluid=fluid)
    blocks.loc[name, 'N'] = G * (nodes.loc[node1, 'H'] - nodes.loc[node2, 'H'])
    pass
def heat_exch(name, node1, node2,**out):
    P1 = nodes.loc[node1, 'P']
    fluid = nodes.loc[node1, 'fluid']
    G = nodes.loc[node1, 'G']
    P2 = P1 - out['dP'] if 'dP' in out else P1
    if 'Q' in out:
        Node(node2, G=G, P=P2, H=nodes.loc[node1, 'H']+out['Q']/G, fluid=fluid)
    elif 'T' in out:
        Node(node2, G=G, P=P2, T=out['T'], fluid=fluid)
    elif 'x' in out:
        Node(node2, G=G, P=P2, Q=out['x'], fluid=fluid)
    blocks.loc[name,'Q'] = abs(G*(nodes.loc[node2, 'H']-nodes.loc[node1, 'H']))
    pass
def comb_stoic(name, node11, node12,node2):
    H11 = nodes.loc[node11]['H']
    P11 = nodes.loc[node11]['P']
    F11 = nodes.loc[node11]['fluid']
    F12 = nodes.loc[node12]['fluid']
    H12 = nodes.loc[node12]['H']
    Gox = nodes.loc[node11]['G']
    Gf = nodes.loc[node12]['G']
    P2 = P11
    Qname = ['Methane_h', 'H2_h', 'CO_h', 'Methane_l', 'H2_l', 'CO_l']
    Qnum = [55515100, 141783257, 10103390, 50030044, 119957537, 10103390]
    Qc = dict(zip(Qname, Qnum))
    Mname = ['Methane', 'H2', 'CO', 'H2O', 'N2','O2','CO2', 'Ar', F11,F12]
    M = dict(zip(Mname, [prop('M', x) * 1000 for x in Mname]))
    f11num = re.sub('<[^>]+>', ' ', '<'+ F11.replace(']','<').replace('[','>')+'a>').split(' ')[1:-1]
    f11name = re.sub("\[[^]]*\]", '', F11.replace('REFPROP::','')).split('&')
    f11 = dict(zip(f11name, f11num))
    f12num = re.sub('<[^>]+>', ' ', '<' + F12.replace(']', '<').replace('[', '>') + 'a>').split(' ')[1:-1]
    f12name = re.sub("\[[^]]*\]", '', F12.replace('REFPROP::', '')).split('&')
    f12 = dict(zip(f12name, f12num))
    m11num = [prop('M', x)*1000*float(f11[x])/M[F11] for x in f11name]
    m11 = dict(zip(f11name,m11num))
    m12num = [prop('M', x)*1000*float(f12[x])/M[F12] for x in f12name]
    m12 = dict(zip(f12name,m12num))
    G_O2need = {}
    G_CO2frFuel = {}
    G_H2OfrFuel = {}
    G_CH4 = Gf*m12['Methane']
    G_O2need['Methane'] = G_CH4 * (2 * M['O2'] / M['Methane'])
    G_CO2frFuel['Methane'] = G_CH4 * (M['CO2'] / M['Methane'])
    G_H2OfrFuel['Methane'] = G_CH4 * (2*M['H2O'] / M['Methane'])
    G_H2 = Gf*m12['H2']
    G_O2need['H2'] = G_H2 * (0.5 * M['O2'] / M['H2'])
    G_H2OfrFuel['H2'] = G_H2 * (M['H2O'] / M['H2'])
    G_CO = Gf*m12['CO']
    G_O2need['CO2'] = G_CO * (0.5 * M['O2'] / M['CO'])
    G_CO2frFuel['CO2'] = G_CO * (M['CO2'] / M['CO'])
    G_O2in = m11['O2'] * Gox
    G_O2 = G_O2in - sum(G_O2need.values())
    G_CO2 = m11['CO2']*Gox + sum(G_CO2frFuel.values()) if 'CO2' in f11name else sum(G_CO2frFuel.values())
    G_H2O = m11['H2O']*Gox + sum(G_H2OfrFuel.values()) if 'H2O' in f11name else sum(G_H2OfrFuel.values())
    G_N2 = m11['N2']*Gox if 'N2' in f11name else 0
    G_Ar = m11['Ar']*Gox if 'Ar' in f11name else 0
    m_num = [x/(Gox + Gf) for x in [G_O2, G_CO2, G_H2O, G_N2, G_Ar]]
    m = dict(zip(['O2','CO2','H2O','N2','Ar'],m_num))
    mole_mix = sum([m[x]/M[x] for x in ['O2','CO2','H2O','N2','Ar']])
    w_num = [m[x]/M[x]/mole_mix for x in ['O2','CO2','H2O','N2','Ar']]
    w = dict(zip(['O2','CO2','H2O','N2','Ar'],w_num))
    fluid = "REFPROP::N2[" + str(w['N2']) + "]&CO2[" + str(w['CO2']) + "]&H2O[" + str(
        w['H2O']) + "]&O2[" + str(w['O2']) + "]&Ar[" + str(w['Ar'])+ "]"
    Qh = [m12[x] * Qc[x+'_h'] for x in f12name]
    Ql = [m12[x] * Qc[x + '_l'] for x in f12name]
    H2 = (Gox * H11 + Gf * (H12 + sum(Qh))) / (Gox + Gf)
    Node(node2,G=Gox + Gf,P=P2,H=H2,fluid=fluid)
    blocks.loc[name, 'Q'] = Gf * (sum(Ql))
    pass




# def heat(name, node11, node12, node21, node22, T12):
#     n = 20
#     fluid1 = nodes.loc[node11]['fluid']
#     fluid2 = nodes.loc[node21]['fluid']
#     H11 = nodes.loc[node11]['H']
#     T11 = nodes.loc[node11]['T']
#     P11 = nodes.loc[node11]['P']
#     P21 = nodes.loc[node21]['P']
#     T21 = nodes.loc[node21]['T']
#     G1 = nodes.loc[node11]['G']
#     H12 = prop("H", "T", T12, "P", P11, fluid1)
#     G2 = nodes.loc[node21]['G']
#     H21 = nodes.loc[node21]['H']
#     step = (H11 - H12) / (n-1)
#     t1 = np.zeros(n)
#     t2 = np.zeros(n)
#     Q = np.zeros(n)
#     h11 = H11
#     h21 = H21
#     t1[0] = T11
#     t2[-1] = T21
#     for i in range(n-1):
#         h12 = h11 - step
#         t1[i+1] = prop('T', 'H', h12, 'P', P11, fluid1)
#         Q[i+1] = Q[i] + G1*(h11-h12)
#         h11 = h12
#     for i in range(n-1):
#         h22 = h21 + (Q[-1-i] - Q[-2-i])/G2
#         t2[-2-i] = prop('T', 'H', h22, 'P', P21, fluid2)
#         h21 = h22
#     T22 = t2[0]
#     H22 = h22
#     S12 = prop('S', 'H', H12, 'P', P11, fluid1)
#     Q12 = prop('Q', 'H', H12, 'P', P11, fluid1)
#     S22 = prop('S', 'H', H22, 'P', P21, fluid2)
#     Q22 = prop('Q', 'H', H22, 'P', P21, fluid2)
#     nodes.loc[node12] = [T12, P11, H12, S12, Q12, G1, fluid1]
#     nodes.loc[node22] = [T22, P21, H22, S22, Q22, G2, fluid2]
#     blocks.loc[name, 'Q'] = Q
#     blocks.loc[name, 'T1'] = t1
#     blocks.loc[name, 'T2'] = t2
#     pass
#
# def cond(name, node1, node2):
#     P = nodes.loc[node1]['P']
#     H1 = nodes.loc[node1]['H']
#     fluid = nodes.loc[node1]['fluid']
#     G = nodes.loc[node1]['G']
#     T2 = prop('T', 'Q', 0, 'P', P, fluid)
#     H2 = prop('H', 'Q', 0, 'P', P, fluid)
#     S2 = prop('S', 'Q', 0, 'P', P, fluid)
#     nodes.loc[node2] = [T2, P, H2, S2, 0, G, fluid]
#     blocks.loc[name, 'Q'] = G*(H1 - H2)
#     pass