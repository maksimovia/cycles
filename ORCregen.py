from CoolProp.CoolProp import PropsSI as prop
from data import nodes, blocks
from modules import comp, heat


# Исходные данные блоков
Ppump = 4e6
KPDpump = 0.8
KPDturb = 0.9
dTh = 10
dTr = 5

# Ввод исходных данных в (3)
X3 = 'R236ea'
G3 = 9999
P3 = 245000
T3 = prop("T", "P", P3, "Q", 0, X3)
H3 = prop("H", "P", P3, "Q", 0, X3)
S3 = prop("S", "P", P3, "Q", 0, X3)
Q3 = prop("Q", "P", P3, "Q", 0, X3)
nodes.loc['3'] = [T3, P3, H3, S3, Q3, G3, X3]

# Ввод исходных данных в (1)
X1 = 'REFPROP::N2[0.77]&O2[0.14]&H2O[0.06]&CO2[0.03]'
G1 = 500
T1 = 200+273.15
P1 = 1e5
H1 = prop("H", "P", P1, "T", T1, X1)
S1 = prop("S", "P", P1, "T", T1, X1)
Q1 = prop("Q", "P", P1, "T", T1, X1)
nodes.loc['1'] = [T1, P1, H1, S1, Q1, G1, X1]
Tout = 80+273.15

i=0

comp('PUMP', '3', '4', Ppump, KPDpump)
if i == 0:
    nodes.loc['5'] = nodes.loc['4']
else:
    print(0)
heat('HEATER', '1', '2', '5', '6', Tout, dTh)

print(nodes.loc[:,"T":"G"])
