import pandas as pd

nodes = pd.DataFrame(index=['1', '2', '3', '4', '5', '6', '7', '8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23','24','25'],
                     columns=['T', 'P', 'H', 'S', 'Q', 'G', 'fluid'])
blocks = pd.DataFrame(index=['PUMP', 'REGEN', 'TURB', 'HEAT', 'COND','PP','EVAP','ECON','GPK'],
                      columns=['N', 'Q', 'T1', 'T2'])
nodes = pd.DataFrame(index=['1g', '2g', '3g', '4g', '5g', '6g', '7g', '8g','9g', '10g','1','2','3','4','5','6','7','8','9','10','11', '12', '13', '14', '15', '16', '17', '18', '19'],
                     columns=['T', 'P', 'H', 'S', 'Q', 'G', 'fluid'])