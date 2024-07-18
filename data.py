import pandas as pd

nodes = pd.DataFrame(index=['1','2','3','4','5','6','7','8'], columns=["T", "P", "H", "S", "Q", "G", "fluid"])
blocks = pd.DataFrame(index=["PUMP", "REGEN", "TURB", "HEAT", "COND"], columns=["N", "Q", "dT"])
