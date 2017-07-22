import Calc
import numpy as np
from time import *

A = np.array([[8,5,5],[3,4,4],[1,5,3]], dtype = float)
B = np.array([[8],[3],[12]], dtype = float)
A2 = np.copy(A)

print("行列 A")
print(A)
print("\n行列 B")
print(B)




start = time()    # 開始時刻

calc = Calc.Calc(A)
#print("行列式\n",calc.determinant(A))

print("\n解")
print(calc.cramer(B))

goal = time()
t = goal - start   # 終了時刻
print("実行時間 ",t)

"""
print("----------------------------------")
start = time()    # 開始時刻

print("行列式\n",np.linalg.det(A))

goal = time()   # 終了時刻
t = goal - start
print("実行時間 ",t)
"""
