import Calc
import numpy as np
from time import *

n = int(input(">> "))
A = np.arange(n * n).reshape(n, n) + 1
for i in range(0, A.shape[0]): 
  for j in range(0, A.shape[1]): 
    if i < j: 
      A[j][i] = A[i][j]

A = A/10
#B = np.arange(n)
A2 = np.copy(A)

print("行列 A")
print(A)
#print("\n行列 B")
#print(B)



f = open('time.csv', 'a')
f.write("%s," % n)

start = time()    # 開始時刻

calc = Calc.Calc(A)
print("行列式\n",calc.determinant(A))
print("行列式\n",np.linalg.det(A))

#print("\n解")
#print(calc.cramer(B))

calc.eigenvalue()


goal = time()
t = goal - start   # 終了時刻

print("固有値")
print(calc.EVAL)
print("固有ベクトル")
print(calc.EVEC.T)

print(t)
f.write(str(t))
f.write(",")


print("----------------------------------")
start = time()    # 開始時刻

la, v = np.linalg.eig(A2)    # 行列Aの固有値・固有ベクトル
Q,R = np.linalg.qr(A2)
print("固有値")
print(la)
print("固有ベクトル：\n"+str(v))


goal = time()   # 終了時刻
t = goal - start
print(t)

f.write(str(t))
f.write(",\n")
f.close()
