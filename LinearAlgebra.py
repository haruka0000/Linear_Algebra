#####################################################
# 以下のプログラムを作成し，検証を行う．
# 1. 3次正方行列の行列式を求める．
# 2. 3次の連立一次方程式の解を求める．
# 3. 3次の正方行列の固有値，固有ベクトルを求める．
# 4. n次の場合の1. から 3. を求める．
#   (nは10以上，可能であればn = 100~1000を解いてみる)
#####################################################

import Calc
import numpy as np
import copy

A = np.array([[8,5,5],[3,4,4],[1,5,3]], dtype = float)
B = np.array([[8],[3],[12]], dtype = float)

#A = np.array([[1,0,0,4],[0,3,0,0],[0,0,5,0],[0,0,0,3]], dtype = float)
#B = np.array([[3],[0],[0],[0]], dtype = float)
#A = np.array([[3,-1,2],[-1,5,-3],[1,-1,3,]], dtype = float)
#B = np.array([[-7],[35],[-19]], dtype = float)


A2 = np.copy(A)

calc = Calc.Calc(A)
print("行列 A")
print(A)
print("\n行列 B")
print(B)
print("\n解")
print(calc.cramer(B))

#calc.qr()
calc.eigenvalue_qr()

print("\nQR分解")
print(calc.Q)
print(calc.R)
print("固有値")
print(calc.EVAL)
calc.eigenvector_qr()
print(calc.EVEC)

print("----------------------------------")
la, v = np.linalg.eig(A2)    # 行列Aの固有値・固有ベクトル
print("\nQR分解")
Q,R = np.linalg.qr(A2)
print(Q)
print(R)
print("固有値：\n"+str(la))
print("固有ベクトル：\n"+str(v))
