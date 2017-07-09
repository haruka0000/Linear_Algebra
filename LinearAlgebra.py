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

A = np.array([[8,5,5],[3,4,4],[1,5,3]], dtype = float)
B = np.array([[8],[3],[12]], dtype = float)

#A = np.array([[1,0,0,4],[0,3,0,0],[0,0,5,0],[0,0,0,3]], dtype = float)
#B = np.array([[3],[0],[0],[0]], dtype = float)
#A = np.array([[3,-1,2],[-1,5,-3],[1,-1,3,]], dtype = float)
#B = np.array([[-7],[35],[-19]], dtype = float)

calc = Calc.Calc()
print("行列 A")
print(A)
print("\n行列 B")
print(B)
print("\n解")
print(calc.cramer(A, B))
print("\nQR分解")
qr = calc.qr(A)
print(qr[0])
print(qr[1])

print("固有値")
calc.eigenvalue_qr(A)
