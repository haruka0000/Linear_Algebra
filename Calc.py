import copy
import numpy as np
import math

class Calc():

  def determinant(self, A):
  # 前進消去
    A = np.copy(A)
    n = A.shape[0]
    p = np.arange(n)	# [0,1,2,...,n-1]
    det = 1.0
    
    for k in range(n-1):
      
      # ピボット選択
      pivot_idx = p[k]
      pivot_max = 0
      for i in range(k, n):
        v = abs(A[p[i], k])
        if v > pivot_max:
          pivot_max = v
          pivot_idx = i
    
      # 実際には誤差があるので pivot_max == 0.0よりは
      # pivot_max < (小さい値) とした方が良いかも。
      if pivot_max == 0.0:
        return 0.0
      
      # ピボット行の交換
      if p[k] != pivot_idx:
        p[k], p[pivot_idx] = p[pivot_idx], p[k]
        det *= -1 # ピボット交換では符号を変える
     
      pivot = A[p[k],k]
      det *= pivot # 対角成分を掛け合わせる
      for i in range(k+1, n):
        l = A[p[i], k]/pivot
        for j in range(k+1, n):
          A[p[i], j] -= l * A[p[k], j]
    det *= A[p[n-1], n-1]	 # これを忘れずに
    return det

  
  def cramer(self, A, B):
    clone_A = copy.deepcopy(A)
    x = []
    i_A = 1/self.determinant(A)		# 1/|A|

    for i in range(0, A.shape[0]):
      clone_A = copy.deepcopy(A)
      for j in range(0, A.shape[0]):
        clone_A[j][i] = B[j]
      x.append(i_A* self.determinant(clone_A))
    
    return x


  def trans_mat(self, A): 
    clone_A = copy.deepcopy(A)
    
    for i in range(0, A.shape[0]):
      for j in range(0, A.shape[1]):
        clone_A[j][i] = A[i][j]
        
    return clone_A

  
  def qr(self, A):
    n = A.shape[0]
    Q = np.zeros((n,n), dtype = float)
    R = np.zeros((n,n), dtype = float)
    X = np.zeros(n, dtype = float)
    
    for i in range(0, n, 1):
      for j in range(0, n, 1):
        X[j] = A[j][i]

      for k in range(0, i, 1):
        t = 0.0
        for j in range(0, n, 1):
          t += A[j][i] * Q[j][k]
        R[k][i] = t
        R[i][k] = 0.0
        for j in range(0, n, 1):
          X[i] -= t * Q[j][k]
      s = 0.0
      for j in range(0, n, 1):
        s += X[j] * X[j]
      R[i][i] = math.sqrt(s)
      for j in range(0, n, 1):
        Q[j][i] = X[j] / R[i][i]

    #print(X)
    return Q,R

  def multiMat(self, ANS, A, B):
    n = A.shape[0]

    for i in range(0, n, 1):
      for j in range(0, n, 1):
        s = 0.0
        for k in range(0, n, 1):
          s += A[i][k] * A[k][j]
        ANS[i][j] = s

    return ANS


  def eigenvalue_qr(self, A):
    n = A.shape[0]

    for i in range(1, 1000 ,1):
      QR = self.qr(A)   #QR[0] = Q, QR[1] = R
      self.multiMat(A, QR[0], QR[1])
      print(A)
      print()

      e = 0.0
      for j in range(1, n, 1):
        for k in range(1, n, 1):
          e += abs(A[j][k])
      
      if(e < 0.00000000001):
        break


