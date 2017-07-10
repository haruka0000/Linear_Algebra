import copy
import numpy as np
import math

class Calc:
  def __init__(self, A):
    self.A = A
    self.A_0 = A;
    self.n = A.shape[0]
    self.Q = np.zeros((self.n,self.n), dtype = float)
    self.R = np.zeros((self.n,self.n), dtype = float)

  def determinant(self, mat):
  # 前進消去
    n = self.n
    p = np.arange(n)	# [0,1,2,...,n-1]
    det = 1.0
    
    for k in range(n-1):
      
      # ピボット選択
      pivot_idx = p[k]
      pivot_max = 0
      for i in range(k, n):
        v = abs(mat[p[i], k])
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
     
      pivot = mat[p[k],k]
      det *= pivot # 対角成分を掛ける
      for i in range(k+1, n):
        l = mat[p[i], k]/pivot
        for j in range(k+1, n):
          mat[p[i], j] -= l * mat[p[k], j]
    det *= mat[p[n-1], n-1]
    return det

  
  def cramer(self,B):
    x = []
    
    clone_A = np.copy(self.A)
    i_A = 1/self.determinant(clone_A)		# 1/|A|
    n = self.n

    for i in range(0, n):
      clone_A = np.copy(self.A)
      for j in range(0, n):
        clone_A[j][i] = B[j]
      x.append(i_A * self.determinant(clone_A))
    
    return x


  def trans_mat(self): 
    clone_A = np.copy(self.A)
    
    for i in range(0, A.shape[0]):
      for j in range(0, A.shape[1]):
        clone_A[j][i] = A[i][j]
        
    return clone_A

  
  def qr(self):
    n = self.n
    Q = np.copy(self.Q)
    R = np.copy(self.R)
    X = np.zeros(n, dtype = float)
    
    for i in range(0, n, 1):
      for j in range(0, n, 1):
        X[j] = self.A[j][i]

      for k in range(0, i, 1):
        t = 0.0
        for j in range(0, n, 1):
          t += self.A[j][i] * self.Q[j][k]
        self.R[k][i] = t
        self.R[i][k] = 0.0
        for j in range(0, n, 1):
          X[j] -= t * self.Q[j][k]

      s = 0.0
      for j in range(0, n, 1):
        s += X[j] * X[j]
      self.R[i][i] = math.sqrt(s)
      for j in range(0, n, 1):
        self.Q[j][i] = X[j] / self.R[i][i]


  def multiMat(self, B, C):
    n = self.n
    ANS = np.zeros((n,n), dtype = float)

    for i in range(0, n, 1):
      for j in range(0, n, 1):
        s = 0.0
        for k in range(0, n, 1):
          s += B[i][k] * C[k][j]
        ANS[i][j] = s

    return ANS


  def eigenvalue_qr(self):
    n = self.n

    for i in range(1, 1000 ,1):
      self.qr()
      self.A = self.multiMat(self.R, self.Q)
      e = 0.0
      for j in range(1, n, 1):
        for k in range(0, j, 1):
          e += abs(self.A[j][k])
      
      if(e < 0.00000000001):
        break
    
  def eigenvector(self):
    


  # 対角要素を表示
  def disp_eigenvalue(self):
    n = self.n
    ev = np.zeros(n, dtype = float)
 
    for i in range(0, n, 1):
      ev[i] = self.A[i][i]
    
    return ev
