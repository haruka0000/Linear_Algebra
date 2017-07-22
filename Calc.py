import copy
import numpy as np
import math

class Calc:
  def __init__(self, A):
    self.A = A
    self.A_0 = np.copy(A);
    self.n = A.shape[0]
    self.Q = np.zeros((self.n,self.n), dtype = float)
    self.R = np.zeros((self.n,self.n), dtype = float)
    self.EVAL = np.zeros(self.n, dtype = float)
    self.EVEC = np.zeros((self.n, self.n), dtype = float)


  def determinant(self, A):
    mat = np.copy(A)
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


  def trans_mat(self, mat_o): 
    mat = np.copy(mat_o)

    for i in range(0, mat.shape[0]):
      for j in range(0, mat.shape[1]):
        mat[j][i] = mat_o[i][j]
        
    return mat
 

  def qr(self):
    n = self.n
    X = np.zeros(n, dtype = float)

    for i in range(0, n, 1):
      
      X = self.A[:,i]

      for j in range(0, i, 1):
        t = 0.0
        for k in range(0, n, 1):
          t += self.A[k][i] * self.Q[k][j]
        self.R[j][i] = t
        self.R[i][j] = 0.0
        for k in range(0, n, 1):
          X[k] -= t * self.Q[k][j]

      
      s = 0.0
      
      for j in range(0, n, 1):
        s += X[j] * X[j]
      
      self.R[i][i] = math.sqrt(s)

      for j in range(0, n, 1):
        self.Q[j][i] = X[j] / self.R[i][i]




  def eigenvalue(self):
    n = self.n
    print(self.A)
    tmp = self.trans_mat(np.copy(self.A))
    print(tmp)
    flag = False
    for i in range(0,n):
      for j in range(0,n):
        if self.A[i][j] != tmp[i][j]:
          flag = True
          break

    if flag:
      
      if abs(self.determinant(self.A_0)):
        for i in range(1, 1000 ,1):
          self.qr()
          self.A = np.dot(self.R, self.Q)
          e = 0.0
          for j in range(1, n, 1):
            for k in range(0, j, 1):
              e += abs(self.A[j][k])
      
          if e < 0.00000000001:
            break

        for i in range(0, n, 1):
          self.EVAL[i] = self.A[i][i]
      
        self.eigenvector_qr()
    
      else:
        print("非正則かつ非対称行列であるため、計算できません(T_T)")
        exit()

    else:
      # 対称行列のときjacobi法
      self.jacobi()

  def eigenvector_qr(self):
    print("qr_vec")
    I = np.eye(self.n, dtype = float)

    for i in range(0, self.n, 1):
      U = np.zeros(self.n, dtype = float)
      U[0] = 1

      B = self.A_0 - self.EVAL[i] * I

      for k in range(200):
        U = np.dot(np.linalg.inv(B), U)
        U = U / np.linalg.norm(U)   # 発散を防ぐためにノルムで割る
        e = 0.0
        for j in range(0, self.n, 1):
          e += abs(U[j])
      
        if e < 0.00000000001:
          break
      self.EVEC[i] = U


  def mat_inv(self,mat):
    new_mat = np.eye(self.n, dtype = float)
    for i in range(0, self.n, 1):
      buf = 1 / mat[i][i];
      for j in range(0, self.n, 1):
        mat[i][j] *= buf
        new_mat[i][j] *= buf

      for j in range(0, self.n, 1):
        if(i != j):
          buf = mat[j][i]
          for k in range(0, self.n, 1):
            mat[j][k] -= mat[i][k] * buf
            new_mat[j][k] -= new_mat[i][k] * buf
    return new_mat


  def jacobi(self):
    print("jacobi")
    n = self.n
    A = np.copy(self.A_0)
    B = np.eye(n, dtype = float)
    EVEC = np.eye(n, dtype = float)

    for k in range(0, 200, 1):
      p = 0
      q = 0
      MAX,p,q = self.getMax(A)
      
      theta = 0.0
     
      
      if abs(A[p][p] - A[q][q]) < 0.00000000001:
        # self.A[p][p] == self.A[q][q]のとき回転角はπ/4
        theta = math.pi / 4.0
        if A[p][p] < 0:
          theta = -theta
      else:
        theta = 0.5 * math.atan(2.0 * A[p][q] / (A[p][p] - A[q][q]))
      
      co = math.cos(theta)
      si = math.sin(theta)
      
      B = np.eye(n, dtype = float)

      B[p][p] = co
      B[p][q] = si
      B[q][p] = -si
      B[q][q] = co

      A = np.dot(np.dot(B, A), B.T)
      EVEC = np.dot(EVEC, B.T)

      if MAX < 0.00000000001:
        break

    for i in range(n):
      self.EVAL[i] = A[i][i]
    self.EVEC = EVEC.T



  def getMax(self, mat):
    maxval = 0
    max_p = 0
    max_q = 0

    for i in range(self.n):
      for j in range(i+1,self.n):
        if maxval < abs(mat[i][j]) and i != j:
          maxval = abs(mat[i][j])
          max_p = i
          max_q = j

    if max_p > max_q:
      max_p, max_q = max_q, max_p
    return maxval, max_p, max_q


