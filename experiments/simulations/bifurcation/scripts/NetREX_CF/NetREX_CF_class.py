# There is a bug in the original jupyternotebook. In particular, X's and Y's definitions are exchanged.
# I fixed it below

import numpy as np
import pandas as pd
import sys
from implicit import als
import math
import cvxopt
from scipy.linalg import fractional_matrix_power, block_diag
from pylab import *
import argparse
#import mosek
from scipy.sparse import csr_matrix
import copy
import cvxpy as cvx
import anndata
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import os 
import progressbar

os.environ["MKL_NUM_THREADS"] = '1'


def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def MatKeepRowTop(Mat, TopNum = 3):
    MatOut = np.zeros(Mat.shape)
    NumRow, NumCol = Mat.shape
    for i in range(NumRow):
        tmprow = Mat[i,:]
        topid = np.argsort(-tmprow)
        Num0 = np.count_nonzero(tmprow)
        Num = np.min([TopNum, Num0])
        MatOut[i,topid[:Num]] = 1.
    return MatOut

def MatKeepRowTopVal(Mat, TopNum = 3):
    MatOut = np.zeros(Mat.shape)
    NumRow, NumCol = Mat.shape
    for i in range(NumRow):
        tmprow = Mat[i,:]
        topid = np.argsort(-tmprow)
        Num0 = np.count_nonzero(tmprow)
        Num = np.min([TopNum, Num0])
        MatOut[i,topid[:Num]] = tmprow[topid[:Num]]
    return MatOut

def MatKeepRowTopAbsVal(Mat, TopNum = 3):
    MatOut = np.zeros(Mat.shape)
    NumRow, NumCol = Mat.shape
    for i in range(NumRow):
        tmprow = Mat[i,:]
        topid = np.argsort(-np.abs(tmprow))
        Num0 = np.count_nonzero(tmprow)
        Num = np.min([TopNum, Num0])
        MatOut[i,topid[:Num]] = tmprow[topid[:Num]]
    return MatOut

def MatKeepColTop(Mat, TopNum = 3):
    MatOut = np.zeros(Mat.shape)
    NumRow, NumCol = Mat.shape
    for i in range(NumCol):
        tmprow = Mat[:,i]
        topid = np.argsort(-tmprow)
        Num0 = np.count_nonzero(tmprow)
        Num = np.min([TopNum, Num0])
        MatOut[topid[:Num],i] = 1.
    return MatOut

def MatKeepColTopVal(Mat, TopNum = 3):
    MatOut = np.zeros(Mat.shape)
    NumRow, NumCol = Mat.shape
    for i in range(NumCol):
        tmprow = Mat[:,i]
        topid = np.argsort(-tmprow)
        Num0 = np.count_nonzero(tmprow)
        Num = np.min([TopNum, Num0])
        MatOut[topid[:Num],i] = 1.
    MatOut = MatOut*Mat
    return MatOut
    

class CFNetREX: 
    #init
    def __init__(self, E, C, lambda_, alpha, beta, reg=1., lambda_a = 1., lambda_s = 1., samplesize=100, eta=.1, iternum = 100, p = 100):
        #parameters
        #self.theta = theta # scalar > 0
        self.alpha = alpha
        self.xbound = 1.0
        self.ybound = 1.0
        self.iternum = iternum
        self.reg = reg
        self.lambda_ = lambda_ #parameter matrix
        self.lambda_a = lambda_a
        self.lambda_s = lambda_s
        self.Aregular  = 1.
        self.samplesize = samplesize
        #read in exp, read prior data 
        #self.ReadExp(ExpFile) #get ExpMat from this
        self.ExpMat = E #expression for now
        self.C = MatKeepRowTopVal(C,200)#matrix of confidences (weights)
        self.P = (self.C != 0).astype(float) 
        self.NumGene = E.shape[0]
        self.NumTF = C.shape[1]
        self.NumExp = E.shape[1]
        self.BestRankMat = np.zeros(C.shape)
        #parameter matrix eta
        #self.eta = self.lambda_*self.C + theta
        self.i_ = range(self.NumGene) #range of i
        self.j_ = range(self.NumTF) #range of j 
        #controls amount of factors in xi/yj
        self.p = p
        self.beta = beta
        self.eta = self.P*eta + (1-self.P)*(self.lambda_*self.alpha*(1.+self.alpha/self.beta) + eta)
        
        #set solver
        #self.solver = solver
        
    #Put read in exp here later
    
    #initialization
    def initialization(self):
        #!!! other parameters!!!
        #Use implicit package for CF problems to initialize X & Y: 
        implicit_model = als.AlternatingLeastSquares(factors=self.p, regularization=.1, iterations=100)
        #add to get new weight matrix
        C_adj = self.C#1 + 40*self.C
    
        
        #train model on weight matrix
        #implicit_model.fit(csr_matrix(C_cutoff))
        print("Initialization begin...")
        implicit_model.fit(50*csr_matrix(MatKeepRowTopVal(C_adj,100)))#
        
        
        #initial X (matrix of xi vectors)
        self.X = implicit_model.user_factors.T
        XX = copy.deepcopy(implicit_model.item_factors.T)
        XColNorm = np.linalg.norm(self.X, axis=0)
        self.X = self.X #/ XColNorm
        
        #initial Y (matrix of yj vectors)
        self.Y = implicit_model.item_factors.T
        YY = copy.deepcopy(implicit_model.user_factors.T)
        YColNorm = np.linalg.norm(self.Y, axis=0)
        self.Y = self.Y #/ YColNorm
        
        self.C = self.C + 1.
        #self.eta = self.lambda_*(np.power(self.alpha+self.C,2)/self.alpha - 1 - self.C)#self.lambda_*self.C + self.theta
    
        #try using NCA to get initial A
        self.NCA_()
        print("Initialization done...")
        
    
    #try to program iterative NCA
    def NCA_(self):
        #set initial S to C
        RankMat = self.X.T.dot(self.Y)
        RankMat_RT20 = MatKeepRowTop(RankMat, 100) + MatKeepColTop(RankMat, 500)#np.where(RankMat>0.5, 1.0, 0.)#MatKeepRowTop(RankMat, 100) + MatKeepColTop(RankMat, 100)
        #self.S = RankMat_RT20 #copy.deepcopy(self.C)
        #binarize S to use for finding A 
        S_init = (RankMat_RT20 > 0).astype(float) #* (2*np.random.rand(self.NumGene, self.NumTF)-1.)*.1#(self.S > 5.).astype(int).astype(float)
        for k in range(20):#progressbar.progressbar(range(20)):
            #print(k)
            A_init = np.linalg.inv(S_init.T.dot(S_init)+ 0.1*np.eye(self.NumTF)).dot(S_init.T).dot(self.ExpMat)
            # update S
            for i in range(self.NumGene):
                Indn0 = np.where(S_init[i,:]!=0)[0]
                Et = self.ExpMat[i,:]
                At = A_init[Indn0,:]
                St = Et.dot(At.T).dot(np.linalg.inv(At.dot(At.T)+0.1*np.eye(At.shape[0])))
                S_init[i,Indn0] = St 
            
        self.A = A_init
        self.S = S_init
        
    #try to program iterative NCA
    def NCA_Input(self, SMat):
        #binarize S to use for finding A 
        S_init = np.where(SMat!=0, 1.0, 0)#(RankMat_RT20 > 0).astype(float) #* (2*np.random.rand(self.NumGene, self.NumTF)-1.)*.1#(self.S > 5.).astype(int).astype(float)
        for k in progressbar.progressbar(range(20)):
            #print(k)
            A_init = np.linalg.inv(S_init.T.dot(S_init)+ 0.1*np.eye(self.NumTF)).dot(S_init.T).dot(self.ExpMat)
            # update S
            for i in range(self.NumGene):
                Indn0 = np.where(S_init[i,:]!=0)[0]
                Et = self.ExpMat[i,:]
                At = A_init[Indn0,:]
                St = Et.dot(At.T).dot(np.linalg.inv(At.dot(At.T)+0.1*np.eye(At.shape[0])))
                S_init[i,Indn0] = St 
            
        self.A = A_init
        self.S = S_init
        
    #Define the objective function 
    def ObjFunction(self):
        S_0 = (self.S != 0).astype(float) 
        #print("nonzeros in S_0 is %d" % np.count_nonzero(S_0))
        self.a_ = self.lambda_*(self.beta*S_0 + self.alpha)#self.C
        W = np.sqrt(self.a_) 
        #W2 = np.sqrt(self.lambda2_*self.C)
        Val = np.linalg.norm(self.ExpMat - self.S.dot(self.A), 'fro')**2 +np.linalg.norm(W*(S_0 + (1-S_0)*self.P - self.X.T.dot(self.Y)),
                    'fro')**2 + np.sum(self.eta*S_0) + self.lambda_a*np.linalg.norm(self.A, 'fro')**2+self.lambda_s*np.linalg.norm(self.S, 'fro')**2 + self.reg*np.linalg.norm(self.X,'fro')**2 + self.reg*np.linalg.norm(self.Y,'fro')**2#np.linalg.norm(np.sqrt(self.eta*S_0), 'fro')**2
        return Val
    
    def PALM_A(self):
        mu_A = np.linalg.norm(self.S.T.dot(self.S)+self.lambda_a*np.eye(self.NumTF), 'fro')
        #solution for A
        Tmp = self.S.T.dot(self.S)+ (self.lambda_a+(mu_A/2.0))*np.identity(self.NumTF)
        self.A = np.linalg.pinv(Tmp).dot(( self.ExpMat.T.dot(self.S) + (mu_A/2.0) *self.A.T).T) 
        
    def PALM_S(self):
        mu_S = np.linalg.norm(self.A.dot(self.A.T)+self.lambda_s*np.eye(self.NumTF), 'fro')
        Uk = self.S - (2./mu_S)*(self.S.dot(self.A).dot(self.A.T)+ self.lambda_s*self.S - self.ExpMat.dot(self.A.T))
        #ct = ((2./mu_S)*(self.lambda_*((self.C+self.alpha)*np.power(1-self.X.T.dot(self.Y),2))+self.eta))
        #indt = np.where(ct<=0)
        #print(indt)
        #c = np.sqrt((2./mu_S)*(self.lambda_*((self.beta*(1.-2*self.X.T.dot(self.Y)))+(self.C+self.alpha)*np.power(1-self.X.T.dot(self.Y),2))+self.eta))
        R = self.X.T.dot(self.Y)
        ctmp = self.lambda_*((self.alpha+self.beta)*(1.-self.P)*(1+2.*(self.P-R)) + self.beta*np.power((self.P-R),2)) + self.eta
        c = np.sqrt((2./mu_S)*ctmp)
        #self.C   
        
        #print()

        #solution for S   
        #Sequal = np.where(np.absolute(Uk) == c, c*np.random.randint(2, size=(self.NumGene,self.NumTF)), 0)
        Sgreater = np.where(np.absolute(Uk) > c, Uk, 0)        
        self.S = Sgreater#Sequal + Sgreater

        
    def CG(self, A, b, x, imax=800, epsilon = 0.000001):

        i = 0
        r = b - A @ x
        d = r.copy()
        delta_new = r.T @ r
        delta_0 = delta_new
        #print(delta_new)
        while i < imax and delta_new > epsilon**2 * delta_0:
            q = A @ d
            alpha = float(delta_new / (d.T @ q))
            x = x + alpha * d
            if i%50 == 0:
                r = b - A @ x
            else:
                r = r - alpha * q
            delta_old = delta_new
            delta_new = r.T @ r
            beta = float(delta_new / delta_old)
            d = r + beta * d
            i = i + 1
            #print(0.5* x.T@A@x + b.T@x)
            #print(i)
            #steps = np.append(steps, np.asarray(x), axis=1)
        return x.T

            
           
    def PALM_x(self):
        #define variables for matrix formulation
        S_0 = (self.S != 0).astype(float)
        a_ = self.lambda_*(self.beta*S_0 + self.alpha)#self.C
        #Valnew = self.ObjFunction()
        #print("BeforeX: *Obj %f" % (Valnew))
        RandomInd = np.arange(self.NumGene)
        np.random.shuffle(RandomInd)
        if self.samplesize < self.NumGene:
            SRandomInd = RandomInd[0:self.samplesize]
            #print(SRandomInd)
        else:
            SRandomInd = RandomInd 
        for i in SRandomInd:#progressbar.progressbar(SRandomInd): #self.i_: 
            #i = RandomInd[ii]
            x = self.X[:,i]
            x = x[:, np.newaxis]
            A_i = np.diag([a_[i,j] for j in self.j_])
            A_i_rowvec = np.array([a_[i,j] for j in self.j_])
            #C_i_rowvec = np.array([self.lambda2_*self.C[i,j] for j in self.j_])
            S_i = np.array([S_0[i,j] + (1.-S_0[i,j])*self.P[i,j] for j in self.j_]) #
            S_i = S_i[newaxis,:]
            #P_i = np.array([self.P[i,j] for j in self.j_]) 
            #P_i = P_i[newaxis,:]
            mu_x = np.linalg.norm(2*((self.Y*A_i_rowvec).dot(self.Y.T)+self.reg*np.identity(self.p)), 'fro')                            
            phi = ((mu_x/2)+self.reg)*np.identity(self.p) + (self.Y*A_i_rowvec).dot(self.Y.T) #+ 0.01*np.identity(self.p) 
            varphi = 2*(S_i*A_i_rowvec).dot(self.Y.T) + mu_x*x.T 
            #phi = (phi + phi.T)/2
            
            xt = np.random.rand(self.p,1)
            varphi = varphi.T#[:, np.newaxis]
            self.X[:,i] = self.CG(2*phi, varphi, xt)#self.solve_qcqp(phi, varphi)
        
    def PALM_y(self):
        #matrix formulation
        S_0 = (self.S != 0).astype(float)
        a_ = self.lambda_*(self.beta*S_0 + self.alpha)#self.C
        #Valnew = self.ObjFunction()
        #print("BeforeY: *Obj %f" % (Valnew))
        #sys.stdout.flush()
        RandomInd = np.arange(self.NumTF)
        np.random.shuffle(RandomInd)
        if self.samplesize < self.NumTF:
            SRandomInd = RandomInd[0:self.samplesize]
            #print(SRandomInd)
        else:
            SRandomInd = RandomInd     
        for j in SRandomInd:#progressbar.progressbar(SRandomInd): #self.j_:
            y = self.Y[:,j]
            y = y[:, np.newaxis]
            A_j = np.diag([a_[i,j] for i in self.i_])
            A_j_rowvec = np.array([a_[i,j] for i in self.i_])
            #C_j_rowvec = np.array([self.lambda2_*self.C[i,j] for i in self.i_])
            S_j = np.array([S_0[i,j]+(1.-S_0[i,j])*self.P[i,j] for i in self.i_])#S_0[i,j]
            S_j = S_j[:, np.newaxis]
            #P_j = np.array([self.P[i,j] for i in self.i_])
            #P_j = P_j[:, np.newaxis]
            mu_y = np.linalg.norm(2*((self.X*A_j_rowvec).dot(self.X.T)+self.reg*np.identity(self.p)), 'fro')
            phi = ((mu_y/2)+self.reg)*np.identity(self.p) + (self.X*A_j_rowvec).dot(self.X.T) 
            varphi = 2*(S_j.T*A_j_rowvec).dot(self.X.T) + mu_y*y.T

            yt = np.random.rand(self.p,1)
            varphi = varphi.T#[:, np.newaxis]
            self.Y[:,j] = self.CG(2*phi, varphi, yt)#self.solve_qcqp(phi, varphi)
            #print("Norm Y is %f" %(np.linalg.norm(self.Y[:,j])))
            #print(self.Y[:,j])
            #Valnew = self.ObjFunction()
            #print("Y >    : *Obj %f" % (Valnew))
            #sys.stdout.flush()
            if np.isnan(np.sum(self.Y[:,j])):
                print("nan values appears at ")
                print(j)
                print(self.Y[:,j])
                sys.stdout.flush()
                

            
    #algorithm
    def NetREX_Plus_algorithm(self):
        self.initialization()
        conv_point = 1 #for now?  

        Valold = self.ObjFunction()
        print("Itr   0: ObjVal: %f" % (Valold))
        Indt = 0
        for k in range(self.iternum): 
            
            self.PALM_A()
            Valnew = self.ObjFunction()
            #print("After A: *Obj %f" % (Valnew), end= ' ')
            sys.stdout.flush()
            
            
            self.PALM_S()
            Valnew = self.ObjFunction()
            #print("After S: *Obj %f" % (Valnew), end= ' ')
            sys.stdout.flush()
            
            
            self.PALM_x()
            Valnew = self.ObjFunction()
            #print("After X: *Obj %f" % (Valnew), end= ' ')
            sys.stdout.flush()
            
            
            self.PALM_y()
            Valnew = self.ObjFunction()
            #print("After Y: *Obj %f" % (Valnew))
            sys.stdout.flush()
            
            RankMat = self.X.T.dot(self.Y)

            
        
            Valnew = self.ObjFunction()
            print("Itr   %d: ObjVal: %f" % (k+1, Valnew))
            sys.stdout.flush()
            check_pt = np.abs(Valold-Valnew)
            if (check_pt < conv_point) :
                print("Converge!")
                sys.stdout.flush()
                break
            Valold = Valnew
        return Valold
    
    def ConfEdge(self):
        ConfMean = np.zeros(self.S.shape)
        EdgeInd = np.nonzero(self.S)
        NumN0 = np.count_nonzero(self.S)
        ConfGroup = np.zeros((20,NumN0))
        #print(RankGroup)
        for ii in range(20):
            print(ii)
            RandSampleId = np.random.randint(self.NumExp,size=(self.NumExp))
            RandSampleId = RandSampleId[0:int(self.NumExp/2.)]
            #print(len(set(RandSampleId)))
            ExpSampled = self.ExpMat[:,RandSampleId]
            TFASampled = self.A[:,RandSampleId]
            Conf = np.zeros(NumN0)
            for i in range(NumN0):
                IndGene = EdgeInd[0][i]
                IndTF = EdgeInd[1][i]
                Par = copy.deepcopy(self.S[IndGene,:])
                #print(Par)
                Sigma_full = np.linalg.norm(ExpSampled[IndGene,:] - Par.dot(TFASampled))**2
                Par[IndTF] = 0.
                #print(Par)
                Sigma_miss = np.linalg.norm(ExpSampled[IndGene,:] - Par.dot(TFASampled))**2
                Conf[i] = 1. - (Sigma_full / Sigma_miss)
                if Conf[i] < 0:
                    #print("Conf[i]<0, something is wrong!!!")
                    Conf[i] = 0.
            #print(Conf)
            RankId = np.argsort(Conf)
            #ConfRank = np.zeros(NumN0)
            #ConfRank[RankId] = np.arange(NumN0) + 1.
            #print(ConfRank)
            ConfGroup[ii,:] = Conf
        #print(RankGroup)
        RankGroupMean = np.mean(ConfGroup, axis=0)
        ConfMean[EdgeInd[0], EdgeInd[1]] = RankGroupMean 
        self.ConfEdgeMat = ConfMean
    
        return ConfMean
        

class CFNetREX_Consensus: 
    def __init__(self, E, C, lambda_, alpha, beta, reg=1., lambda_a = 1., lambda_s = 1., samplesize=50, eta=.1, iternum = 100, p = 100, runs = 10):
        self.E = E
        self.C = C
        self.lambda_ = lambda_
        self.alpha = alpha
        self.beta = beta
        self.reg = reg
        self.lambda_a = lambda_a
        self.lambda_s = lambda_s
        self.samplesize = samplesize
        self.eta = eta
        self.iternum = iternum
        self.p = p
        self.runs = runs
    
    def ConsensusRun(self):
        Rind = np.random.permutation(self.E.shape[1])
        E_Input = self.E#[:,Rind[0:int(self.E.shape[1]*0.2)]]
        self.CFNetNEX_RUNs = []
        for i in range(self.runs):
            self.CFNetNEX_RUNs.append(CFNetREX(E_Input, self.C, self.lambda_, self.alpha, self.beta, self.reg, self.lambda_a, self.lambda_s, self.samplesize, self.eta, self.iternum, self.p))
        
        for i in range(self.runs):
            self.CFNetNEX_RUNs[i].NetREX_Plus_algorithm()
            self.CFNetNEX_RUNs[i].ConfEdge()
