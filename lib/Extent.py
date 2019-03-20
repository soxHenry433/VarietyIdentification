#!/usr/bin/python3.6

import numpy as np
import multiprocessing as mp
from scipy.spatial.distance import pdist,squareform


class Primer(list):

    def __init__(self,GT,SNP_list,Score):
        super().__init__(SNP_list)
        self.score = Score
        self.GT = GT

    def cal_score(GT):
        ll = GT.index.tolist()
        GT = np.transpose(GT)
        primer = [ i in primer for i in ll]
        self.score = dis_fun(GT[:,primer])

    def append(self,append_num,core_num=7):
        SNP_num = len(self.GT.index)
        
        while append_num > 0:
            primer_loc = np.where(np.isin(self.GT.index,self))[0]
            primer_loc_mat = np.tile(primer_loc,(SNP_num,1))
            primer_loc_list = np.column_stack((primer_loc_mat,np.arange(SNP_num).reshape(SNP_num,1))).tolist()

            tGT = np.transpose(self.GT)
            size = SNP_num//core_num
            q = mp.Queue()
            #append_SNP(q,self.GT,self)
            for i in range(0, SNP_num, size):
                myprimer_mat = primer_loc_list[i:i+size]
                p = mp.Process(target=append_SNP,args=(q,tGT,myprimer_mat,1))
                p.start()

            p.join()
            L = []
            for i in range(core_num):
                L.append(q.get())

            Primer_list, Stat_list = zip(*L)
            max_eq = max(Stat_list,key=lambda item: item[1])
            best_primer = self.GT.index[Primer_list[Stat_list.index(max_eq)]]

            print(best_primer)
            print(max_eq)
            print(len(best_primer))

            if max_eq[1] <= self.score:
                print('Further extentention is futile')
                break
            self = Primer(GT = self.GT, SNP_list = best_primer, Score = max_eq[1])
            append_num -= 1

            

def append_SNP(Q,GT,primer_list,method_index=0):
    Score = []
    for i in primer_list:
        Score.append(dis_fun(GT[:,i]))
    max_eq = max(Score,key=lambda item: item[method_index])
    Append_result = (primer_list[Score.index(max_eq)],max_eq)
    Q.put(Append_result)

def dis_fun(GT_mat):
    pd = squareform(pdist(GT_mat,'chebyshev'))
    cc = np.sum(pd,axis=1)
    cc2 = (cc.mean(),(cc== (pd.shape[0]-1)).sum())
    return cc2

