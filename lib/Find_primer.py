#!/usr/bin/python3.6

import re
import numpy as np
from sklearn import cluster as cl
from itertools import product
from scipy.spatial.distance import pdist,squareform
#from sklearn.feature_selection import SelectPercentile
from random import choices
#from scipy.stats.stats import pearsonr   
import multiprocessing as mp
import objgraph

class Primer_SNP():
    
    def __init__(self,KY):
        self.FS = KY
        
    def Get_refined_GT(self,maf_thred = 0.35, missing_thred = 0.05):
        Index = self.FS.Filter_SNP(missing_thred,maf_thred,index_only=True)
        self.GT = self.FS.GT.slice_SNP(Index)
    
    def Resolve_into_gr(self,kmer=8,ld_thred = 0.8):
        km = cl.KMeans(n_clusters = kmer).fit(self.GT)
        gr = list()
        for i in range(kmer):
            gr.append(self.GT.index[km.labels_ == i])
        pruned_gr = list()
        for snp in gr:
            if(len(snp) <= 2):
                pruned_gr.append(snp)
                next
            ii = [i in snp for i in self.GT.index.tolist()]
            gt = self.GT.slice_SNP(ii)
            maskedarr = np.ma.array(gt, mask=(gt == 0.5))
            cov_arr = np.ma.corrcoef(maskedarr,rowvar=True,allow_masked=True)
            np.fill_diagonal(cov_arr,0)
            rm_index = np.where(cov_arr >= ld_thred)
            rm_arr = np.stack(rm_index)
            rm_choices = np.random.choice(2,rm_arr.shape[1])
            rm_indices = np.arange(rm_arr.shape[1])
            rm_snp = np.unique(rm_arr[rm_choices,rm_indices])
            if len(snp) == len(rm_snp):
                pruned_gr.append([snp[0]])
            else:
                pruned_gr.append(np.delete(snp,rm_snp))

        pruned_gr.sort(key = lambda x: len(x), reverse = True )
        p = 1
        for i in pruned_gr:
            print(i)
            if len(i) != 0:
                p = p*len(i)

        len_gr = [ len(i) for i in pruned_gr ]
        print(f'SNP num after ld trimming: {sum(len_gr)}')
        print(f'Total combinations: {p}')
        print(f'Combinations per loop: {p/len_gr[1]}')

        self.gr = pruned_gr

    def Find_best_comb(self,core_num = 8):
        
        for i in range(len(self.gr)):
            self.gr[i] = np.where(np.isin(self.GT.index,self.gr[i]))[0].tolist()

        gr_to_thread = self.gr.pop(0)
        size =len(gr_to_thread)//core_num
        chunks = [gr_to_thread[i:i+size] for i in range(0, len(gr_to_thread), size)]
        gr_to_loop = self.gr.pop(0)

        GT2 = np.transpose(self.GT)
        M = mp.Manager()
        q = M.Queue()
        cc = 0
        for i in gr_to_loop:
            tmp_gr = [[i]] + self.gr
            p = mp.Pool(processes = core_num)
            for j in range(len(chunks)):
                tmp_gr_thread = [chunks[j]] + tmp_gr 
                p.apply_async(func = par_cal_comb, args=(tmp_gr_thread,GT2,q))
            
            p.close()
            p.join()
            cc += 1
            print(f'Progress: {cc/len(gr_to_loop):.1%}')

        L = []
        for i in range(q.qsize()):
            L.append(q.get())
            
        Primer_list, Stat_list = zip(*L)
        max_mean = max(Stat_list,key=lambda item: item[0])
        best_primer = list(Primer_list[Stat_list.index(max_mean)])
        self.primer = self.GT.index[best_primer]
        self.score = max_mean
        print(self.primer)
        print(self.score)


def dis_fun(GT_mat):
    pd = squareform(pdist(GT_mat,'chebyshev'))
    cc = np.sum(pd,axis=1)
    cc2 = (cc.mean(),(cc== (pd.shape[0]-1)).sum())
    return cc2


class GT_table(np.ndarray):

    def __new__(cls,GT_mat):
        obj = np.asarray(GT_mat).view(cls)
        obj.index = np.arange(GT_mat.shape[0])
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.index = getattr(obj, 'index', None)

    def slice_SNP(self, index):
        obj = self[index,]
        obj.index = self.index[index]
        return obj


def par_cal_comb(tmp_gr_thread,GT2,q):
    append_SNP(Q = q,
        GT = GT2,
        primer_list = list(product(*tmp_gr_thread)),
        method_index = 0,
        )

def append_SNP(Q,GT,primer_list,method_index=0):
    Score = []
    for i in primer_list:
        Score.append(dis_fun(GT[:,i]))
    max_eq = max(Score,key=lambda item: item[method_index])
    Append_result = (primer_list[Score.index(max_eq)],max_eq)
    Q.put(Append_result)
