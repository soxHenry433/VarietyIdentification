#!/usr/bin/python3.6

import re
import numpy as np

class GT_SNP():
    
    def __init__(self, dataset_name="Default name"):
        self.name = dataset_name
        self.stat = set()
        print("Stat: ",self.stat)


    def read_vcf(self,VCF,number_of_line_to_read=10000000):

        my_re=re.compile(r"(?<=\t)([0-9\.]+/[0-9\.]+):[0-9,]+:(\d+):?(\d*)")
        with open(VCF) as vcf_open:
            self.GT = list()
            self.DP = list()
            self.GQ = list()
            self.SNP_info = list()

            cc = 0
            for line in vcf_open :

                if cc > number_of_line_to_read: 
                    break

                if re.match("#CHROM",line) is not None:
                    self.ind_name = re.split(r'\t',line)[9:]
                    self.ind_name[-1] = re.sub(r'\n', '', self.ind_name[-1])

                if line[0] != "C":
                    next
                else:
                    matched_grp = my_re.findall(line)
                    gt,dp,gq=zip(*matched_grp)
                    self.GT.append(gt)
                    self.DP.append(dp)
                    self.GQ.append(gq)

                    meta_data = re.split(r'\t',line)[0:5]
                    self.SNP_info.append(meta_data)

                    cc += 1

        self.SNP_info=np.array(self.SNP_info)

        self.GT=np.array(self.GT)
        for pat,sub in [("1/1","1"),("0/0","0")]:
            np.place(self.GT,self.GT==pat,sub)
        np.place(self.GT,np.isin(self.GT,["1","0"],invert=True),"0.5")
        self.GT=self.GT.astype(float)
        self.GT = GT_table(self.GT)

        self.DP=np.array(self.DP)
        self.DP=self.DP.astype(int)

        self.GQ=np.array(self.GQ)
        np.place(self.GQ,self.GQ=='',"-9")
        self.GQ=self.GQ.astype(int)

        self.stat.add("VCF read")
        snp_num, ind_num = self.GT.shape
        print("Datasize:\n\tSNP: ",snp_num,"\n\tIndividual: ",ind_num)
        print("Stat: ",self.stat)
        

    def Filter_GT(self,dp_thred = 3 ,gq_thred= 88 ):
        if len(self.stat) == 0:
            raise ImportError("No VCF imported yet")
        snp_num,ind_num = self.GT.shape
        filter_arr = np.logical_and(self.DP >= dp_thred,self.GQ >= gq_thred)
        num_missings1 = (self.GT == 0.5).sum()/(snp_num*ind_num)
        print(f'Dataset missing rate: {num_missings1:.2%}')

        self.GT = np.where(filter_arr, self.GT, 0.5)
        self.GT = GT_table(self.GT)
 
        num_missings2 = (self.GT == 0.5).sum()/(snp_num*ind_num)
        print(f'Dataset missing rate after filering: {num_missings2:.2%}')
        self.stat.add("Data_point Filtered") 
        print("Stat: ",self.stat)


    def Filter_SNP(self,missing_thred=0.1, maf_thred=0.4,index_only=True):
        if len(self.stat) == 0:
            raise ImportError("No VCF imported yet")

        print(f'SNP number before filter: {self.GT.shape[0]}')
        print(f'Missing threshold: {missing_thred}\nMAF threshold: {maf_thred}')

        filter_index = np.apply_along_axis(filter_snp_fun,1,self.GT,missing_thred,maf_thred)
        print(f'SNP number after filter: {np.sum(filter_index)}')

        #filter_index = np.where(filter_index)[0]

        if not index_only:
            self.GT = self.GT.slice_SNP(filter_index)
            self.stat.add("SNP Filtered") 
            print("Stat: ",self.stat)

        return filter_index


    def Filter_ind(self,missing_thred=0.75):
        if len(self.stat) == 0:
            raise ImportError("No VCF imported yet")
        maf_thred=0
        print(f'Ind number before filter: {self.GT.shape[1]}')
        print(f'Missing threshold: {missing_thred}\nMAF threshold: {maf_thred}')

        filter_index = np.apply_along_axis(filter_snp_fun,0,self.GT,missing_thred,maf_thred)
        print(f'Ind number after filter: {np.sum(filter_index)}')

        self.GT = self.GT[:,filter_index]
        print(f'{np.array(self.ind_name)[~filter_index]} is/are removed')

        self.stat.add("Ind Filtered") 
        print("Stat: ",self.stat)



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


def filter_snp_fun(gt,missing_thred=0.1,maf_thred=0.01):
    gt_uniq = dict.fromkeys([0.5,0,1], 0)

    for i in gt:
        gt_uniq[i] += 1

    missing_rate = gt_uniq[0.5]/len(gt)
    maf = min([gt_uniq[1],gt_uniq[0]])/max(1,(len(gt)-gt_uniq[0.5]))

    if missing_rate < missing_thred and maf >= maf_thred:
        return(True)
    else:
        return(False)


#看不懂如何做subclass of np.ndarray，抄下來ㄉ，求講解QQ
# -> cls 的作用為何?
# -> __array_finalize__ 官網看不懂
#要如何不要那麼多self (例如self.GT,self.DP...)，覺得python class內的scope有點小麻煩
#他好像有點慢...QQ自定義的GT_table會嚴重影響效能嗎?






