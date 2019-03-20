#!/usr/bin/python3.6

from lib import Filter_SNP as FS
from lib import Find_primer as FP
from lib import Extent as EXT
import multiprocessing as mp
import pickle

vcf_file="MSG.vcf"


def Filter_SNP_fun(VCF_dir):
    print("_"*40+"SNP Filtering"+"_"*40)
    KY=FS.GT_SNP("KY")
    KY.read_vcf(VCF = VCF_dir, number_of_line_to_read = 10000000)
    KY.Filter_GT(dp_thred = 3 ,gq_thred= 88 )
    KY.Filter_ind(missing_thred = 0.5)
    Index=KY.Filter_SNP(missing_thred = 0.5, maf_thred = 0.01 ,index_only=False)
    return KY

def Find_primer_fun(KY):
    print("_"*40+"Find primer"+"_"*40)
    KY_FP = FP.Primer_SNP(KY)
    KY_FP.Get_refined_GT(maf_thred = 0.35, missing_thred = 0.05)
    KY_FP.Resolve_into_gr(kmer = 9, ld_thred = 0.9)
    KY_FP.Find_best_comb(core_num = 8)
    return KY_FP

def Retrieve_primer():
    with open('store.pckl', 'rb') as f:
        my_save = pickle.load(f)
    my_save[0].index = my_save[1]
    ARG = {'GT': my_save[0],'SNP_list': my_save[2],'Score': my_save[3]}
    return ARG

def Extent_primer(arg):
    print("_"*40+"Extent primer"+"_"*40)
    snp_primer = EXT.Primer(**arg)
    snp_primer.append(append_num = 10, core_num = 8)
    
if __name__ == '__main__':
    omit_preprocess = False
    if not omit_preprocess:
        KY = Filter_SNP_fun(vcf_file)
        KY_FP = Find_primer_fun(KY)
        To_save = [KY.GT,KY.GT.index,KY_FP.primer,KY_FP.score]
        with open('store.pckl', 'wb') as f:
            pickle.dump(To_save, f)

    arg = Retrieve_primer()
    Extent_primer(arg)

    #Primer = Extent_primer(KY_FP)




