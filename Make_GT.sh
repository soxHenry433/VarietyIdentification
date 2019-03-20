#!/bin/bash

VCF="/mnt/e/MSG/afterNGS/Rawdata/TASSEL/BWAmem_TO1000_Q3.vcf"
GT="/mnt/e/ubuntu/SNPset/GT"

sudo chmod 777 $VCF

plink --vcf $VCF --make-bed --out $GT --allow-extra-chr
plink --bfile $GT --maf 0.01 --geno 0.5 --make-bed --out $GT"2" --allow-extra-chr
plink --bfile $GT"2" --indep-pairwise 50 1 0.9 --allow-extra-chr

awk 'BEGIN{
    getline SNP < "plink.prune.in"
}
/#/{
    print
}
!/^SC/{
    if($3 == SNP){
        print $0
        getline SNP < "plink.prune.in"
    }
}
' $VCF > MSG.vcf
wc MSG.vcf
rm GT*

