#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 23:32:22 2020

@author: panzhang
"""

import os
import pandas as pd
import numpy as np
import re

            
def filter_name(df, ENS_id):
    colName= df.columns.tolist()
    data_id = df.index.tolist()
    df['id']=df.index
    df['new_id']=df['id'].apply(lambda x: x.split('_')[0])
    df.index = df['new_id']
    new_df = df.loc[ENS_id]
    new_df.index = new_df["id"]
    new_df=new_df[colName]
    return (new_df)       

os.chdir('/Users/panzhang/Desktop/AS/TWN')
geneAnnot = pd.read_table("gene_annot_protein_coding.txt", sep="\t",index_col=0)
geneAnnot["gene_id"] = geneAnnot.index
geneAnnot.drop_duplicates("gene_id", False, inplace=True )
geneAnnot = geneAnnot.drop(columns="gene_id")
geneAnnot.to_csv("geneAnnot.csv")


transcriptAnnot = pd.read_table("transcript_annot_protein_coding.txt", sep="\t",index_col=0)
transcriptAnnot["transcript_id"] = transcriptAnnot.index
transcriptAnnot.index = transcriptAnnot["gene_id"]
geneIds = geneAnnot.index.tolist()
transcriptAnnot_new = transcriptAnnot.loc[geneIds]
transcriptAnnot =  transcriptAnnot_new 
transcriptAnnot.index = transcriptAnnot["transcript_id"] 
transcriptAnnot= transcriptAnnot.drop(columns="transcript_id")
transcriptAnnot.to_csv("transcriptAnnot.csv")


geneCount = pd.read_table("OIH_NTG_gene_count.txt", sep=" ",index_col=0)
ENS_id = geneAnnot.ensembl_gene_id.tolist()
geneCount = filter_name(geneCount, ENS_id)


geneTPM = pd.read_table("OIH_NTG_gene_TPM.txt", sep=" ",index_col=0)
geneTPM.columns = geneCount.columns
geneTPM = filter_name(geneTPM, ENS_id)

            

isoCount = pd.read_table("OIH_NTG_isoforms_count.txt", sep=" ",index_col=1)
isoCount = filter_name(isoCount, ENS_id)
isoCount.index = isoCount["transcript_id"]
isoCount = isoCount.drop(columns="transcript_id")
isoCount.columns = geneCount.columns

isoTPM = pd.read_table("OIH_NTG_isoforms_TPM.txt", sep=" ",index_col=1)
isoTPM = filter_name(isoTPM, ENS_id)
isoTPM.index = isoTPM["transcript_id"]
isoTPM = isoTPM.drop(columns="transcript_id")
isoTPM.columns = geneCount.columns


#expressed in 75% of samples   
def filter_df(df,i):
    df = df[df.gtSum> i]
    df = df[df.gcSum> i]
    df= df.drop(columns = ["gtSum","gcSum"])
    #df = np.log(df+1)
#    df_var.columns = ["dfVar"]
#    df = pd.concat([df, df_var],axis=1)
#    df = df.sort_values(by="dfVar", ascending=False)
#    df= df.drop(columns = ["dfVar"])
    return(df)

#filter by count > 6 and TPM > 6 in at least 15 samples    
def filter_gene(geneCount, geneTPM, t= "gene"):
    gcSum=(geneCount > 6).sum(axis=1).to_frame()
    gtSum=(geneTPM > 1).sum(axis=1).to_frame()
    gtSum.columns = ["gtSum"]
    gcSum.columns = ["gcSum"]
    geneCount = pd.concat([geneCount, pd.concat([gcSum,gtSum],axis=1)],axis=1)
    geneTPM = pd.concat([geneTPM, pd.concat([gcSum,gtSum],axis=1)],axis=1)   
    geneCount = filter_df(geneCount, 15)
    geneTPM = filter_df(geneTPM, 15)
    if t == "gene": 
        return (geneCount, geneTPM)
    elif t== "iso":
        return geneCount.index.tolist()
        

NA_geneCount = geneCount[['OIH_M1_NA', 'OIH_M2_NA', 'OIH_M3_NA', 'OIH_M4_NA', 'OIH_M5_NA',
       'CON_V1_NA', 'CON_V2_NA', 'CON_V3_NA', 'CON_V4_NA', 'CON_V5_NA',
       'N1_NA', 'N2_NA', 'N3_NA', 'N4_NA', 'N5_NA', 
       'V6_NA', 'V7_NA', 'V8_NA', 'V9_NA', 'V10_NA']]
TG_geneCount = geneCount[['OIH_M1_TG', 'OIH_M2_TG', 'OIH_M3_TG', 'OIH_M4_TG', 'OIH_M5_TG',
       'CON_V1_TG', 'CON_V2_TG', 'CON_V3_TG', 'CON_V4_TG', 'CON_V5_TG',
       'N1_TG', 'N2_TG', 'N3_TG','N4_TG', 'N5_TG', 
       'V6_TG','V7_TG', 'V8_TG', 'V9_TG', 'V10_TG']]

NA_geneTPM = geneTPM[['OIH_M1_NA', 'OIH_M2_NA', 'OIH_M3_NA', 'OIH_M4_NA', 'OIH_M5_NA',
       'CON_V1_NA', 'CON_V2_NA', 'CON_V3_NA', 'CON_V4_NA', 'CON_V5_NA',
       'N1_NA', 'N2_NA', 'N3_NA', 'N4_NA', 'N5_NA', 
       'V6_NA', 'V7_NA', 'V8_NA', 'V9_NA', 'V10_NA']]
TG_geneTPM = geneTPM[['OIH_M1_TG', 'OIH_M2_TG', 'OIH_M3_TG', 'OIH_M4_TG', 'OIH_M5_TG',
       'CON_V1_TG', 'CON_V2_TG', 'CON_V3_TG', 'CON_V4_TG', 'CON_V5_TG',
       'N1_TG', 'N2_TG', 'N3_TG','N4_TG', 'N5_TG', 
       'V6_TG','V7_TG', 'V8_TG', 'V9_TG', 'V10_TG']]

NA_gene = filter_gene(NA_geneCount,NA_geneTPM )
TG_gene = filter_gene(TG_geneCount,TG_geneTPM )

def changeIndex_TE(filtered_TE):
    filtered_TE['id']=filtered_TE.index
    filtered_TE['new_id']=filtered_TE['id'].apply(lambda x: x.split('_')[1])
    filtered_TE.index = filtered_TE['new_id']
    filtered_TE= filtered_TE.drop(columns = ["id", "new_id"])
    return(filtered_TE)
    
NA_geneCount = changeIndex_TE(NA_gene[0])
NA_geneTPM = changeIndex_TE(NA_gene[1])

TG_geneCount = changeIndex_TE(TG_gene[0])
TG_geneTPM = changeIndex_TE(TG_gene[1])

NA_geneCount.to_csv("filtered_NA_geneCount.csv") 
NA_geneTPM.to_csv("filtered_NA_geneTPM.csv") 

TG_geneCount.to_csv("filtered_TG_geneCount.csv") 
TG_geneTPM.to_csv("filtered_TG_geneTPM.csv") 

NA_isoCount = isoCount[['OIH_M1_NA', 'OIH_M2_NA', 'OIH_M3_NA', 'OIH_M4_NA', 'OIH_M5_NA',
       'CON_V1_NA', 'CON_V2_NA', 'CON_V3_NA', 'CON_V4_NA', 'CON_V5_NA',
       'N1_NA', 'N2_NA', 'N3_NA', 'N4_NA', 'N5_NA', 
       'V6_NA', 'V7_NA', 'V8_NA', 'V9_NA', 'V10_NA']]
TG_isoCount = isoCount[['OIH_M1_TG', 'OIH_M2_TG', 'OIH_M3_TG', 'OIH_M4_TG', 'OIH_M5_TG',
       'CON_V1_TG', 'CON_V2_TG', 'CON_V3_TG', 'CON_V4_TG', 'CON_V5_TG',
       'N1_TG', 'N2_TG', 'N3_TG','N4_TG', 'N5_TG', 
       'V6_TG','V7_TG', 'V8_TG', 'V9_TG', 'V10_TG']]
NA_isoTPM = isoCount[['OIH_M1_NA', 'OIH_M2_NA', 'OIH_M3_NA', 'OIH_M4_NA', 'OIH_M5_NA',
       'CON_V1_NA', 'CON_V2_NA', 'CON_V3_NA', 'CON_V4_NA', 'CON_V5_NA',
       'N1_NA', 'N2_NA', 'N3_NA', 'N4_NA', 'N5_NA', 
       'V6_NA', 'V7_NA', 'V8_NA', 'V9_NA', 'V10_NA']]
TG_isoTPM = isoCount[['OIH_M1_TG', 'OIH_M2_TG', 'OIH_M3_TG', 'OIH_M4_TG', 'OIH_M5_TG',
       'CON_V1_TG', 'CON_V2_TG', 'CON_V3_TG', 'CON_V4_TG', 'CON_V5_TG',
       'N1_TG', 'N2_TG', 'N3_TG','N4_TG', 'N5_TG', 
       'V6_TG','V7_TG', 'V8_TG', 'V9_TG', 'V10_TG']]

NA_iso = filter_gene(NA_isoCount,NA_isoTPM, "iso" )
TG_iso = filter_gene(TG_isoCount,TG_isoTPM, "iso" )


isoRatio = pd.read_table("OIH_NTG_isoforms_IR.txt", sep=" ",index_col=1)
isoRatio = filter_name(isoRatio, ENS_id)
isoRatio["gene_id"] = isoRatio.index 
isoRatio.index = isoRatio["transcript_id"]
isoRatio = isoRatio.drop(columns="transcript_id")
NA_isoRatio = isoRatio[['OIH_M1_NA', 'OIH_M2_NA', 'OIH_M3_NA', 'OIH_M4_NA', 'OIH_M5_NA',
       'CON_V1_NA', 'CON_V2_NA', 'CON_V3_NA', 'CON_V4_NA', 'CON_V5_NA',
       'IsoPct', 'IsoPct.1', 'IsoPct.2', 'IsoPct.3', 'IsoPct.4',
       'IsoPct.10','IsoPct.11', 'IsoPct.12', 'IsoPct.13', 'IsoPct.14']]
NA_isoRatio.columns = NA_isoTPM.columns
NA_isoRatio["gene_id"]= isoRatio.gene_id.tolist()
NA_isoRatio = NA_isoRatio.loc[NA_iso ]

TG_isoRatio =isoRatio[['OIH_M1_TG', 'OIH_M2_TG', 'OIH_M3_TG', 'OIH_M4_TG', 'OIH_M5_TG',
       'CON_V1_TG', 'CON_V2_TG', 'CON_V3_TG', 'CON_V4_TG', 'CON_V5_TG',
       'IsoPct.5','IsoPct.6', 'IsoPct.7', 'IsoPct.8', 'IsoPct.9',
       'IsoPct.15','IsoPct.16', 'IsoPct.17', 'IsoPct.18', 'IsoPct.19']] 
TG_isoRatio.columns = TG_isoTPM.columns
TG_isoRatio["gene_id"]= isoRatio.gene_id
TG_isoRatio = TG_isoRatio.loc[TG_iso]

#n_sample =20 
#Filter all the isoform of a gene with main isoform ratio > 95%
#The isoform with the lowest ratio is excluded to resolve the dependency among isoform ratio patterns within a gene
def IR_filter(IR, n_sample):
    df = pd.DataFrame(columns=IR.columns.tolist())
    IR["Mean"] = IR.iloc[:, 0:n_sample].mean(axis=1)
    IR = IR[(IR["Mean"] >0) & (IR["Mean"] <100) ]
    genes = IR.gene_id.unique().tolist()
    IR_grouped = IR.groupby("gene_id")
    for i in genes:
        gene_group = IR_grouped.get_group(i)
        if len(gene_group.index) >1:
            if (gene_group["Mean"].max() <95):
                minInd = gene_group.loc[gene_group["Mean"] == gene_group.loc[:,"Mean"].min()].index.tolist()
                gene_group = gene_group.drop(index = minInd)
                gene_group = gene_group.drop(["Mean"], axis=1)    
                df = df.append(gene_group)
    return(df)
    
filtered_NA_isoRatio = IR_filter(NA_isoRatio, 20)
filtered_TG_isoRatio = IR_filter(TG_isoRatio, 20)


def changeIndex_IR(filtered_IR):
    filtered_IR['id']=filtered_IR.index
    filtered_IR['new_id']=filtered_IR['id'].apply(lambda x: x.split('_')[0])
    filtered_IR.index = filtered_IR['new_id']
    filtered_IR= filtered_IR.drop(columns = ["gene_id","id", "new_id"])
    return(filtered_IR)
    
filtered_NA_isoRatio = changeIndex_IR(filtered_NA_isoRatio ) 
filtered_TG_isoRatio = changeIndex_IR(filtered_TG_isoRatio )
filtered_NA_isoRatio.to_csv("filtered_NA_isoRatio.csv") 
filtered_TG_isoRatio.to_csv("filtered_TG_isoRatio.csv")    

###################TSN#########################
def filter_gene(geneCount, geneTPM, t= "gene"):
    gcSum=(geneCount > 6).sum(axis=1).to_frame()
    gtSum=(geneTPM > 1).sum(axis=1).to_frame()
    gtSum.columns = ["gtSum"]
    gcSum.columns = ["gcSum"]
    geneCount = pd.concat([geneCount, pd.concat([gcSum,gtSum],axis=1)],axis=1)
    geneTPM = pd.concat([geneTPM, pd.concat([gcSum,gtSum],axis=1)],axis=1)   
    geneCount = filter_df(geneCount, 35)
    geneTPM = filter_df(geneTPM, 35)
    if t == "gene": 
        return (geneTPM)
genes_TSN = filter_gene(geneCount,geneTPM )
genes_TSN = changeIndex_TE(genes_TSN)
genes_TSN.to_csv("genes_TSN_35.csv") 

