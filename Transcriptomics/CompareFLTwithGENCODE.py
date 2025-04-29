#!/usr/bin/env python3

"""
Created on Wed Jun 29 11:14:43 2022
@author: Mengran Wang
In this code, we compare the full-length transcripts(FLT) from PacBio with GENCODE annotations and generate a comprehensive transcript table with ORF info. 
For FLT, before comparing with GENCODE, we also adjust the FLT by merging and get the longest transcript for transcripts with the same junction structure.
"""

import os
import sys
import re
import pandas as pd
import numpy as np
import pip
from itertools import chain

GENCODE_transcript_fa = pd.read_csv('/data/bio-wangmr/Genome/human/gencode.v37.pc_transcripts.txt',sep='\t',header=None)
GENCODE_aa_fa = pd.read_csv('/data/bio-wangmr/Genome/human/gencode.v37.pc_translations.txt',sep='\t',header=None)
ISO_fa = pd.read_csv('/data/bio-wangmr/ISO-seq/10_GENCODE_ORF/ISO_ovlp_GENCODE.bed',sep='\t',header=None)

############################################################## Functions start ##############################################################
# check and change the work directory
def set_workDir(workdir):
    os.getcwd() # check work directory
    os.chdir(workdir) # Change the current working directory
    # Print the current working directory
    print("Current working directory: {0}".format(os.getcwd()))

# function to return key for any value
def get_key(val, my_dict):
    for key, value in my_dict.items():
         if val == value:
             return key
    return "key doesn't exist"

# input the row(record) from genepred dataframe,
# reuturn list of junction sturctures (represented by location of junction sites)
def get_junc_structure(t_s, t_e):
    # t = df_PB.iloc[0]
    start = t_s.split(',')[ : -1]
    end = t_e.split(',')[ : -1]
    # link start and end to form [start1, end1, start2, end2, ...]
    junc_struc = [None]*(len(start)+len(end))
    junc_struc[::2] = start
    junc_struc[1::2] = end
    # remove the first start and the last end of the junc_struc (not care about non-slice-junction boundary)
    return junc_struc[1: -1]

# compare the junction structures of two transcripts
def is_share_structure(t1, t2): # t1 shorter than t2
    flag = False
    # judge which transctipt is longer
    # if len(t1) <= len(t2):
    #     list1 = t1
    #     list2 = t2
    # else:
    #     list1 = t2
    #     list2 = t1
    if len(t1) < len(t2):
        for i in range(len(t2) - len(t1) + 1):
            if t2[i: i+len(t1)] == t1:
                flag = True
                break
    return flag

# group and cluster the partially overlap strings
def group_partials(strings):
# check for subsets
    for i in range(len(strings)):
       for j in range(len(strings)):
          if i==j: continue # same index
          if (set(strings[i].split()) & set(strings[j].split())) == set(strings[i].split()): # if subset
             strings[i]="" # clear string
    b = []
    for x in strings:  # each string in a
       if len(x) > 0: # if not empty
          b.append(x)  # add to final list
    strings = b
    return strings

# For two transcripts with the same junction structure, merge and get the longest transcript
# compare the PB transcripts for the same gene, get longest unique transcripts
def get_longest_unique_transcripts(gene_id, gpd_PB):
    gpd_PB_gene = gpd_PB[gpd_PB['gene_ID'] == gene_id]
    # gpd_PB_gene = df_PB[df_PB['gene_ID'] == 'ENSG00000173482.17']
    longest_transID = []
    PB_strc_all = dict()
    # loop for PacBio transcripts
    for t in gpd_PB_gene.itertuples():
        t_start = t.exon_start; t_end = t.exon_end
        # get PacBio structure
        PB_strc = get_junc_structure(t_start, t_end)
        PB_strc_all[t.trans_ID] = ' '.join(PB_strc) # {trans_ID: structure}
        # PB_strc_all.append(' '.join(PB_strc))
        # compare and cluster structures
        longest = group_partials(list(PB_strc_all.values()))
    # get the trans_ID for longest and unqiue transcripts
    for l in longest:
        ID = get_key(l, PB_strc_all)
        longest_transID.append(ID)
    # return dataframe rows which contain longest transcripts
    return gpd_PB_gene[gpd_PB_gene['trans_ID'].isin(longest_transID)]

# first element of each sublist in a list of lists
def Extract(lst):
    return [item[0] for item in lst]

# compare the transcripts for the same gene
def compare_PBwithGENCODE(gene_id, gpd_PB, gpd_GENCODE):
    gpd_PB_gene = gpd_PB[gpd_PB['gene_ID'] == gene_id]
    gpd_GENCODE_gene = gpd_GENCODE[gpd_GENCODE['gene_ID'] == gene_id]
    #gpd_PB_gene = df_PB[df_PB['gene_ID'] == 'ENSG00000166813.17']
    #gpd_GENCODE_gene = df_G[df_G['gene_ID'] == 'ENSG00000166813.17']
    subs_info = []
    # loop for PacBio transcripts
    for t in gpd_PB_gene.itertuples():
        # subs_info_2 = []
        t_start = t.exon_start; t_end = t.exon_end
        # get PacBio structure
        PB_strc = get_junc_structure(t_start, t_end)

        found = False
        for T in gpd_GENCODE_gene.itertuples():
            T_start = T.exon_start; T_end = T.exon_end
            # get GENCODE structure
            G_strc = get_junc_structure(T_start, T_end)
            # compare PB with GENCODE
            #Flag = compare_transcript(PB_strc, G_strc)
            if is_share_structure(PB_strc, G_strc):
                found = True
                subs_info.append([t.trans_ID, T.trans_ID]) # [original, updated]
            else:
                #subs_info.append([t.trans_ID, t.trans_ID]) # keep the original
                subs_info.append(['Absent_in_PB', T.trans_ID]) # add the different GENCODE transcript
        # judge if all GENCODE transcripts have different stucture to PB transcript
        if found == False:
            subs_info.append([t.trans_ID, t.trans_ID])
        # subs_info = pd.concat([subs_info,pd.DataFrame(subs_info_2)])
    return pd.DataFrame(subs_info)

def get_fullSpliceMatch(sqanti_classification):
    # get recoreds which categorized as full-splice-match.
    df_FSM = sqanti_classification[(sqanti_classification['structural_category'] == 'full-splice_match')]
    # select FSM isoforms with TPM > 1
    #df_FSM_1 = df_FSM[df_FSM['iso_exp'] > 1]
    # output data frames
    info_FSM = pd.DataFrame()
    info_FSM = pd.concat([df_FSM['isoform'], df_FSM['associated_gene'],df_FSM['associated_transcript'],df_FSM['ORF_seq']], axis=1)
    return info_FSM

# get the ORF sequences from GENCODE_only transcripts
def get_ORF_for_GENCODE(trans_ID, df_G_ORF):
    df_out = df_G_ORF[df_G_ORF['trans_ID'].isin(trans_ID) == True]
    return df_out
# get the gene ID using transcript ID
def get_geneID_from_transID(trans_ID, df_G_expressed):
    df_out = df_G_expressed[df_G_expressed['trans_ID'].isin(trans_ID) == True]
    df_out = df_out[['trans_ID','gene_ID']]
    return df_out
# get the gene ID and ORF seq from GENCODE transcript ID.
def get_info_from_GENCODE_transID(trans_ID, info_sub1):
    trans_ID = list(info_sub1['ID_corrected'])
    # get the ORF info for trans_ID
    df1 = get_ORF_for_GENCODE(list(trans_ID), df_G_ORF)
    # get GENCODE gene ID for trans_ID
    df2 = get_geneID_from_transID(list(trans_ID), df_G_expressed)
    # merge df1 and df2
    df_merge = df2.merge(df1, on = 'trans_ID')
    return df_merge
# get the gene ID, trans ID and ORF seq, for genes with at least one PB transcript.
def get_info_geneWithPBs(info_all):
    # for transcripts that absent in PB, get the info from GENCODE
    info_sub1 = info_all[info_all['ID_before'] == 'Absent_in_PB'] # 26688 rows
    df_GENCODE_only_p2 = get_info_from_GENCODE_transID(list(info_sub1['ID_corrected']), info_sub1)
    df_GENCODE_only_p2['PacBio_ID'] = list(['Absent_in_PB'] * len(df_GENCODE_only_p2.index))
    df_GENCODE_only_p2['category'] = list(['GENCODE_only'] * len(df_GENCODE_only_p2.index))
    df_GENCODE_only_p2 = df_GENCODE_only_p2[['PacBio_ID','trans_ID','gene_ID','category','ORF_seq']]

    # for transcripts that exited in PB, and partially overlapped with GENCODE, get the info from GENCOE
    info_sub2 = info_all[(info_all['ID_before'] != 'Absent_in_PB') & (info_all['ID_before'] != info_all['ID_corrected'])] # 2461 rows
    df_GENCODE_corrected = get_info_from_GENCODE_transID(list(info_sub2['ID_corrected']), info_sub2)
    df_GENCODE_corrected = df_GENCODE_corrected.merge(info_sub2, left_on='trans_ID', right_on='ID_corrected')
    df_GENCODE_corrected['category'] = list(['GENCODE_from_PBcorrected'] * len(df_GENCODE_corrected.index))
    df_GENCODE_corrected = df_GENCODE_corrected[['ID_before','trans_ID','gene_ID','category','ORF_seq']]
    df_GENCODE_corrected.set_axis(['PacBio_ID','trans_ID','gene_ID','category','ORF_seq'], axis=1, inplace=True)

    # for transcripts that only in PB, get the info from PB dataset
    info_sub3 = info_all[info_all['ID_before'] == info_all['ID_corrected']] # 11552 rows
    df_PB_only = info_sub3.merge(df_PB, left_on = 'ID_before', right_on = 'trans_ID')
    df_PB_category_sub = df_PB_category[['isoform','ORF_seq']] # to get the ORF info
    df_PB_only = df_PB_only.merge(df_PB_category_sub, left_on = 'ID_before', right_on = 'isoform')
    df_PB_only['trans_ID'] = list(['Absent_in_GENCODE'] * len(df_PB_only.index))
    df_PB_only['category'] = list(['PB_only'] * len(df_PB_only.index))
    df_PB_only = df_PB_only[['ID_before', 'trans_ID', 'gene_ID','category','ORF_seq']]
    df_PB_only.set_axis(['PacBio_ID','trans_ID','gene_ID','category','ORF_seq'], axis=1, inplace=True)

    # merge all three parts
    df_geneWithPBs = pd.concat([df_GENCODE_only_p2, df_GENCODE_corrected, df_PB_only])
    return df_geneWithPBs

set_workDir('/Users/amdreamer/Documents/博后期间课题/ISO-seq/V2')

# df_PB = pd.read_csv('pbmm2Aln.collapsed.sqanti3_classification.filtered_lite.genePred2',sep="\t", header=None)
# df_G = pd.read_csv('gencode.v40.comprehensive.annotation.genePred2',sep="\t", header=None)
df_PB = pd.read_csv('pbmm2Aln.collapsed.sqanti3_classification.filtered_lite.genePred2',sep="\t", header=None)
df_G_expressed = pd.read_csv('gencode.v40.comprehensive.annotation.TPM_morethan1.genePred',sep="\t", header=None)
df_PB_category = pd.read_csv('pbmm2Aln.collapsed.sqanti3_classification.filtered_lite_classification.txt', sep="\t")
df_G_ORF = pd.read_csv('gencode.v40.pc_translations_tidy.txt', sep = "\t", header = None)

df_PB.set_axis(['trans_ID','chr','strand','TSS_start','TSS_end','CDS_start','CDS_end','exon_num','exon_start','exon_end','gene_ID'], axis=1, inplace=True)
df_G_expressed.set_axis(['trans_ID','chr','strand','TSS_start','TSS_end','CDS_start','CDS_end','exon_num','exon_start','exon_end','gene_ID','gene_name'], axis=1, inplace=True)
df_G_ORF.set_axis(['trans_ID','ORF_seq'], axis=1, inplace=True)

# len(df_PB) # 35960  len(df_G) # 246624

# first we need to remove the version ".xx" from gene ID and transcript ID
df_PB['gene_ID_2'] = df_PB['gene_ID'].str.split('.').str[0]; del df_PB['gene_ID']; df_PB['gene_ID'] = df_PB['gene_ID_2']; del df_PB['gene_ID_2']
df_PB_category['associated_gene_2'] = df_PB_category['associated_gene'].str.split('.').str[0]; del df_PB_category['associated_gene']; df_PB_category['associated_gene'] = df_PB_category['associated_gene_2']
df_PB_category['associated_transcript_2'] = df_PB_category['associated_transcript'].str.split('.').str[0]; del df_PB_category['associated_transcript']; df_PB_category['associated_transcript'] = df_PB_category['associated_transcript_2']

#df_PB2 = df_PB.iloc[0:1000]
#df_G2 = df_G.iloc[0:1000]

# get all PB isoforms which is full splice match to GENCODE, also are protein-coding isoforms.
df_PB_category = df_PB_category[df_PB_category['coding'] == 'coding'] # 33813 rows
df_PB_FSM = get_fullSpliceMatch(df_PB_category)
df_PB_FSM.set_axis(['PacBio_ID','gene_ID','trans_ID','ORF_seq'], axis=1, inplace=True)
# #  we need to remove the version ".xx" from gene ID and transcript ID
#df_PB_FSM['gene_ID_2'] = df_PB_FSM['gene_ID'].str.split('.').str[0]; del df_PB_FSM['gene_ID']; df_PB_FSM['gene_ID'] = df_PB_FSM['gene_ID_2']; del df_PB_FSM['gene_ID_2']
#df_PB_FSM['trans_ID_2'] = df_PB_FSM['trans_ID'].str.split('.').str[0]; del df_PB_FSM['trans_ID']; df_PB_FSM['trans_ID'] = df_PB_FSM['trans_ID_2']; del df_PB_FSM['trans_ID_2']


# get the structure(genepred) of all PB isoforms which is not full splice match to GENCODE
df_PB_notFSM = df_PB[df_PB['trans_ID'].isin(list(df_PB_FSM['PacBio_ID']))==False]

# to remove the retundancy, before compare the structure with PacBio and GENCODE,
# we need to remove the gencode transcript IDs that had the full-slice-match PacBio transcripts.
df_G_expressed_NoFSM = df_G_expressed[df_G_expressed['trans_ID'].isin(list(df_PB_FSM['trans_ID']))==False]

# loop for all FSM isoforms, keep all of them, write to the final result
#for gene_id in set(df_PB_FSM['gene_ID']):
    # get the longest and unique PB transcripts
#    PB_longest_FSM = get_longest_unique_transcripts(gene_id, df_PB_FSM_gdp)
    # filter the short transcripts and only keep the longest one

info_all = pd.DataFrame()
# loop for all genes in PacBio dataframes
for gene_id in set(df_PB_notFSM['gene_ID']):
    # get the longest and unique PB transcripts
    PB_longest = get_longest_unique_transcripts(gene_id, df_PB_notFSM)
    # compare longest and unique PB trans with GENCODE
    info_gene_id = compare_PBwithGENCODE(gene_id, PB_longest, df_G_expressed_NoFSM)
    # append the info for one gene to info_all
    info_all = pd.concat([info_all, info_gene_id])

# duplicated rows may generated, remove them
info_all.drop_duplicates(inplace=True)
# colnames the info_all
info_all.set_axis(['ID_before','ID_corrected'], axis=1, inplace=True)

# for genes that with at least one PB transcripts(means these genes are considered in PB correction process), add their info, labeled by conditions.
df_correction = get_info_geneWithPBs(info_all)

# for genes that are not considered in PB correction process,
# Part 1: for PB transcripts that are FSM with GENCODE, add their info, labeled as 'PB_&_GENCODE'
df_PB_FSM['category'] = list(['PB_&_GENCODE'] * len(df_PB_FSM.index))
df_PB_FSM = df_PB_FSM[['PacBio_ID','trans_ID','gene_ID','category','ORF_seq']]

# Part 2: for genes that without any of PB transcripts, add their info directly, labeled as "GENCODE_only"
## genes that occurs in info all, means these genes are considered in PB correction process.
gene_contains = df_G_expressed_NoFSM[df_G_expressed_NoFSM['trans_ID'].isin(list(info_all['ID_corrected']))==True]['gene_ID'].to_list()
## all genes - genes occurs = genes absent
gene_absent = set(df_G_expressed_NoFSM['gene_ID'].to_list()).difference(set(gene_contains))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:14:43 2022

@author: amdreamer
"""
import os
import sys
import pandas as pd
import numpy as np
from itertools import chain

# check and change the work directory
def set_workDir(workdir):
    os.getcwd() # check work directory
    os.chdir(workdir) # Change the current working directory
    # Print the current working directory
    print("Current working directory: {0}".format(os.getcwd()))

# function to return key for any value
def get_key(val, my_dict):
    for key, value in my_dict.items():
         if val == value:
             return key
    return "key doesn't exist"

# input the row(record) from genepred dataframe,
# reuturn list of junction sturctures (represented by location of junction sites)
def get_junc_structure(t_s, t_e):
    # t = df_PB.iloc[0]
    start = t_s.split(',')[ : -1]
    end = t_e.split(',')[ : -1]
    # link start and end to form [start1, end1, start2, end2, ...]
    junc_struc = [None]*(len(start)+len(end))
    junc_struc[::2] = start
    junc_struc[1::2] = end
    # remove the first start and the last end of the junc_struc (not care about non-slice-junction boundary)
    return junc_struc[1: -1]

# compare the junction structures of two transcripts
def is_share_structure(t1, t2): # t1 shorter than t2
    flag = False
    # judge which transctipt is longer
    # if len(t1) <= len(t2):
    #     list1 = t1
    #     list2 = t2
    # else:
    #     list1 = t2
    #     list2 = t1
    if len(t1) < len(t2):
        for i in range(len(t2) - len(t1) + 1):
            if t2[i: i+len(t1)] == t1:
                flag = True
                break
    return flag

# group and cluster the partially overlap strings
def group_partials(strings):
# check for subsets
    for i in range(len(strings)):
       for j in range(len(strings)):
          if i==j: continue # same index
          if (set(strings[i].split()) & set(strings[j].split())) == set(strings[i].split()): # if subset
             strings[i]="" # clear string
    b = []
    for x in strings:  # each string in a
       if len(x) > 0: # if not empty
          b.append(x)  # add to final list
    strings = b
    return strings

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 11:14:43 2022

@author: amdreamer
"""
import os
import sys
import pandas as pd
import numpy as np
from itertools import chain

# check and change the work directory
def set_workDir(workdir):
    os.getcwd() # check work directory
    os.chdir(workdir) # Change the current working directory
    # Print the current working directory
    print("Current working directory: {0}".format(os.getcwd()))

# function to return key for any value
def get_key(val, my_dict):
    for key, value in my_dict.items():
         if val == value:
             return key
    return "key doesn't exist"

# input the row(record) from genepred dataframe,
# reuturn list of junction sturctures (represented by location of junction sites)
def get_junc_structure(t_s, t_e):
    # t = df_PB.iloc[0]
    start = t_s.split(',')[ : -1]
    end = t_e.split(',')[ : -1]
    # link start and end to form [start1, end1, start2, end2, ...]
    junc_struc = [None]*(len(start)+len(end))
    junc_struc[::2] = start
    junc_struc[1::2] = end
    # remove the first start and the last end of the junc_struc (not care about non-slice-junction boundary)
    return junc_struc[1: -1]

# compare the junction structures of two transcripts
def is_share_structure(t1, t2): # t1 shorter than t2
    flag = False
    # judge which transctipt is longer
    # if len(t1) <= len(t2):
    #     list1 = t1
    #     list2 = t2
    # else:
    #     list1 = t2
    #     list2 = t1
    if len(t1) < len(t2):
        for i in range(len(t2) - len(t1) + 1):
            if t2[i: i+len(t1)] == t1:
                flag = True
                break
    return flag

# group and cluster the partially overlap strings
def group_partials(strings):
# check for subsets
    for i in range(len(strings)):
       for j in range(len(strings)):
          if i==j: continue # same index
          if (set(strings[i].split()) & set(strings[j].split())) == set(strings[i].split()): # if subset
             strings[i]="" # clear string
    b = []
    for x in strings:  # each string in a
       if len(x) > 0: # if not empty
          b.append(x)  # add to final list
    strings = b
    return strings

# For two transcripts with the same junction structure, merge and get the longest transcript
# compare the PB transcripts for the same gene, get longest unique transcripts
def get_longest_unique_transcripts(gene_id, gpd_PB):
    gpd_PB_gene = gpd_PB[gpd_PB['gene_ID'] == gene_id]
    # gpd_PB_gene = df_PB[df_PB['gene_ID'] == 'ENSG00000173482.17']
    longest_transID = []
    PB_strc_all = dict()
    # loop for PacBio transcripts
    for t in gpd_PB_gene.itertuples():
        t_start = t.exon_start; t_end = t.exon_end
        # get PacBio structure
        PB_strc = get_junc_structure(t_start, t_end)
        PB_strc_all[t.trans_ID] = ' '.join(PB_strc) # {trans_ID: structure}
        # PB_strc_all.append(' '.join(PB_strc))
        # compare and cluster structures
        longest = group_partials(list(PB_strc_all.values()))
    # get the trans_ID for longest and unqiue transcripts
    for l in longest:
        ID = get_key(l, PB_strc_all)
        longest_transID.append(ID)
    # return dataframe rows which contain longest transcripts
    return gpd_PB_gene[gpd_PB_gene['trans_ID'].isin(longest_transID)]

# first element of each sublist in a list of lists
def Extract(lst):
    return [item[0] for item in lst]

# compare the transcripts for the same gene
def compare_PBwithGENCODE(gene_id, gpd_PB, gpd_GENCODE):
    gpd_PB_gene = gpd_PB[gpd_PB['gene_ID'] == gene_id]
    gpd_GENCODE_gene = gpd_GENCODE[gpd_GENCODE['gene_ID'] == gene_id]
    #gpd_PB_gene = df_PB[df_PB['gene_ID'] == 'ENSG00000166813.17']
    #gpd_GENCODE_gene = df_G[df_G['gene_ID'] == 'ENSG00000166813.17']
    subs_info = []
    # loop for PacBio transcripts
    for t in gpd_PB_gene.itertuples():
        # subs_info_2 = []
        t_start = t.exon_start; t_end = t.exon_end
        # get PacBio structure
        PB_strc = get_junc_structure(t_start, t_end)

        found = False
        for T in gpd_GENCODE_gene.itertuples():
            T_start = T.exon_start; T_end = T.exon_end
            # get GENCODE structure
            G_strc = get_junc_structure(T_start, T_end)
            # compare PB with GENCODE
            #Flag = compare_transcript(PB_strc, G_strc)
            if is_share_structure(PB_strc, G_strc):
                found = True
                subs_info.append([t.trans_ID, T.trans_ID]) # [original, updated]
            else:
                #subs_info.append([t.trans_ID, t.trans_ID]) # keep the original
                subs_info.append(['Absent_in_PB', T.trans_ID]) # add the different GENCODE transcript
        # judge if all GENCODE transcripts have different stucture to PB transcript
        if found == False:
            subs_info.append([t.trans_ID, t.trans_ID])
        # subs_info = pd.concat([subs_info,pd.DataFrame(subs_info_2)])
    return pd.DataFrame(subs_info)

def get_fullSpliceMatch(sqanti_classification):
    # get recoreds which categorized as full-splice-match.
    df_FSM = sqanti_classification[(sqanti_classification['structural_category'] == 'full-splice_match')]
    # select FSM isoforms with TPM > 1
    #df_FSM_1 = df_FSM[df_FSM['iso_exp'] > 1]
    # output data frames
    info_FSM = pd.DataFrame()
    info_FSM = pd.concat([df_FSM['isoform'], df_FSM['associated_gene'],df_FSM['associated_transcript'],df_FSM['ORF_seq']], axis=1)
    return info_FSM

# get the ORF sequences from GENCODE_only transcripts
def get_ORF_for_GENCODE(trans_ID, df_G_ORF):
    df_out = df_G_ORF[df_G_ORF['trans_ID'].isin(trans_ID) == True]
    return df_out
# get the gene ID using transcript ID
def get_geneID_from_transID(trans_ID, df_G_expressed):
    df_out = df_G_expressed[df_G_expressed['trans_ID'].isin(trans_ID) == True]
    df_out = df_out[['trans_ID','gene_ID']]
    return df_out

# get the gene ID and ORF seq from GENCODE transcript ID.
def get_info_from_GENCODE_transID(trans_ID, info_sub1):
    trans_ID = list(info_sub1['ID_corrected'])
    # get the ORF info for trans_ID
    df1 = get_ORF_for_GENCODE(list(trans_ID), df_G_ORF)
    # get GENCODE gene ID for trans_ID
    df2 = get_geneID_from_transID(list(trans_ID), df_G_expressed)
    # merge df1 and df2
    df_merge = df2.merge(df1, on = 'trans_ID')
    return df_merge
# get the gene ID, trans ID and ORF seq, for genes with at least one PB transcript.
def get_info_geneWithPBs(info_all):
    # for transcripts that absent in PB, get the info from GENCODE
    info_sub1 = info_all[info_all['ID_before'] == 'Absent_in_PB'] # 26688 rows
    df_GENCODE_only_p2 = get_info_from_GENCODE_transID(list(info_sub1['ID_corrected']), info_sub1)
    df_GENCODE_only_p2['PacBio_ID'] = list(['Absent_in_PB'] * len(df_GENCODE_only_p2.index))
    df_GENCODE_only_p2['category'] = list(['GENCODE_only'] * len(df_GENCODE_only_p2.index))
    df_GENCODE_only_p2 = df_GENCODE_only_p2[['PacBio_ID','trans_ID','gene_ID','category','ORF_seq']]

    # for transcripts that exited in PB, and partially overlapped with GENCODE, get the info from GENCOE
    info_sub2 = info_all[(info_all['ID_before'] != 'Absent_in_PB') & (info_all['ID_before'] != info_all['ID_corrected'])] # 2461 rows
    df_GENCODE_corrected = get_info_from_GENCODE_transID(list(info_sub2['ID_corrected']), info_sub2)
    df_GENCODE_corrected = df_GENCODE_corrected.merge(info_sub2, left_on='trans_ID', right_on='ID_corrected')
    df_GENCODE_corrected['category'] = list(['GENCODE_from_PBcorrected'] * len(df_GENCODE_corrected.index))
    df_GENCODE_corrected = df_GENCODE_corrected[['ID_before','trans_ID','gene_ID','category','ORF_seq']]
    df_GENCODE_corrected.set_axis(['PacBio_ID','trans_ID','gene_ID','category','ORF_seq'], axis=1, inplace=True)

    # for transcripts that only in PB, get the info from PB dataset
    info_sub3 = info_all[info_all['ID_before'] == info_all['ID_corrected']] # 11552 rows
    df_PB_only = info_sub3.merge(df_PB, left_on = 'ID_before', right_on = 'trans_ID')
    df_PB_category_sub = df_PB_category[['isoform','ORF_seq']] # to get the ORF info
    df_PB_only = df_PB_only.merge(df_PB_category_sub, left_on = 'ID_before', right_on = 'isoform')
    df_PB_only['trans_ID'] = list(['Absent_in_GENCODE'] * len(df_PB_only.index))
    df_PB_only['category'] = list(['PB_only'] * len(df_PB_only.index))
    df_PB_only = df_PB_only[['ID_before', 'trans_ID', 'gene_ID','category','ORF_seq']]
    df_PB_only.set_axis(['PacBio_ID','trans_ID','gene_ID','category','ORF_seq'], axis=1, inplace=True)

    # merge all three parts
    df_geneWithPBs = pd.concat([df_GENCODE_only_p2, df_GENCODE_corrected, df_PB_only])
    return df_geneWithPBs
############################################################## Functions end ##############################################################

set_workDir('./ISO-seq/V2')

############################################################## Analysis start ##############################################################
df_PB = pd.read_csv('pbmm2Aln.collapsed.sqanti3_classification.filtered_lite.genePred2',sep="\t", header=None)
df_G_expressed = pd.read_csv('gencode.v40.comprehensive.annotation.TPM_morethan1.genePred',sep="\t", header=None)
df_PB_category = pd.read_csv('pbmm2Aln.collapsed.sqanti3_classification.filtered_lite_classification.txt', sep="\t")
df_G_ORF = pd.read_csv('gencode.v40.pc_translations_tidy.txt', sep = "\t", header = None)

df_PB.set_axis(['trans_ID','chr','strand','TSS_start','TSS_end','CDS_start','CDS_end','exon_num','exon_start','exon_end','gene_ID'], axis=1, inplace=True)
df_G_expressed.set_axis(['trans_ID','chr','strand','TSS_start','TSS_end','CDS_start','CDS_end','exon_num','exon_start','exon_end','gene_ID','gene_name'], axis=1, inplace=True)
df_G_ORF.set_axis(['trans_ID','ORF_seq'], axis=1, inplace=True)

# len(df_PB) # 35960  len(df_G) # 246624

# first we need to remove the version ".xx" from gene ID and transcript ID
df_PB['gene_ID_2'] = df_PB['gene_ID'].str.split('.').str[0]; del df_PB['gene_ID']; df_PB['gene_ID'] = df_PB['gene_ID_2']; del df_PB['gene_ID_2']
df_PB_category['associated_gene_2'] = df_PB_category['associated_gene'].str.split('.').str[0]; del df_PB_category['associated_gene']; df_PB_category['associated_gene'] = df_PB_category['associated_gene_2']
df_PB_category['associated_transcript_2'] = df_PB_category['associated_transcript'].str.split('.').str[0]; del df_PB_category['associated_transcript']; df_PB_category['associated_transcript'] = df_PB_category['associated_transcript_2']

#df_PB2 = df_PB.iloc[0:1000]
#df_G2 = df_G.iloc[0:1000]

# get all PB isoforms which is full splice match to GENCODE, also are protein-coding isoforms.
df_PB_category = df_PB_category[df_PB_category['coding'] == 'coding'] # 33813 rows
df_PB_FSM = get_fullSpliceMatch(df_PB_category)
df_PB_FSM.set_axis(['PacBio_ID','gene_ID','trans_ID','ORF_seq'], axis=1, inplace=True)
# #  we need to remove the version ".xx" from gene ID and transcript ID
#df_PB_FSM['gene_ID_2'] = df_PB_FSM['gene_ID'].str.split('.').str[0]; del df_PB_FSM['gene_ID']; df_PB_FSM['gene_ID'] = df_PB_FSM['gene_ID_2']; del df_PB_FSM['gene_ID_2']
#df_PB_FSM['trans_ID_2'] = df_PB_FSM['trans_ID'].str.split('.').str[0]; del df_PB_FSM['trans_ID']; df_PB_FSM['trans_ID'] = df_PB_FSM['trans_ID_2']; del df_PB_FSM['trans_ID_2']


# get the structure(genepred) of all PB isoforms which is not full splice match to GENCODE
df_PB_notFSM = df_PB[df_PB['trans_ID'].isin(list(df_PB_FSM['PacBio_ID']))==False]

# to remove the retundancy, before compare the structure with PacBio and GENCODE,
# we need to remove the gencode transcript IDs that had the full-slice-match PacBio transcripts.
df_G_expressed_NoFSM = df_G_expressed[df_G_expressed['trans_ID'].isin(list(df_PB_FSM['trans_ID']))==False]

# loop for all FSM isoforms, keep all of them, write to the final result
#for gene_id in set(df_PB_FSM['gene_ID']):
    # get the longest and unique PB transcripts
#    PB_longest_FSM = get_longest_unique_transcripts(gene_id, df_PB_FSM_gdp)
    # filter the short transcripts and only keep the longest one

info_all = pd.DataFrame()
# loop for all genes in PacBio dataframes
for gene_id in set(df_PB_notFSM['gene_ID']):
    # get the longest and unique PB transcripts
    PB_longest = get_longest_unique_transcripts(gene_id, df_PB_notFSM)
    # compare longest and unique PB trans with GENCODE
    info_gene_id = compare_PBwithGENCODE(gene_id, PB_longest, df_G_expressed_NoFSM)
    # append the info for one gene to info_all
    info_all = pd.concat([info_all, info_gene_id])

# duplicated rows may generated, remove them
info_all.drop_duplicates(inplace=True)
# colnames the info_all
info_all.set_axis(['ID_before','ID_corrected'], axis=1, inplace=True)

# for genes that with at least one PB transcripts(means these genes are considered in PB correction process), add their info, labeled by conditions.
df_correction = get_info_geneWithPBs(info_all)

# for genes that are not considered in PB correction process,
# Part 1: for PB transcripts that are FSM with GENCODE, add their info, labeled as 'PB_&_GENCODE'
df_PB_FSM['category'] = list(['PB_&_GENCODE'] * len(df_PB_FSM.index))
df_PB_FSM = df_PB_FSM[['PacBio_ID','trans_ID','gene_ID','category','ORF_seq']]

# Part 2: for genes that without any of PB transcripts, add their info directly, labeled as "GENCODE_only"
## genes that occurs in info all, means these genes are considered in PB correction process.
gene_contains = df_G_expressed_NoFSM[df_G_expressed_NoFSM['trans_ID'].isin(list(info_all['ID_corrected']))==True]['gene_ID'].to_list()
## all genes - genes occurs = genes absent
gene_absent = set(df_G_expressed_NoFSM['gene_ID'].to_list()).difference(set(gene_contains))
trans_absent = df_G_expressed_NoFSM[df_G_expressed_NoFSM['gene_ID'].isin(gene_absent) == True]['trans_ID'].to_list()
# merge the trans_ID, gene_ID, and ORF seq for genes that absent in
df1 = get_ORF_for_GENCODE(list(trans_absent), df_G_ORF)
df2 = get_geneID_from_transID(list(trans_absent), df_G_expressed)
df_GENCODE_only = df2.merge(df1, on = 'trans_ID')
df_GENCODE_only['PacBio_ID'] = list(['Absent_in_PB'] * len(df_GENCODE_only.index))
df_GENCODE_only['category'] = list(['GENCODE_only'] * len(df_GENCODE_only.index))
df_GENCODE_only = df_GENCODE_only[['PacBio_ID','trans_ID','gene_ID','category','ORF_seq']]

# in summary, we combinded all three parts, which include all expressed protein coding genes in our cell line. (non-coding and low expressed are filtered)
df_all = pd.concat([df_correction, df_PB_FSM, df_GENCODE_only]) # 49421 transcripts


# write to out put
pd.DataFrame(info_all).to_csv('PB_withGENCODEcorrected_2.txt', index=None, sep="\t")
df_all.to_csv('protein_coding_transcript_info.txt', index=None, sep="\t")

# get unqiue ORF sequences
ORF_list = set(df_all['ORF_seq'].to_list()) #  32612 ORFs
ORF_all = []; joinInfo_ORF_all = []; ORF_category = []
for ORF in ORF_list:
    # get the transcripts that match the ORF
    # df_all[df_all['ORF_seq'] == ORF]
    # set a blank list to store info for each ORF
    info_ORF = []
    # for each ORF, iter rows(transcript)
    for index, row in df_all[df_all['ORF_seq'] == ORF].iterrows():
        #print(row['PacBio_ID'],row['trans_ID'],row['gene_ID'],row['category'])
        # join info for each transcript
        joinInfo_transcript = ','.join([row['PacBio_ID'],row['trans_ID'],row['gene_ID'],row['category']])
        # append the info for each trans into info_ORF
        info_ORF.append(joinInfo_transcript)
    # join info for each ORF
    joinInfo_ORF = '|'.join(info_ORF)
    # judge the ORF category with priority: PB_&_GENCODE > GENCODE_from_PBcorrected > PB_only > GENCODE_only
    if joinInfo_ORF.find('PB_&_GENCODE') != -1:
        cate = 'PB_&_GENCODE'
    elif (joinInfo_ORF.find('PB_only') != -1) & (joinInfo_ORF.find('GENCODE_only') != -1): # means found PB and GENCODE for the same ORF
        cate = 'PB_&_GENCODE'
    elif joinInfo_ORF.find('GENCODE_from_PBcorrected') != -1:
        cate = 'GENCODE_from_PBcorrected'
    elif joinInfo_ORF.find('PB_only') != -1:
        cate = 'PB_only'
    else:
        cate = 'GENCODE_only'

    # store ORF and Info_ORF into two list
    ORF_all.append(ORF); joinInfo_ORF_all.append(joinInfo_ORF); ORF_category.append(cate)

# output ORF-centric results
df_all_ORF = pd.DataFrame({'transcript_info':joinInfo_ORF_all, 'ORF_category':ORF_category, 'ORF_seq':ORF_all})

df_all_ORF.to_csv('protein_coding_transcript_ORFinfo.txt', index=None, sep="\t")
