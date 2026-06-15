#!/usr/bin/env python3
"""
Created on Thu Aug 22 09:57:08 2024

@author: iaengels
"""
#Classicol_version_2_0_0.py
import argparse
import time
import os
import numpy as np
import pandas as pd
import csv
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio import Align
import warnings
warnings.filterwarnings("ignore")     
warnings.filterwarnings("ignore", category=DeprecationWarning)
import re
import plotly
import plotly.graph_objs as go
import plotly.figure_factory as ff
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
import maxquant
import random
from scipy.spatial.distance import braycurtis
from Bio.Align import substitution_matrices
import plotly.express as px
import pathlib
from datetime import date
import warnings
import typing
warnings.filterwarnings("ignore")     
warnings.filterwarnings("ignore", category=DeprecationWarning)
# from numba import jit, types
# from numba.typed import List
from plotly.subplots import make_subplots


def retreive_uncertain() -> dict[str, list[str]]:
    return {'B':['D','N'],'Z':['Q','E'],'X':['A','C','D','E','F','G','H','I','L','K','M','N','P','Q','R','S','T','V','W']}

def retreive_AA_codes()-> dict[str, float]:
    return {"A" : 71.03711,     "C" : 103.00919,    "D" : 115.02694,    "E" : 129.04259,    "F" : 147.06841,
            "G" : 57.02146,     "H" : 137.05891,    "L" : 113.08406,    "M" : 131.04049,    "K" : 128.09496,
            "N" : 114.04293,    "P" : 97.05276,     "Q" : 128.05858,    "S" : 87.03203,     "R" : 156.10111, 
            "T" : 101.04768,    "V" : 99.06841,     "W" : 186.07931,    "Y" : 163.06333,    "I" : 113.08406,
            "U" : 168.05,       "&" : 1000}

def crap_f(path: pathlib.Path, manual_fasta: pathlib.Path | None) -> dict[Seq, str]:
    print('making database')
    crap = {}
    files_own = []
    for fastafile in os.walk(path /'BoneDB'):
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta') or i.endswith('.fa'):
                files_own.append(path/'BoneDB'/i)
    if add_fasta != None:
        for fastafile in os.walk(add_fasta):
            for i in fastafile[-1]:
                if i.endswith('.txt') or i.endswith('.fasta') or i.endswith('.fa'):
                    files_own.append(add_fasta/i)
    done = ''
    for i in files_own:
        print(f"Processing '{i.name}'")
        for record in SeqIO.parse(i, "fasta"):
            if str(record.description) in done:
                continue
            if 'J' in record.seq or 'O' in record.seq:
                print('skip',record.description)
                continue
            add = record.seq
            while add in crap:
                add = Seq(str(add)+'&') #if seq the same as relative
            done += str(record.description)
            crap[add]=record.description
    return crap

def do_unimod(path,ptms):
    ptms = list(ptms)
    raw_unimod = pd.DataFrame()
    ptm_list = []
    ptms.append('Deamidated (NQ)')
    ambig_PTM = []
    for ptm_in in ptms:
        if len(ptm_in)>0:
            for ptm in ptm_in.split('; '):#multiple per peptide possible, do not miss any
                if len(ptm.split(' '))==3:#remove numbers
                    ptm = ' '.join(ptm.split(' ')[1:])
                if ')' in ptm:#some PTMs are not return on a specific residue, like Asp->Glu
                    temp=ptm.split('(')
                    if 'N-term' in temp[1].split(')')[0]:
                        ptm_list.append((temp[0].replace(' ',''),'N-term'))
                    if 'C-term' in temp[1].split(')')[0]:
                        ptm_list.append((temp[0].replace(' ',''),'C-term'))
                    for a in temp[1].split(')')[0]:
                        ptm_list.append((temp[0].replace(' ',''),a))
                else:
                    ambig_PTM.append(ptm.split(' ')[-1])#in case of numbers
    ptm_list = set(ptm_list)
    ambig_PTM = set(ambig_PTM)
    with open(path/'MISC'/'unimod.txt') as r:
        for line in r:
            line = line.replace('=', ',')
            line = line.strip()
            array = np.array(line.split(',')).reshape(1,-1)
            df_add = pd.DataFrame(array)
            raw_unimod = pd.concat([raw_unimod, df_add], ignore_index=True)
    AAs = []
    for aa,locaa in raw_unimod[[1,4]].values:
        
        if '[' in aa:
            i = aa.split(']')
        else:
            AAs.append(aa)
            continue
        if ''.join(c for c in i[1] if c.isdigit()==False) != i[1]:
            add = ''.join(c for c in i[1] if c.isdigit()==False)+'_!'#can change ! with locaa
            while add in AAs:
                add = add+'!'
            AAs.append(add)
        else:
            AAs.append(i[1])
    raw_unimod[1]=np.array(AAs)
    temp_unimod = pd.DataFrame()
    temp_unimod['PTM']=raw_unimod[1].values[2:]
    temp_unimod['delta_mass']=raw_unimod[2].values[2:]
    unimod = {}
    for i,u in temp_unimod[['PTM', 'delta_mass']].values:
        unimod[i]=u

    unimod['AA']='0'
    for i, u in AA_codes.items():
        unimod[i]=str(u)
    unimod_db = raw_unimod.iloc[2:]
    unimod_db = unimod_db.drop([0,3], axis=1)
    unimod_db.columns = ['PTM', 'mass', 'AA', 'type']
    unimod_db['mass']=np.array([float(m) for m in unimod_db['mass'].values])
    unimod_db['PTM']=np.array([num.split(']')[-1] for num in unimod_db['PTM'].values])
    unimod_db=unimod_db[unimod_db['PTM'] != '']
    
    unimod_db=unimod_db[unimod_db['type'] != 'Manual']
    drop = [False if ''.join(c for c in element if c.isdigit()==False) != element else True for element in unimod_db['PTM'].values]
    unimod_db['digit']=np.array(drop)
    unimod_db=unimod_db[unimod_db['digit'] == True]
    
    drop = []
    for p,aa in unimod_db[['PTM','AA']].values:
        added=False
        for ptm,aas in set(ptm_list): #We only check for ptms that mascot searched for
            if ptm==p and aas==aa:
                drop.append(True)
                added=True
                break
        if added == False:
            if p in ambig_PTM:
                drop.append(True)
                added=True
                break
        if added==False:
            drop.append(False)
    unimod_db['digit']=np.array(drop)
    unimod_db=unimod_db[unimod_db['digit'] == True]
    return unimod_db, unimod

def find_mass(peptide,AA_codes) -> float: #sums the masses of the amino acid sequence inputed 
    return float(sum(AA_codes[AA] for AA in peptide))    

def make_matrix(codes):
    doubles = [] #making all pairs of amino acids, here the ptms are also included
    for element1 in list(codes.keys()):
        for element2 in codes.keys():
            doubles.append(element1+'|'+element2)
    
    names = list(codes.keys())+doubles#Add together all amino acids, ptms, doubles and if wanted triples, although the triples take long to compute.
    reduced_matrix = []
    for element in names:
        m_el1=find_mass(element.split('|'),codes)
        for element2 in names:
            m_el2=find_mass(element2.split('|'),codes)
            if -0.01<=(m_el2-m_el1)<=0.01:
                if element != element2:
                    reduced_matrix.append((element2,element,m_el2-m_el1))        
    return reduced_matrix

def load_files_peaks(path,name_file,AA_codes):
    print('open file {}'.format(name_file))
    df = pd.read_csv(name_file,header=0)
    peps = []
    pvm = []
    pvm_loc = []
    pst = []
    protein = []
    df = df.fillna('')
    for i,p,scan,pro in df[['Peptide','AScore','Scan','Accession']].values:
        pst.append('Peaks_result_peptide_'+str(scan))
        protein.append(pro)
        i = ''.join([el for el in i if el.isalpha()==True])
        peps.append(i)
        if len(p)==0:
            pvm.append('')
            pvm_loc.append('0.'+'0'*len(i)+'.0')
            continue
        p = p.split(';')
        temp_ptm = {}
        aa_loc = {}
        pvm_pos = [0]*len(i)
        nterm = 0
        cterm = 0
        for ptm in p:
            ptm = ptm.split(':')
            
            if ptm[0]=='N-term':
                nterm=1
                continue
            if ptm[0]=='C-term':
                cterm=1
                continue
            locs = int([el for el in ptm[0] if el.isdigit()==True][0])
            ptm_1 = ptm[1]
            if '(' not in ptm_1:
                ptm_1 += ' ('+i[locs-1]+')'
            if ptm_1 in temp_ptm:
                temp_ptm[ptm_1]=temp_ptm[ptm_1]+1
                aa_loc[ptm_1].extend([ptm_1])
            else:
                temp_ptm[ptm_1] = 1
                aa_loc[ptm_1] = [ptm_1]
            pvm_pos[locs-1]=1
        pvm_loc.append(str(nterm)+'.'+''.join([str(el) for el in pvm_pos])+'.'+str(cterm))
        combo_ptm = []
        for k,v in temp_ptm.items():
            combo_ptm.append(str(v)+' '+k)
        pvm.append(';'.join(combo_ptm))
            
    df['pep_seq']=peps
    df['pep_var_mod']=pvm
    df['pep_var_mod_pos']=pvm_loc
    df['prot_desc']=protein
    df['pep_scan_title']=pst
    df = df.drop_duplicates()
    
    
    
    df_4_uni= df
    
    df2=df
    unimod_db, unimod = do_unimod(path,df_4_uni['pep_var_mod'].values)
    ids = {}
    for p, a,m in unimod_db[['PTM','AA','mass']].values:
        if a=='N-term':
            a='!'
        elif a=='C-term':
            a='*'
        add = '?'
        while add+a in AA_codes.keys():
            add+='?'
        if a=='!' or a=='*':
            AA_codes[add+a]=float(m)
        else:
            AA_codes[add+a]=AA_codes[a]+float(m)
        ids[add+a]=p
    
    df = df[['prot_desc','pep_seq','pep_var_mod','pep_var_mod_pos','pep_scan_title']]
    return df2, unimod_db, ids,AA_codes

def load_files_winnow(path,name_file,AA_codes,overall_ptm):
    print('open file {}'.format(name_file))
    df = pd.read_excel(name_file, header=0)
    df = df.fillna('')
    df['prot_desc']=['unknown']*len(df)
    name = list(df['experiment_name'].values)[0]
    df['pep_scan_title']=[name+'_'+str(scan) for scan in df['scan_number'].values]
    #turn modification in the style of mascot
    #'pep_var_mod_pos' 0.001100.0
    #'pep_var_mod' Oxidation (P); Deamidated(NQ)
    new_mod = []
    new_mod_pos = []
    for p in df['preds'].values:
        if '(' not in p:
            new_mod.append('')
            new_mod_pos.append('0.'+'0'*len(p)+'.0')
        else:
            pvm = []
            adj_loc = 1
            ptm_found = False
            new_ptm = ''
            for loc,i in enumerate(p):
                if '('==i:
                    adj_loc -= 1
                    pvm.append(loc+adj_loc)
                    ptm_found = True
                elif i==')':
                    adj_loc -= 1
                    ptm_found = False
                    if new_ptm not in overall_ptm:
                        overall_ptm[new_ptm]=input('{} has been found as PTM, please provide the Unimod annotation for this PTM to continue')
                    new_ptm = overall_ptm[new_ptm]
                    pvm.append(new_ptm)
                    new_ptm=''
                elif ptm_found==True:
                    new_ptm += i
            new_pvm = {}
            for i in range(0,len(pvm),2):
                new_pvm[int(pvm[i])]=pvm[i+1]
            temp_mod={}
            adj_mod_pos=['0']*len(p)
            n_term='0.'
            c_term='.0'
            for k,v in new_pvm.items():
                if k==0:
                    n_term='1.'
                    v=v+' (N-term)'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
                elif k==-1:
                    c_term='.1'
                    v=v+' (C-term)'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
                else:
                    adj_mod_pos[k-1]='1'
                    v=v+' ('+p[k-1]+')'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
            adj_mod = '; '.join([str(v)+' '+k if v>1 else k for k,v in temp_mod.items()])
            new_mod.append(adj_mod)
            new_mod_pos.append(n_term+''.join(adj_mod_pos)+c_term)
    df['pep_var_mod']=new_mod
    df['pep_var_mod_pos']=new_mod_pos
    
    df['pep_seq']=df['preds']
    df_4_uni= df
    
    df2=df
    unimod_db, unimod = do_unimod(path,df_4_uni['pep_var_mod'].values)
    ids = {}
    for p, a,m in unimod_db[['PTM','AA','mass']].values:
        if a=='N-term':
            a='!'
        elif a=='C-term':
            a='*'
        add = '?'
        while add+a in AA_codes.keys():
            add+='?'
        if a=='!' or a=='*':
            AA_codes[add+a]=float(m)
        else:
            AA_codes[add+a]=AA_codes[a]+float(m)
        ids[add+a]=p
    
    df = df[['prot_desc','pep_seq','pep_var_mod','pep_var_mod_pos','pep_scan_title']]
    return df2, unimod_db, ids,AA_codes, overall_ptm

def load_files_mascot(path, name_file,AA_codes):
    print('open file {}'.format(name_file))
    found = False
    header_row = 0
    while found == False:
        try:
            df = pd.read_csv(name_file, header=header_row)
            if 'pep_seq' in df.columns:
                found = True
            else:
                header_row += 1
        except:
            header_row += 1
            if header_row >1000:
                break
    df = df.fillna('')
    charges = df['pep_exp_z'].values
    df = df[['prot_desc','pep_seq','pep_var_mod','pep_var_mod_pos','pep_scan_title']]
    df['pep_scan_title']=[num.replace('~', '"') for num in df['pep_scan_title'].values]
    df_4_uni= df
    
    df2=df
    unimod_db, unimod = do_unimod(path,df_4_uni['pep_var_mod'].values)
    ids = {}
    for p, a,m in unimod_db[['PTM','AA','mass']].values:
        if a=='N-term':
            a='!'
        elif a=='C-term':
            a='*'
        add = '?'
        while add+a in AA_codes.keys():
            add+='?'
        if a=='!' or a=='*':
            AA_codes[add+a]=float(m)
        else:
            AA_codes[add+a]=AA_codes[a]+float(m)
        ids[add+a]=p
    return df2, unimod_db, ids,AA_codes

def maxquant_bulk(name_file):
    print('Checking for bulk maxquant')
    #name_file = name of the file benchmark
    MQ_file = maxquant.io.read_maxquant(name_file)
    bulk = set(list(MQ_file['Raw file'].values))
    return list(bulk)

def load_files_maxquant(path, name_file,experiment,AA_codes):
    print('open maxquant file for experiment {}'.format(experiment))
    #name_file = name of the file benchmark
    MQ_file = maxquant.io.read_maxquant(name_file)
    MQ_file = MQ_file[['Sequence','Modified sequence','Raw file','Modifications']][MQ_file['Raw file']==experiment]
    MQ_file.columns = ['pep_seq','pep_var_mod_pos','pep_scan_title','mods']
    MQ_file = MQ_file.fillna('')
    MQ_file['prot_tax_str']=['no species']*len(MQ_file)
    MQ_file['prot_desc'] = ['no description']*len(MQ_file)
    MQ_file['prot_seq']=['no seq']*len(MQ_file)
    MQ_file['pep_scan_title']=[title+'_'+str(nr) for nr,title in enumerate(MQ_file['pep_scan_title'].values)]
    add_mods = []
    change_position = []
    for m,mp in MQ_file[['mods','pep_var_mod_pos']].values:
        if m == 'Unmodified':
            m = ''
        m =str(m)
        if 'Hydroxyproline' in m:
            m=m.replace('Hydroxyproline', 'Oxidation (P)')
        if ',' in m:
            m=m.replace(',','; ')
        if 'Deamidation' in m:
            m=m.replace('Deamidation','Deamidated')
        if 'Glu->pyro-Glu' in m:
            m=m.replace('Glu->pyro-Glu','Glu->pyro-Glu (E)')
        if 'Gln->pyro-Glu' in m:
            m=m.replace('Gln->pyro-Glu','Gln->pyro-Glu (Q)')
        add_mods.append(m)
        temp = ''
        mp = mp.replace('_(','&')
        for t in mp.split('('):
            if ')' in t:
                t=t.split(')')
                for t2 in t:
                    if t2.lower()==t2:
                        temp+='&'
                        
                    else:
                        temp += t2
            else:
                temp+=t
        mp=temp
        temp = ''
        for i in range(0,len(mp)):
            if i == 0 and mp[i]=='_':
                temp += '0.'
            elif i == 0 and mp[i]=='&':
                temp += '1.'
            elif i==len(mp)-1 and mp[i]=='_':
                temp+='.0'
            elif mp[i]=='&':
                temp = temp[:-1]
                temp += '1'
            else:
                temp += '0'
        change_position.append(temp)
    MQ_file['pep_var_mod_pos_old']=MQ_file['pep_var_mod_pos'].values
    MQ_file['pep_var_mod_pos']=change_position
    
    MQ_file['pep_var_mod']=add_mods
    df = MQ_file
    #pep_var_mod aanpassen
    #add pep_var_mod_pos
    
    df_4_uni= df[['prot_tax_str','prot_desc','prot_seq','pep_seq','pep_var_mod','pep_var_mod_pos']]
    df2=MQ_file
    protein = {}
    unimod_db, unimod = do_unimod(path,df_4_uni['pep_var_mod'].values)
    ids = {}

    for p, a,m in unimod_db[['PTM','AA','mass']].values:
        if a=='N-term':
            a='!'
        elif a=='C-term':
            a='&'
        add = '?'
        while add+a in AA_codes.keys():
            add+='?'
        if a=='!' or a=='&':
            AA_codes[add+a]=float(m)
        else:
            AA_codes[add+a]=AA_codes[a]+float(m)
        ids[add+a]=p
    
    return df2, unimod_db, ids,AA_codes

def animals_from_db_input_mix(sequence_db):
    sp = []
    for sequence,name in sequence_db.items():
        if 'OS=' in name:
            anim = name.split('OS=')[-1]
            anim = anim.split(' OX=')[0]
            if anim not in sp:
                sp.append(anim)
        elif '[' in name:
            anim = name.split('[')[1]
            anim = anim.split(']')[0]
            if anim not in sp:
                sp.append(anim)
        elif '|' in name:
            anim = name.split('|')[1]
            if anim not in sp:
                sp.append(anim)
        else:
            print('{} has no species'.format(name))
    return sp

def recover_taxa(ncbi,t_df,recov,tax_to_lin):
    
    if recov == 'individual':
        df_temp = t_df[(t_df['Scientific name']==ncbi)]
        if len(df_temp)==0:
            df_temp =  t_df[(t_df['Common name']==ncbi) | (t_df['Synonyms'].str.contains(ncbi)) | (t_df['Other Names'].str.contains(ncbi))]
        if len(df_temp)>1:
            if ncbi in df_temp['Scientific name']:
                df_temp = df_temp[[df_temp['Scientific name']==ncbi]]
            else:
                df_temp = df_temp.iloc[[0]]
        elif len(df_temp)==0:
            return False
        lin = []
        species_add = [(list(df_temp['Rank'].values)[0],ncbi)]
        previous = ''
        for lineage,rank in df_temp[['Lineage','Rank']].values:
            for l in lineage.split(', '):
                if (l,previous) in tax_to_lin:
                    if tax_to_lin[(l,previous)] == '':
                        lin.append(('unranked',l)) 
                    else:
                        lin.append((tax_to_lin[(l,previous)],l)) 
                    previous = l
                else:
                    lin.append(('unknown',l))
                    previous = l
        return species_add+lin[::-1]
    else:
        print('Making higher taxonomy')
        all_org = []
        for k,v in recov.items():
            if k.isdigit()==False and 'unclassified' not in k and 'sp.' not in k:
                org = [el[1] for el in v if el[1] not in all_org and el[1]!=k]
                all_org = list(set(all_org + org))
        orders = []
        for k,v in recov.items():
            org = [el[1] for el in v if el[0]=='order']
            if len(org)>0:
                orders = list(set(orders+org))
        df_temp=t_df[(t_df['Scientific name'].isin(all_org)) | (t_df['Lineage'].str.contains('|'.join(orders), regex=True))]
        #need to find higher taxons
        #need to find missing species under higher taxa uptill order level
        other_missing_taxa = []
        for sp,lineage,rank in df_temp[['Scientific name','Lineage','Rank']].values:
            lin = []
            species_add = [(rank,sp)]
            previous = ''
            for l in lineage.split(', '):
                if l not in recov:
                    other_missing_taxa.append(l)
                if (l,previous) in tax_to_lin:
                    if tax_to_lin[(l,previous)] == '':
                        lin.append(('unranked',l))
                    else:
                        lin.append((tax_to_lin[(l,previous)],l)) 
                    previous = l
                else:
                    lin.append(('unknown',l))
                    previous = l
            recov[sp]=species_add+lin[::-1]
        df_temp=t_df[t_df['Scientific name'].isin(other_missing_taxa)]
        #need to find higher taxons
        #need to find missing species under higher taxa uptill order level
        for sp,lineage,rank in df_temp[['Scientific name','Lineage','Rank']].values:
            lin = []
            species_add = [(rank,sp)]
            previous = ''
            for l in lineage.split(', '):
                if (l,previous) in tax_to_lin:
                    if tax_to_lin[(l,previous)] == '':
                        lin.append(('unranked',l))
                    else:
                        lin.append((tax_to_lin[(l,previous)],l)) 
                    previous = l
                else:
                    lin.append(('unknown',l))
                    previous = l
            recov[sp]=species_add+lin[::-1]
    return recov

def animals_from_db_input(sequence_db, lim_t,demo,path):
    if lim_t == None:
        print('searching against all species in the database')
    else:
        print('Making selection of {} and random other species'.format(lim_t))
        lim_t = lim_t.split('/')
        lim_sp = []
        for x in lim_t:
            lim_sp.append(' '.join(x.split('_')))
        lim_t=lim_sp
    print('Loading taxonomy')
    taxonomy_loc = path/'MISC'
    taxa_file = []
    for file in taxonomy_loc.iterdir():
        if file.suffix == '.tsv':
            taxa_file.append((file, None))
    print(f'Opening taxonomy file: {taxa_file[0][0]}')
    taxonomy_df = pd.read_csv(taxa_file[0][0],sep = '\t')
    taxonomy_df = taxonomy_df.fillna('')
    tax_to_lin = {(sc_name,lin.split(', ')[-1]):rank for sc_name,lin,rank in taxonomy_df[['Scientific name','Lineage','Rank']].values if sc_name != ''}
    print('Taxonomy loading finished')
    input_animals = []
    skip_animals = []
    Class = {}
    taxonomy = {}
    species = []
    for sequence,name in sequence_db.items():
        if 'OS=' in name:
            anim = name.split('OS=')[-1]
            anim = anim.split(' OX=')[0]
        elif '[' in name:
            anim = name.split('[')[1]
            anim = anim.split(']')[0]
        elif '|' in name:
            anim = name.split('|')[1]
        else:
            print('{} has no species'.format(name))
            continue
        if anim not in species:
            species.append(anim)
    temp_df = taxonomy_df.copy()
    temp_df = temp_df[(temp_df['Scientific name'].str.contains('|'.join(species), regex=True)) | 
                      (temp_df['Common name'].str.contains('|'.join(species), regex=True)) |
                      (temp_df['Other Names'].str.contains('|'.join(species), regex=True)) |
                      (temp_df['Synonyms'].str.contains('|'.join(species), regex=True))]
    for anim in species:
        if lim_t == None or anim in ['Pseudomonas aeruginosa','Sus scrofa']:
            ncbi_animal = anim
            found = False
            input_animals.append(anim)
            while found == False:
                taxon = recover_taxa(ncbi_animal,temp_df,'individual',tax_to_lin)
                if taxon == False:
                    
                    if ' ' not in ncbi_animal:
                        found= True
                        print('Species {} has no taxonomy'.format(anim))
                        skip_animals.append(anim)
                        break
                    ncbi_animal = ' '.join(ncbi_animal.split(' ')[:-1])
                    continue
                found = True
                taxonomy[anim] = taxon
        else:
            ncbi_animal = anim
            found = False
            while found == False:
                taxon = recover_taxa(ncbi_animal,temp_df,'individual',tax_to_lin)
                if taxon == False:
                    
                    if ' ' not in ncbi_animal:
                        found= True
                        print('Species {} has no taxonomy'.format(anim))
                        skip_animals.append(anim)
                        break
                    ncbi_animal = ' '.join(ncbi_animal.split(' ')[:-1])
                    continue
                found = True
                taxonomy[anim]=taxon
                added_to_input = False
                added_to_random = False
                for tax in taxon:
                    for lt in lim_t:
                        if lt in tax:
                            added_to_input = True
                            print('Adding {} because it is within {}'.format(anim,lim_t))
                            input_animals.append(anim)
                            added_to_random = True
                    if 'class' in tax and added_to_random == False:
                        added_to_random = True
                        if tax[1] not in Class.keys():
                            Class[tax[1]]=[]
                        Class[tax[1]]=Class[tax[1]]+[anim]
                if added_to_input == False:
                    skip_animals.append(anim)
      
    if lim_t != None and demo==False:
        for k,v in Class.items():
            print(k)
            random.shuffle(v)
            random_animals =v[:10]
            input_animals = list(set(input_animals)|set(random_animals))
            skip_animals = list(set(skip_animals)^set(random_animals))
    taxonomy = recover_taxa(input_animals,taxonomy_df,taxonomy,tax_to_lin)
    skip_animals = list(set(skip_animals + [el for el in input_animals if el not in taxonomy.keys()]))
    input_animals = list(set([el for el in input_animals if el not in skip_animals]))
    return input_animals, skip_animals,taxonomy

def find_mass_matches(sequence, p_mass,pep,unimod_masses,AA_codes,uncertain):
    sequence +='&&&'
    if any(el in sequence for el in uncertain.keys()):
        ambig = True
    else:
        ambig = False
    start = 0
    max_mass = max(unimod_masses)
    min_mass = min(unimod_masses)
    end = len(sequence)
    possible = []
    temp = len(pep)-1
    tryptic = True if pep.endswith('R') or pep.endswith('K') else False
    while start+temp <= end-2:
        while tryptic==True and sequence[start+temp-1] != pep[-1]:#if tryptic look for the next tryptic pepitide in the sequence matching the mass
            temp += 1
            if temp-start>len(pep)+1:
                start += 1
                temp = len(pep)-1
            if start+temp >= end-1:
                break
        
        test_seq = sequence[start:start+temp]#stepwise window slide
        if len(set(test_seq)^set(pep))>3:
            start +=1
            temp = len(pep)-1
            continue
        if len(test_seq)>len(pep)+1:#if test peptide too large, reset slide
            start += 1
            temp = len(pep)-1
            continue
        unknown = False
        if ambig == True:
            if any(el in test_seq for el in uncertain.keys()):#see if there is a B,Z or X in the test_seq
                keep_seq = test_seq
                other_seq = ''.join([el for el in test_seq if el in uncertain.keys()])#extract the B,Z,X
                test_seq = ''.join([el for el in test_seq if el not in uncertain.keys()])#Retain amino acids
                unknown = True
        test = find_mass(test_seq,AA_codes)#test mass will be extact OR extact without BZX is present
        if unknown == True:#find all potential peptide candidates that have a mass match when BZX
            if len(other_seq)>2 or other_seq.count('X')>1:#Too ambiguous!
                start += 1
                temp = len(pep)-1
                continue
            missing_mass = p_mass-test#find ambiguous mass
            missing_mass = [missing_mass-num for num in unimod_masses]#include ptms
            other_seq = [uncertain[el] for el in other_seq]#potential AA masses
            poss_seqs = []
            for l in range(0,len(other_seq)):
                if len(poss_seqs)==0:
                    poss_seqs = [num for num in other_seq[l]]
                else:
                    poss_seqs = [num+el for num in poss_seqs for el in other_seq[l]]
            locx,locy = start,start+temp
            for x in set(poss_seqs):
                if any(num-0.015<=(test-p_mass+find_mass(x,AA_codes))<=num+0.015 for num in unimod_masses):
                    new_seq = ''.join([num if num in AA_codes.keys() else '!' for num in keep_seq])
                    add_seq = ''
                    count_x = 0
                    for m in new_seq:
                        if m=='!':
                            add_seq+=x[count_x]
                            count_x +=1
                        else:
                            add_seq+=m
                    possible.append((add_seq,(locx,locy)))
            temp += 1
        elif test-p_mass<min_mass or test-p_mass>max_mass or len(test_seq)<len(pep)-1:#extend test_seq if no more than length 2 difference, mass can be bigger because negative ptm mass
            temp +=1
        elif any(num-0.015<=test-p_mass<=num+0.015 for num in unimod_masses): #go if same mass or with a new ptm or without an existing one
            possible.append((test_seq,(start,start+temp)))
            start += 1
            temp = len(pep)-1
        else:#go on with the slide, reset
            start += 1
            temp = len(pep)-1
    return possible

def assign_pairs(index_to2, seq_real,check_seq):
    seq_real += '-'
    check_seq += '-'
    seq_real_new = []
    for i in index_to2:
        seq_real_new.append((seq_real[i],(i-1,i)))
        seq_real_new.append((seq_real[i],(i,i+1)))
        if i >1 and i<len(seq_real)-3:
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i+1]=='-':
                extra = 1
            if check_seq[i+1+extra] == '-':
                extra2 =1
            if check_seq[i+2+extra] == '-':
                extra3 =1
            seq_real_new.append((seq_real[i]+seq_real[i+1+extra],(i,i+1+extra+extra2)))#recht 1
            seq_real_new.append((seq_real[i]+seq_real[i+2+extra],(i,i+2+extra+extra3)))#rechts 2
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i-1]=='-':
                extra = 1
            if check_seq[i-1-extra]=='-':
                extra2 = 1
            if check_seq[i-2-extra]=='-':
                extra3 = 1
            seq_real_new.append((seq_real[i-1-extra]+seq_real[i],(i-extra-1-extra2,i))) #links 1
            seq_real_new.append((seq_real[i-2-extra]+seq_real[i],(i-extra-2-extra3,i))) #links 2
        elif i==0:
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i+1]=='-':
                extra = 1
            if check_seq[i+1+extra] == '-':
                extra2 =1
            if check_seq[i+2+extra] == '-':
                extra3 =1
            seq_real_new.append((seq_real[i]+seq_real[i+1+extra],(i,i+1+extra+extra2))) #rechts 1
            seq_real_new.append((seq_real[i]+seq_real[i+2+extra],(i,i+2+extra+extra3))) #rechts 2
        elif i == len(seq_real)-1:
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i-1]=='-':
                extra = 1
            if check_seq[i-1-extra]=='-':
                extra2 = 1
            if check_seq[i-2-extra]=='-':
                extra3 = 1
            seq_real_new.append((seq_real[i-1-extra]+seq_real[i],(i-extra-1-extra2,i))) #links 1
            seq_real_new.append((seq_real[i-2-extra]+seq_real[i],(i-extra-2-extra3,i))) #links 2
        elif i==1:
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i+1]=='-':
                extra = 1
            if check_seq[i+1+extra] == '-':
                extra2 =1
            if check_seq[i+2+extra] == '-':
                extra3 =1
            seq_real_new.append((seq_real[i]+seq_real[i+1+extra],(i,i+1+extra+extra2))) #rechts 1
            seq_real_new.append((seq_real[i]+seq_real[i+2+extra],(i,i+2+extra+extra3))) #rechts 2
            if seq_real[0]!='-':
                seq_real_new.append((seq_real[i-1]+seq_real[i],(i-1,i))) #links 1
        elif i == len(seq_real)-2:
            extra=0
            extra2 = 0
            extra3 = 0
            if seq_real[i-1]=='-':
                extra = 1
            if check_seq[i-1-extra]=='-':
                extra2 = 1
            if check_seq[i-2-extra]=='-':
                extra3 = 1
            seq_real_new.append((seq_real[i-1-extra]+seq_real[i],(i-extra-1-extra2,i))) #links 1
            seq_real_new.append((seq_real[i-2-extra]+seq_real[i],(i-extra-2-extra3,i))) #links 2
            extra=0
            if seq_real[-1]!='-':
                seq_real_new.append((seq_real[i]+seq_real[i+1+extra],(i,i+1+extra)))#recht 1
    return seq_real_new

def program(seq_db,seq_real,mass_matrix): #compare in-silico peptide with the found peptide
    #code below looks for where the differences are between the sequences, and includes adjecent amino acids to check aswel
    seq_db_new=[]

    index_to2 = [loc for loc in range(0,len(seq_db)) if seq_db[loc]=='-']
    index_to1 =  [loc for loc in range(0,len(seq_real)) if seq_real[loc]=='-']
    if len(index_to2)>4 or len(index_to1)>4: #more than 5 different locations is too much
        return False,[],[], []
    
    seq_real_new= assign_pairs(index_to2, seq_real,seq_db)
    seq_db_new= assign_pairs(index_to1, seq_db,seq_db)
    
    seq_1=[num for num,loc in seq_db_new]
    seq_2=[num for num,loc in seq_real_new]
    
    adaptation_db=[]
    adaptation_real=[]
    combo=[]
    
    for l, r,m in mass_matrix: #check for all diferences if they explain isobaric changes
        l1 = l.replace('|','')
        r1 = r.replace('|','')
        if '?' in l1:
            l1 = l1.replace('?','')
        if '?' in r1:
            r1 = r1.replace('?','')
        if l1 not in seq_1:
            continue
        if l1 in seq_1 and r1 in seq_2:
            adaptation_real.append(r1)#find isobaric
            adaptation_db.append(l1)
            combo.append((l,r))
    
    adaptation_real = [num for num in seq_real_new if num[0] in adaptation_real]
    adaptation_db = [num for num in seq_db_new if num[0] in adaptation_db]
    return True, adaptation_real, adaptation_db, combo

def do_alignment(s1,s2):#observed, db
    #align the sequences based on perfect matching, is quicker and better for our purposes than global alignment van de Bio package
    s1 = [num for num in s1]
    s2 = [num for num in s2]
    s1_align = ''
    s2_align = ''
    loc_s2 = 0
    loc_s1 = 0
    while loc_s2 < min(len(s2),len(s1)) and loc_s1 < min(len(s2),len(s1)):
        if s1[loc_s1]==s2[loc_s2]:
            s1_align += s1[loc_s1]
            s2_align += s2[loc_s2]
        elif loc_s1 < len(s1)-1 and loc_s2  < len(s2)-1:
            if s1[loc_s1+1] == s2[loc_s2]:
                s1_align += s1[loc_s1]+s1[loc_s1+1]
                s2_align += '-'+s2[loc_s2]
                loc_s1+=1
            elif s1[loc_s1] == s2[loc_s2+1]:
                s2_align += s2[loc_s2]+s2[loc_s2+1]
                s1_align += '-'+s1[loc_s1]
                loc_s2+=1
            else:
                s1_align += s1[loc_s1]+'-'
                s2_align += '-'+s2[loc_s2]   
        else:
            s1_align += '-'+s1[loc_s1]
            s2_align += s2[loc_s2]+'-'
        loc_s1 += 1
        loc_s2 += 1
    while loc_s1 != len(s1):
        s1_align += s1[loc_s1]
        s2_align += '-'
        loc_s1 +=1
    while loc_s2 != len(s2):
        s2_align += s2[loc_s2]
        s1_align += '-'
        loc_s2 +=1
    # remove_num = []
    # for i in range(0,len(s1_align)):
    #     if s1_align[i]=='-' and s2_align[i]=='-':
    #         remove_num.append(i)
    return (s1_align,s2_align)#(''.join([el for num,el in enumerate(s1_align) if num not in remove_num]), ''.join([el for num,el in enumerate(s2_align) if num not in remove_num]))

def locate_switches(adapt_observed,adapt_db, seq_observed, seq_db,combo):
    #of all possibilities found, check if the isobaric switch can occur at the location AND chose the smallest isobaric switch
    covering_seq = ''
    for i in range(0,len(seq_db)): #make one lin sequence that is like a stitched version of both sequences
        if seq_db[i]=='-':
            covering_seq+=seq_observed[i]
        else:
            covering_seq+=seq_db[i]
    adapt_observed=sorted(adapt_observed, key=lambda x:len(x[0])) #1 amino acid switches are preferred to multiple
    adapt_db = sorted(adapt_db, key=lambda x:len(x[0]))
    index_to1 =  [loc for loc in range(0,len(seq_observed)) if seq_observed[loc]=='-'] #all indexes that are gaps in seq_observed
    possible = []
    include_ptm=[]
    fixed_it = {}
    for i in index_to1: #find an explanation for each of the gaps
        fixed = False
        fixed_it[i]=False
        for n in adapt_db: #for each of the isobaric switches in adapt_db look if they fit the gap
            if i in n[1] and fixed == False:
                new_loc = [x for x in n[1] if x != i]
                for t in adapt_observed:
                    if new_loc[0] in t[1] and -1 not in t[1] and len(covering_seq) not in t[1]:
                        annot_check = [True if ''.join([x for x in num[0] if x.isalpha()])==n[0] and ''.join([x for x in num[1] if x.isalpha()])==t[0] else False for num in combo]
                        if True not in annot_check:
                            continue
                        annot =[num for loc_annot, num in enumerate(combo) if annot_check[loc_annot]==True]
                        if n[1]==t[1] and len(n[0])==1:#1 VS 1
                            test1 = n[0]+t[0]
                            test2= t[0]+n[0]
                            if test1[0]==covering_seq[new_loc[0]] and test1[1]==covering_seq[i]:
                                fixed=True
                                fixed_it[i]=True
                                possible.append('from '+annot[0][0]+' to '+annot[0][1])
                                if '?' in annot[0][0] or '?' in annot[0][1]:
                                    include_ptm.append((annot[0][0],annot[0][1], (new_loc[0],i)))
                                break
                            else:
                                 test2[0]==covering_seq[new_loc[0]] and test2[1]==covering_seq[i]   
                                 fixed = True
                                 fixed_it[i]=True
                                 possible.append('from '+annot[0][1]+' to '+annot[0][0])
                                 if '?' in annot[0][0] or '?' in annot[0][1]:
                                     include_ptm.append((annot[0][1],annot[0][0],(new_loc[0],i)))
                                 break
                        else:#>=1 VS >1
                            temp = ''
                            for z in range(0,len(covering_seq)):
                                if z in n[1]:
                                    if len(n[0])==1:
                                        search = 0
                                    else:
                                        search = n[1].index(z)
                                    temp+= n[0][search]
                                elif z in t[1]:
                                    if len(t[0])==1:
                                        search = 0
                                    else:
                                        search = t[1].index(z)
                                    temp += t[0][search]
                                else:
                                    temp += covering_seq[z]
                            if temp==covering_seq:
                                possible.append('from '+annot[0][1]+' to '+annot[0][0])
                                fixed = True
                                fixed_it[i]=True
                                if '?' in annot[0][0] or '?' in annot[0][1]:
                                    include_ptm.append((annot[0][1],annot[0][0], (new_loc[0],i)))
                                break
    if False in fixed_it.values() or len(possible)==0: #means that the sequences was not filled in properly
        return [],False,''
    #return the good annotations
    return possible, True, include_ptm

def find_ptm_location(ptm, seq,unimod_db,mascot_pos,ids,ptm2=[]):
    if len(mascot_pos)>0:
        mascot_pos=mascot_pos.split('.')[1]
    else:
        mascot_pos='0'*len(seq)
    loc_str = ''
    minus=[]
    extra_mass = 0
    if len(ptm2)>0:
        for i in ptm2:
            if '?' in i[1]:
                temp = seq[min(i[2])-2:max(i[2])]
                temp_lo = ''
                for z in i[1]:
                    if z=='?':
                        temp_lo += '?'
                    elif '?' in temp_lo:
                        lo = z
                        temp_lo += z
                        
                        if lo not in temp:
                            return False, False
                        found = temp.index(lo)
                        loc_str+=str(found+min(i[2])-1)+'|'+ids[temp_lo]+'|'
                        
                        extra_mass += float(unimod_db['mass'][(unimod_db['PTM']==ids[temp_lo])&(unimod_db['AA']==z)].values)
                        temp_lo = ''
            if '?' in i[0]:
                temp_lo = ''
                for z in i[0]:
                    if z == '?':
                        temp_lo += z
                    elif '?' in temp_lo:
                        temp_lo+=z
                        minus.append((ids[temp_lo],temp_lo))
                        temp_lo = ''
    if len(ptm)>0:
        ptm = ptm.split(';')
        for i in ptm:
            count = [x for x in i if x.isdigit()==True]
            if len(count)>0:
                count = int(count[0])
            else:
                count = 1
            temp = i.split(' (')
            for aa, p in unimod_db[['AA', 'PTM']].values:
                if p in temp[0] and aa in temp[1]:
                    if aa=='N-term':
                        loc_str = '0|'+p+'|'+loc_str
                    elif aa=='C-term':
                        loc_str = loc_str+'-1|'+p+'|'
                    count -= minus.count((p,aa))
                    if count <0:
                        return False, False
                    if count >0:
                        for locs, t in enumerate(seq):
                            if locs > len(mascot_pos)-1: #if there is a difference in length
                                count -= 1
                                continue
                            if t==aa and str(locs+1) not in loc_str.split('|') and count >0 and int(mascot_pos[locs])!=0:
                                loc_str += str(locs+1)+'|'+p+'|'
                                count -= 1
                                extra_mass += float(unimod_db['mass'][(unimod_db['PTM']==p)&(unimod_db['AA']==aa)].values)
    temp = loc_str[:-1].split('|')
    temps = sorted([(int(temp[i]),(temp[i+1])) for i in range(0,len(temp)-1,2) if temp[i]!='-1'], key=lambda x:x[0])
    temps = ['|'.join([str(num[0]),num[1]]) for num in temps]
    if '-1' in temp:
        temps.append('-1|'+temp[temp.index('-1')+1])
    return '|'.join(temps), extra_mass
   
def check_ptms_mascot(ad, ptm, ids):#location of ptm in mascot output? If mascot for example doesn't find a deamidation, it means that no deamidation isobaric switch can occur.
    checks = []
    amount ={}
    temp = ptm.split(';')
    if len(temp)>0:
        for i in temp:
            digit_temp = ''.join([num for num in temp if num.isdigit()==True])
            if len(digit_temp)>0:
                digit_temp=int(digit_temp)
            else:
                digit_temp = 1
            other = ''.join([num for num in temp if num.isdigit()==False])
            for k,v in ids.items():
                if ''.join([num for num in k if num != '?']) in other and v in other:
                    amount[k]=digit_temp
    for i in ad:
        i = i.split('to')[0]
        if '?' not in i:
            checks.append(True)
        else:
            temp = ''
            for t in i:
                if t == '?':
                    temp += t
                elif '?' in temp:
                    temp += t
                    test = ids[temp]
                    if temp not in amount:
                        checks.append(False)
                    elif test in ptm and amount[temp]!=0:
                        checks.append(True)
                        amount[temp]=amount[temp]-1
                    else:
                        checks.append(False)
                    temp = ''
    if False in checks:
        # if len(ptm)>0:
        #     print('found one that is not possible:',ad,ptm)
        # else:
        #     print('found one that is not possible:',ad,'no mascot ptm')
        return False
    return True

def ptm_mass(ptm, unimod_lookup):
    ptm_group = re.compile(r"(?:(\d+)\s+)?([\w\s\-]+)\s+\((\w+)\)")
    total_mass = 0
    ptm_insearch = {}
    if len(ptm)>0:
        for i in ptm.split(';'):
            i = i.strip()
            match = re.match(ptm_group, i)
            if match:
                count = int(match.group(1)) if match.group(1) else 1
                ptm_name = match.group(2).strip()
                aas = match.group(3).strip()
                for aa in aas:
                    key = (ptm_name, aa)
                    if key in unimod_lookup:
                        total_mass += count * unimod_lookup[key]
                        ptm_insearch['_'.join(list(key))]=1
                        break
    return total_mass,ptm_insearch

def massblast(db,ids,PTMs, peptides,unimod_db,mascot_pos,title,mass_matrix, animal,AA_codes,uncertain,unimod_masses,unimod_lookup):
    db = str(db)
    done= []
    final_result = []
    iso_pep = {}
    amount_exact_match = len([el for el in set(peptides) if el in db])
    if amount_exact_match>2:# any(el in db for el in peptides):#need at least 1 decent match attributed to the protein, only isoblast peptides are not trustworthy enough
        if any(el in str(db) for el in uncertain.keys()):
            ambig = True
        else:
            ambig=False
        
        for ilocs,p in enumerate(peptides):
            if len(p)>=len(db)-2 or len(p)<6 or any(el not in AA_codes for el in p):
                continue
            ptm_loc,extra_mass = find_ptm_location(PTMs[ilocs],p,unimod_db,mascot_pos[ilocs],ids)
            if (p,PTMs[ilocs],ptm_loc) in done:
                continue
            if p in db:#exact match
                location = re.finditer(p, str(db))
                locs = []
                for i in location:
                    locs.append(i.span())
                final_result.append((p,p,PTMs[ilocs],locs,'Original',ptm_loc,title[ilocs]))
                done.append((p,PTMs[ilocs],ptm_loc))
                continue
            if len(final_result) > 0 and any(p in num and PTMs[ilocs] in num for num in final_result):
                start_temp = [num for num in final_result if p in num and PTMs[ilocs] in num]
                if len(start_temp)>1:
                    start_temp=[start_temp[0]]
            #we need to do this step to include PSMs of peptides we already allocated
                if len(start_temp)>0:
                    for st in start_temp:
                        ptm_locs = st[5] if p!=st[1] else ptm_loc
                        final_result.append((p, st[1],st[2], st[3],st[4],ptm_locs,title[ilocs]))
                    done.append((p,PTMs[ilocs],ptm_loc))
                    continue
            done.append((p,PTMs[ilocs],ptm_loc))
            ptm_masses,PTMs_insearch = ptm_mass(PTMs[ilocs],unimod_lookup)
            #remake the peptide as it can have other uncertainty
            if ambig ==True:
                exact_match_with_uncertainty = ''
                for l in p:
                    u = ''
                    for k,v in uncertain.items():
                        if l in v:
                            u += k
                    if len(u)>0:
                        u = f"[{l}{u}]"#need for this because the B,Z,X can be in the database sequence
                    else:
                        u = l
                    exact_match_with_uncertainty += u
                #Add all the exact matches, also if X,B or Z in sequence. Only 1 uncertain allowed in output
                match = re.search(exact_match_with_uncertainty, str(db))
                if match and sum(c not in AA_codes for c in match.group()) <= 1:#match all exact matches en remind the location
                    match = match.group()
                    location = re.finditer(exact_match_with_uncertainty,str(db))
                    locs = []
                    for i in location:
                        locs.append(i.span())
                    if any(c not in AA_codes for c in match):
                        label = 'With uncertainty'
                        group_used = p
                        final_result.append((p,group_used,PTMs[ilocs],locs,label,ptm_loc,title[ilocs]))
                        continue
            #start isoblast
            #if no exact match, than start isobaric switches
            if amount_exact_match<=8 and len(db)>=1000:
                continue
            p_mass: float= find_mass(p,AA_codes)+ptm_masses
            if (p,p_mass) in iso_pep:#peptidoforms of same mass will match to same possible peptide matches, no need to find them again
                possible= iso_pep[(p,p_mass)]
            else:
                possible = find_mass_matches(db, p_mass,p,unimod_masses,AA_codes,uncertain)
                iso_pep[(p,p_mass)]=possible
            for test_seq,location in possible:
                seq1 = p #original sequence
                seq2 = test_seq#database sequence
                alignment = do_alignment(seq1, seq2)
                alteration1 = alignment[1]
                alteration2 = alignment[0]
                if alteration1.count('-')<=3:# three isobaric mistakes allowed, else too much possibilities, and overall almost never correct 
                    test_out, adapt_real, adapt_db,combo = program(alteration1, alteration2, mass_matrix)
                    if not test_out:
                        continue
                    adapt, addition_ok,ptms_loc = locate_switches(adapt_real,adapt_db,alteration2,alteration1,combo)
                    if not addition_ok:
                        continue
                    addition_ok = check_ptms_mascot(adapt,PTMs[ilocs],ids)
                    if not addition_ok:
                        continue
                            
                    adapt.append(PTMs[ilocs])
                    ptm_loc,extra_mass = find_ptm_location(PTMs[ilocs],test_seq,unimod_db,mascot_pos[ilocs],ids,ptms_loc)
                    if ptm_loc is not False and abs(p_mass-(find_mass(str(test_seq),AA_codes)+extra_mass))<0.015:
                        final_result.append((p,str(test_seq),adapt,[location],'isobaric',ptm_loc,title[ilocs]))
    return final_result

def thread_worker(sequence, name, df2, unimod_db, mass_matrix, ids,
              check_all_psms, AA_codes, uncertain,
              unimod_masses, unimod_lookup):
    anim = isinput_species(name)
    future = massblast(sequence,ids,df2['pep_var_mod'].values,df2['pep_seq'].values,unimod_db,df2['pep_var_mod_pos'].values,df2['pep_scan_title'].values, mass_matrix, anim,AA_codes,uncertain,unimod_masses,unimod_lookup)
    if len(future) > 0:
        data = [(anim,name,*res) for res in future]
    else:
        data = [False]
    return data

def calc_distance_collagen(sequence_db, df):
    locs_dict = {}
    for k,seq_name in sequence_db.items():
        if seq_name not in df['protein'].values:
            continue
        all_locs = []
        for l in df['location'][df['protein']==seq_name].values:
            for n in l:
                for i in range(n[0],n[1]):
                    all_locs.append(i)
        all_locs = sorted(list(set(all_locs)))
        locs_dict[seq_name]=len(all_locs) #calculate coverage
    
    columns = []
    for animal in df['animal'].values:
        if animal not in columns:
            columns.append(animal)
            
    index = []
    z = []
    for l, n in df[['protein','animal']].values:
        if l not in index:
            temp = [0]*len(columns)
            loc = columns.index(n)
            temp[loc]=locs_dict[l] #should be unique, so easiest to do like this
            z.append(temp)
            index.append(l)
    
    return index, locs_dict, columns, z

def thread_align(i,seq,calculator,seq_keep):
    if i==seq:
        return (0,seq_keep)
    i = Seq(str(i).replace('&',''))
    seq = Seq(str(seq).replace('&',''))
    try:
        aligner = Align.PairwiseAligner()
        align = aligner.align(i,seq)
        align = next(align)
        # align = list(align)[0]
        a=SeqRecord(align[0].replace('-','Z'),id='a')
        b=SeqRecord(align[1].replace('-','B'),id='b')
        align = MultipleSeqAlignment([a, b])
        dm = calculator.get_distance(align)
        return (dm.matrix[1][0], seq_keep)
    except:
        return (10,seq_keep) 
    

def thread_worker2(process_executor,x,m,c,s_db):
    i = x[0]
    seq = x[1]
    future = process_executor.submit(thread_align, i,seq,m,c,s_db)
    
    out = future.result() 
    return out

def make_plots_coverage_per_animal(df, sequence_db, all_animals, data_array,labels,index, locs_dict, columns, z,names, file_name,sample_path,taxons):    
    df_plot=pd.DataFrame(np.array(z),index=index,columns=columns)
    print('calculating taxonomy tree')
           
    taxon_col = columns
    species_in = columns
    c_order = {}
    for n in species_in:
        t = taxons[n]
        count = 0
        for i in t:
            if 'order' in i:
                if count in c_order:
                    c_order[count] =  c_order[count]+1
                else:
                    c_order[count] =1
                break
            count += 1
    #find_most_common location of genus
    c_genus = {}
    for n in species_in:
        t = taxons[n]
        count = 0
        for i in t:
            if 'genus' in i:
                if count in c_genus:
                    c_genus[count] =  c_genus[count]+1
                else:
                    c_genus[count] =1
                break
            count += 1
    #find_most common location of LCA
    in_names = species_in
    lca_all = find_LCA(taxons,in_names)
    c_lca = {}
    for n in species_in:
        t = taxons[n]
        count = 0
        for i in t:
            if lca_all in i:
                if count in c_lca:
                    c_lca[count] =  c_lca[count]+1
                else:
                    c_lca[count] =1
                break
            count += 1
    golca = [[k for k,v in c_genus.items() if v==max(c_genus.values())][0],
             [k for k,v in c_order.items() if v==max(c_order.values())][0],
             [k for k,v in c_lca.items() if v==max(c_lca.values())][0]]
    
    transform_taxa = transform_taxons({},taxons,species_in,golca,lca_all)
    
    tree_matrix = []
    for n in in_names:
        begin = taxons[n]
        tree_row = []
        for n2 in in_names:
            if n==n2:
                tree_row.append(0)
                continue
            lca = find_LCA(taxons,[n,n2])
            t = taxons[n2]
            
            dist_n_lca=transform_taxa[n][[el[1] for el in begin].index(lca)]
            
            dist_n2_lca = transform_taxa[n2][[el[1] for el in t].index(lca)]
            tree_row.append(dist_n_lca+dist_n2_lca)
        tree_matrix.append(tree_row)
    taxon_distance = tree_matrix
    print('end taxonomy tree')
    taxon_distance = np.array(taxon_distance)
    if len(taxon_distance)<2:
        return taxon_distance, [], [],[],taxon_col
    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(data_array, orientation='bottom')
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'
    
    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(taxon_distance, orientation='right',labels=taxon_col)
    dendro_side2 = ff.create_dendrogram(taxon_distance, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
    
    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)
    
    # Create Heatmap
    dendro_leaves = fig['layout']['xaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    dendro_leaves2 = dendro_side2['layout']['yaxis']['ticktext']
    dendro_leaves2 = list(map(int, dendro_leaves2))
    

    heat_data = df_plot.T.values
    heat_data = heat_data[dendro_leaves2,:]
    heat_data = heat_data[:,dendro_leaves]
    
    heat = [go.Heatmap(
        x = dendro_leaves,
        y = dendro_leaves2,
        z = heat_data,
        text=[[num]*len(dendro_leaves) for num in np.array(taxon_col)[dendro_leaves2]],
        colorscale='Hot',
        type='heatmap'
    )]
    
    heat[0]['x'] = fig['layout']['xaxis']['tickvals']
    heat[0]['y'] = dendro_side['layout']['yaxis']['tickvals']
    
    for data in heat:
        fig.add_trace(data)
    
    fig.update_layout({'width':1500, 'height':2000,
                             'showlegend':False, 'autosize':True
                             })
    # Edit xaxis

    fig.update_layout(xaxis={'domain': [.3, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'ticks':"",
                                      'ticktext':np.array(labels)[dendro_leaves]})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .3],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':"",
                                        })
    
    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .7],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks': "",
                            })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.7, 1],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':"",
                                        
                                        })
    file_name_plot = 'Heatmap_'+file_name+'.html'
    fig.update_layout(title = 'Heatmap: '+file_name_plot)
    try:    
        fig.write_html(path /'Output_Classicol'/sample_path/file_name_plot)
    except:
        fig.write_html(path /'Output_Classicol'/sample_path/'Heatmap.html')
    # plotly.offline.plot(fig)
    
    return taxon_distance, df_plot.T, heat_data,names,taxon_col

def clustering(X,names):
    branches = {}
    X_test = linkage(X,'ward')
    i = 0
    one_cluster = False
    while one_cluster == False:
        cluster = list(fcluster(X_test,t=i,criterion='distance'))
        if len(set(cluster))==1:
            one_cluster = True
        if cluster not in branches.values():
            branches[i] = cluster
        i += 0.1
    tree = {}
    for i in sorted(list(branches.keys()))[::-1]:
        cluster_tree = {t:[] for t in branches[i]}
        for l,t in enumerate(branches[i]):
            cluster_tree[t] = sorted(cluster_tree[t]+[names[l]])
        tree[i]=cluster_tree
    return tree

def knapzak(all_input):
    X = all_input[0]
    names = all_input[1]
    Y=all_input[2]
    animals=all_input[3]
    heat = all_input[4]
    tree_sequences = clustering(X,names)
    tree_seqs = sorted([(i,t) for i,t in tree_sequences.items()], key=lambda x:x[0])

    tree_animals = clustering(Y,animals)
    tree_animals = sorted([list(t.values()) for i,t in tree_animals.items()], key=lambda x:x[0])
    t_ani = []
    for element in tree_animals:
        t_ani = t_ani+element
    # tree_animals = sorted([(i,t) for i,t in tree_animals.items()], key=lambda x:x[0])
    
    output = []
    found = False
    save_level=[]
    for level in tree_seqs[::-1]:
        temp_out_count = save_level
        level = level[1]
        if len(level) == len(heat.columns) or found==True:#we are only interested in multiple sequence clusters
            continue
        temp_out = []
        temp_out_to_check = []
        temp_out_count2 = []
        for cluster in level:
            indiv_cluster = level[cluster]
            temp = heat[indiv_cluster]
            temp = temp.loc[~(temp==0).all(axis=1)]
            if len(temp) == len(heat.index):
                continue
            temp_ani = sorted(list(temp.index.values))
            temp_out_count2.append(temp_ani)
            if temp_ani in t_ani:#are they evolutionary close?
                temp_out.append((temp_ani,indiv_cluster))
                temp_out_to_check.append(temp_ani)
        save_level = temp_out_count2 #save only former level
        for i in temp_out_to_check:
            
            count = 0
            for x in temp_out_count+save_level:
                if i == x:
                    count += 1
                elif len(set(i)|set(x))>max(len(x),len(i)):
                    count += 1
            if count == len(temp_out_count+save_level):
                for y in temp_out:
                    if i==y[0] and y[0] not in output:
                        output.append(y[0])
                        found = True
    return output

def filter_unique_clusters(df_distance):
    unique = []
    for i in df_distance.index:
        if min(df_distance.loc[i].values[df_distance.loc[i].values >0]) > 1:
            unique.append((i,i,True))
    return unique

def find_animals(X,names,Y,animals,heat):
    df_distance = pd.DataFrame(np.array(X),columns=names,index=names)
    df_distance_taxon = pd.DataFrame(np.array(Y),columns=animals,index=animals)
    results = []
    end = False
    filter_unique = filter_unique_clusters(df_distance)
    results = results+filter_unique
    f_u = []#num[0] for num in filter_unique]
    heat = heat.drop(axis=1, labels=f_u)
    heat = heat.loc[~(heat==0).all(axis=1)]
    animals = [num for num in animals if num in heat.index]
    names = [num for num in names if num not in f_u]
    X = df_distance[df_distance.index.isin(heat.columns)]
    X = X[heat.columns].values
    Y = df_distance_taxon[df_distance_taxon.index.isin(heat.index)]
    Y = Y[heat.index].values
    
    division = [[X,names,Y,animals,heat]]
    while end == False:
        divide = []
        for div in division:
            if len(div[0])==0:
                end=True
                continue
            
            end = True
            
            results.append((list(div[3]),list(div[1]),False))
        division=divide
    return results

def new_way(all_animals, all_sequences,df_heatmap_values,df_distance,df_distance_taxon,dfs):
   
    if len(all_animals)==1: #if no animals anymore, or in other words we reached the species level, quit
        return [all_animals]
    temp = df_heatmap_values[all_sequences] #temporary dataframe with only peptides of interest
    temp_d = df_distance[all_sequences] #clustering of the dataframe peptides
    temp_d = temp_d.loc[all_sequences] 
    temp_t = df_distance_taxon[all_animals] #taxonomic tree of animals of interest
    temp_t = temp_t.loc[all_animals]
    tree_sequences = clustering(temp_d,temp_d.columns) #recluster the collagens so noo influence of sequences of non-interest
    tree_animals = clustering(temp_t,temp_t.columns)#recluster animals, to quickly find the next branching point
    
    t_ani = []
    for element in sorted(tree_animals.keys())[::-1]: #take next split of the tree
        if len(tree_animals[element])!=1: #the first one will contain all animals, we want the next level
            for value in tree_animals[element].values():
                t_ani.append(value)
            break
    if len(t_ani)==0: #if no animals anymore, or in other words we reached the species level, quit
        return [all_animals]
    
    if len(t_ani)>2:#branched too much, so split at 1 vs rest, with the 1 being the most distant given the found peptides
        index_minidist = []
        minidist = []
        for a in t_ani:
            index_minidist = index_minidist + a
            temp1 = set(dfs['found_match'][dfs['animal'].isin(a)].values)
            temp2 = set(dfs['found_match'][dfs['animal'].isin(a)==False].values)
            diff = len(temp1^temp2)
            minidist.append(diff)
        maximal = max(minidist)
        animal_branch1 = [index_minidist[minidist.index(maximal)]]
        animal_branch2 = list(set(animal_branch1)^set(index_minidist))
        print('split {} to {}'.format(t_ani,[animal_branch1,animal_branch2]))
        t_ani = [animal_branch1,animal_branch2]
        
    #check if there is a difference in types of proteins found
    leave_animal = {}
    for r in t_ani: #per branch of animals
        done = False
        allow = 0
        while done == False:
            allow += 1
            levels = sorted(list(tree_sequences.keys()))[::-1]
            for level_i in levels:
                cluster = tree_sequences[level_i]
                concat_cluster = []
                for x in cluster.values():
                    done = dfs[['animal','protein']][dfs['protein'].isin(x)].values
                    done = {p:a for a,p in done}
                    conc = ''
                    for p in done.keys():
                        if list(done.values()).count(done[p])>1 and p not in conc:
                            temps = dfs[dfs['animal']==done[p]]
                            d_check= set()
                            for pr in temps['protein'].values:
                                d_check = d_check^set(temps['found_match'][temps['protein']==pr].values)
                                if len(d_check) != 0:
                                    conc += pr
                        else:
                            conc += p
                    concat_cluster.append(conc)
                animal_check = True
                for a in r:
                    separate_clusters = [True if el.count(a)<=allow else False for el in concat_cluster]
                    if False in separate_clusters:
                        animal_check = False
                        break
                if animal_check == True:
                    done = True #at this level all sequences in a different cluster
                    leave_animal[tuple(r)]=level_i
                    break 
    if len(leave_animal.values()) == 0:
        minimal_separation = 0
    else:
        minimal_separation = min(list(leave_animal.values())) #now we know the location in the collagen tree where all can be separated
    #Per cluster look at the difference at peptide level 
    animals_A = dfs[(dfs['animal'].isin(t_ani[0]))]
    pep_A = set(animals_A['found_match'].values)
    animals_B = dfs[(dfs['animal'].isin(t_ani[1]))]
    pep_B = set(animals_B['found_match'].values)
    output = {}
    if len(pep_A^pep_B)>0:
        for i in tree_sequences[minimal_separation].values():
            animals_A = dfs[(dfs['animal'].isin(t_ani[0]))&(dfs['protein'].isin(i))]
            pep_A = set(animals_A['found_match'].values)
            animals_B = dfs[(dfs['animal'].isin(t_ani[1]))&(dfs['protein'].isin(i))]
            if (len(animals_B)==0 or len(animals_A)==0) and minimal_separation!=0:#protein missingness included for each taxonomic level, if a group is absent than it cannot be used to distinct based on that protein group
                output[tuple(i)]=3
                continue
            pep_B = set(animals_B['found_match'].values)
            if len(pep_B^pep_A)==0:#the same
                output[tuple(i)]=0
            elif pep_A.issubset(pep_B) == True:# a is subset from b
                output[tuple(i)]=1
            elif pep_B.issubset(pep_A) == True:# b is subset from a
                output[tuple(i)]=2
            else:
                output[tuple(i)]=3 #both have unique sequences
    if len(output)==0:
        output['all'] = 0
    #1 means go on with animals from B
    #0 and 3 means go on with both (MIX)
    #2 means go on with animals from A  
    #4 means no difference in peptides overall
    score_out = []
    A_done = False
    B_done = False
    if 1 in output.values():
        print('do {}'.format(t_ani[1]))
        B_done = True
        new_names = temp.loc[t_ani[1]]
        new_names = new_names.loc[:, (new_names != 0).any(axis=0)]
        new_names = list(new_names.columns)

        score = new_way(list(t_ani[1]), new_names,df_heatmap_values,df_distance,df_distance_taxon,dfs)
        score_out = score_out + score
    if 2 in output.values():
        print('do {}'.format(t_ani[0]))
        A_done = True
        new_names = temp.loc[t_ani[0]]
        new_names = new_names.loc[:, (new_names != 0).any(axis=0)]
        new_names = list(new_names.columns)

        score = new_way(list(t_ani[0]), new_names,df_heatmap_values,df_distance,df_distance_taxon,dfs)
        score_out = score_out + score
    mix_done = False
    if 3 in output.values():
        
        mix_done = True
        for t in t_ani:
            if t == t_ani[1] and B_done == True: #not 2 times the same analysis
                continue
            elif t == t_ani[0] and A_done == True: #not 2 times the same analysis
                continue
            new_names = temp.loc[t]
            new_names = new_names.loc[:, (new_names != 0).any(axis=0)]
            new_names = list(new_names.columns)

            score = new_way(list(t), new_names,df_heatmap_values,df_distance,df_distance_taxon,dfs)
            score_out = score_out + score    
    if mix_done == False and A_done == False and B_done == False and 0 in output.values():#we need to stop this branch here
        print('Found a taxon limit')
        score_out = score_out + [all_animals]
    return score_out

def find_animals2_1(distance,names,taxon_distance,animals,df_heatmap_values,output,dfs):
    output2 = []
    df_distance = distance 
    df_distance_taxon = taxon_distance 
    o1 = []
    o2=[]
    for o in output:
        if True in o:
            continue
        o1 = o1+o[0]
        o2=o2+o[1]
    output = [(o1,o2)]
    for o in output:
        if len(o)>0:
            all_animals = o[0]
            all_sequences = o[1]
            print('scoring ...')
            score = new_way(all_animals,all_sequences,df_heatmap_values,df_distance,df_distance_taxon,dfs)
            #reiterate the score so non of the outcomes are subsets or the others
            keep = {}
            recall = {}
            print('potential candidates are {}'.format(score))
            print('reiteration ...')
            for i in score:
                contain_i = set(dfs['found_match'][dfs['animal'].isin(i)].values)
                re = list(set(dfs['found_match'][dfs['animal'].isin(i)].values))
                keep[tuple(i)]=True
                for t in score:
                    if t==i:
                        continue
                    contain_t = set(dfs['found_match'][dfs['animal'].isin(t)].values)
                    if contain_i.issubset(contain_t)==True and len(contain_i^contain_t)>2:#no accidental subsetters
                        keep[tuple(i)]=False
                    else:
                        re = [el for el in re if el not in contain_t]
                if keep[tuple(i)]==True:
                    recall[tuple(i)]=re
            score = [el for el in keep.keys() if keep[el]==True]
            #adapt for subsetters over multi species
            output_score = []
            for i in score:
                keep = {}
                if len(i)>1:
                    for x in i:
                        x = [x]
                        keep[tuple(x)]=True
                        contain_x = set(dfs['found_match'][dfs['animal'].isin(x)].values)
                        re = list(set(recall[i]))
                        for y in i:
                            y = [y]
                            if x==y:
                                continue
                            contain_y = set(dfs['found_match'][dfs['animal'].isin(y)].values)
                            if contain_x.issubset(contain_y)==True and len(contain_x^contain_y)>2:#no accidental subsetters
                                keep[tuple(x)]=False
                            else:
                                re = [el for el in re if el not in contain_y]
                    new_tuple = tuple([el[0] for el in keep.keys() if keep[el]!=False])
                    if len(new_tuple)>0:
                        recall[tuple(new_tuple)]=re
                    output_score.append(new_tuple)
                else:
                    output_score.append(i)
            score = output_score
            total_amount = len(set(dfs['found_match'].values))
            delete = []
            for i in score:
                contain_i = set(dfs['found_match'][dfs['animal'].isin(i)].values)
                if len(contain_i)<total_amount*0.1:
                    print('deleting {} because it comprises lower than 10% of total peptides'.format(i))
                    delete.append(i)
            score = [num for num in score if num not in delete]
            
            combos = []
            for el in score:
                for num in score:
                    if el != num and tuple(sorted(list(el+num))) not in combos:
                        combos.append(tuple(sorted(list(el+num))))
            print('Catching sneaky ones ...')
            keep = {}
            deleted = []
            all_possible = [el for el in set(dfs['animal'].values) if any(el in element for element in score)]
            dfs = dfs[dfs['animal'].isin(all_possible)]#filter out all the discarded ones
            print('{} taxa left'.format(len(set(dfs['animal'].values))))
            for i in score:
                keep[tuple(i)]=True
                contain_i = set(dfs['found_match'][dfs['animal'].isin(i)].values)
                contain_global = set(dfs['found_match'][dfs['animal'].isin(i)==False].values)
                diff_global = contain_i^contain_global
                if len(diff_global&contain_i)>0:#unique peptide so keep, no need to check against all combos
                    continue
                for x in combos:
                    if any(num in deleted for num in x):#continue if it contains a discarded species
                        continue
                    if any(num in x for num in i): #if a unique peptide, continue#len([num for num in i if num not in list(x)])==0:
                        continue
                    contain_x = set(dfs['found_match'][dfs['animal'].isin(x)].values)
                    if (contain_i.issubset(contain_x) and (len(contain_i)<len(contain_x)*0.85 or sorted(list(contain_x))==sorted(list(contain_i)))):
                        delete = True
                        for animal_check in x:
                            check_temp = set(dfs['found_match'][dfs['animal'].isin([animal_check])].values)
                            if check_temp.issubset(contain_i) or len(contain_i^check_temp)<=2:
                                delete = False
                                break
                        if delete==True:
                            combos = [el for el in combos if any(element in el for element in i)==False]
                            deleted.append(i)
                            keep[tuple(i)]=False
                            print('deleting',i)
                            dfs = dfs[dfs['animal'].isin(i)==False]#make the database shorter so it takes less time, deleting the ones that are discarded
                            break
            output_score=[el for el in keep.keys() if keep[el]!=False]
            if len(output_score)>0:
                output2.append(output_score)
            else:
                output2.append(score)
    return output2, recall,dfs

def make_sunburst(dfs,all_animals,df_output,file_name,ip_animals,df_plot,path,sample_path,taxonomy):
    #sunburst plot to retrace the track 
    all_taxonomy = {k:v for k,v in taxonomy.items() if k in list(set(list(set(df_output['animal'].values))+ip_animals))}
    labels = []
    values = []
    parents = []
    done = []
    values_iso = []
    values_uni = []
    braycurtis_dist = []
    all_peptides = list(set(list(dfs['found_match'].values)))#mascot peptide
    array1 = [1 if pep in all_peptides else 0 for pep in all_peptides]
    seq_concat = '_'.join(all_peptides)
    weights_bc = [1/seq_concat.count(val) for val in all_peptides]
    
    combo_taxon = []
    for key,val in all_taxonomy.items():
        if key in set(df_output['animal'].values):
            combo_taxon = combo_taxon+val
    combo_taxon=set(combo_taxon)
    
    for animal in sorted(list(set(df_output['animal'].values)), key=lambda x:len(x))[::-1]:
        if animal in labels:#means that the subspecies has already been assigned, otherwise it will go in double
            continue
        a_taxon = all_taxonomy[animal]
        labels.append(animal)
        values.append(len(set(dfs['found_match'][dfs['animal']==animal].values)))
        values_iso.append(len(set(dfs['found_match'][(dfs['animal']==animal) & (dfs['type'].str.contains('isobaric')==True)].values)))
        values_uni.append(len(set(dfs['found_match'][(dfs['animal']==animal)].values)^set(dfs['found_match'][dfs['animal']!=animal].values)&set(dfs['found_match'][(dfs['animal']==animal)].values)))
        temp_peptides = set(dfs['mascot_peptide'][dfs['animal']==animal].values)
        array2 = [1 if pep in temp_peptides else 0 for pep in all_peptides]
        score = braycurtis(array1, array2,w=weights_bc)
        braycurtis_dist.append(1-score)#columns need to be the same
        ai = []
        stop = False
        previous_taxon = ''
        for t in a_taxon:
            if stop == True:
                break
            if previous_taxon == t[1]:
                continue
            else:
                previous_taxon = t[1]
            if (t in combo_taxon or t[1] in parents) and t[1]!=animal:
                animals_involved = [key for key,val in all_taxonomy.items() if t in val and key in dfs['animal'].values]
                if sorted(animals_involved) == sorted(ai):
                    if len(animals_involved) == len(set(dfs['animal'])):
                        parents.append('')
                        stop = True
                        continue
                ai = animals_involved
                add = len(set(dfs['found_match'][dfs['animal'].isin(animals_involved)].values))
                add_iso = len(set(dfs['found_match'][(dfs['animal'].isin(animals_involved)) & (dfs['type'].str.contains('isobaric')==True)].values))
                
                add_uni = len(set(dfs['found_match'][(dfs['animal'].isin(animals_involved))].values)^set(dfs['found_match'][(dfs['animal'].isin(animals_involved))==False].values)&set(dfs['found_match'][(dfs['animal'].isin(animals_involved))].values))
                if t[1] in set(df_output['animal'].values):
                    if t[1] in animal and t[0]=='species' and t[1]!=animal:
                        parents.append(t[1])
                        break
                if [t[1]] in done and 'species' not in t:
                    parents.append(t[1])
                    break
                done.append([t[1]])
                parents.append(t[1])
                labels.append(t[1])
                values.append(add) 
                values_iso.append(add_iso)
                values_uni.append(add_uni)
                temp_peptides = set(dfs['found_match'][dfs['animal'].isin(animals_involved)].values)#mascot peptide
                #we want only the shared peptides/level. At species level the uniqueness does count
                unique_peptides = []
                for peptide_test in temp_peptides:
                    test_pep = set(dfs['animal'][(dfs['found_match']==peptide_test)&(dfs['animal'].isin(animals_involved))].values)#mascot_peptide
                    if len(test_pep)!= len(animals_involved):
                        unique_peptides.append(peptide_test)
                temp_peptides = temp_peptides^set(unique_peptides)
                array2 = [1 if pep in temp_peptides else 0 for pep in all_peptides]
                score = braycurtis(array1, array2,w=weights_bc)
                braycurtis_dist.append(1-score)#columns need to be the same
        if len(parents)<len(values):
            parents.append('')
    temporary_dataframe = pd.DataFrame()
    temporary_dataframe['values']=values
    temporary_dataframe['parents']=parents
    temporary_dataframe['labels']=labels
    temporary_dataframe['values_iso']=values_iso
    temporary_dataframe['braycurtis_dist']=braycurtis_dist
    values = []
    parents = []
    labels = []
    values_iso = []
    braycurtis_dist = []
    done = []
    for x in temporary_dataframe.values:
        if (x[1],x[2]) not in done and x[2] not in labels:#account for species and species in there because of subspecies
            values.append(x[0])
            parents.append(x[1])
            labels.append(x[2])
            values_iso.append(x[3])
            braycurtis_dist.append(x[4])
            done.append((x[1],x[2]))
    
    
    fig = go.Figure()
    fig.add_trace(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="remainder",
        meta = values_iso,
        marker=dict(
        colors=braycurtis_dist,
        coloraxis='coloraxis1',
        colorscale='Jet',
        showscale=True,
        cmid=0.5),
    hovertemplate='<b>%{label} </b> <br> Taxon score: %{color:.3f}<br> Total peptide count: %{value:.0f} <br> isoBLAST peptidoforms: %{meta:.0f}'))
    fig.update_layout(title = 'Score of output species to sample '+file_name)
    fig.update_layout(autosize=True,margin = dict(t=30, l=0, r=0, b=0))
    fig.update_layout(coloraxis_colorbar_title='Score')
    try:
        file_name_plot = 'sunburst_'+file_name+'.html'
        fig.write_html(path/'Output_Classicol'/sample_path/file_name_plot)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'sunburst.html')
    # plotly.offline.plot(fig)
    # time.sleep(2)
    f_out = make_output_file(path, df_plot, df_output,file_name, braycurtis_dist, labels,sample_path,all_taxonomy)
    make_sunburst_with_missing2(dfs,ip_animals,list(set(df_output['animal'].values)),path,sample_path,taxonomy)
    return braycurtis_dist, labels, f_out

def make_output_file(path,df_plot, df_output, file_name, bc, l,sample_path,all_taxonomy):
    print('generating output file')
    ranking_taxon = []
    df_distance = pd.DataFrame(np.array(bc).reshape(1,-1),columns=l)
    score_to_taxon = {}
    for animal in set(list(df_output['animal'].values)):
        score = list(df_distance[animal].values)[0]
        if score in score_to_taxon:
            score_to_taxon[score] = [animal]+score_to_taxon[score]
        else:
            score_to_taxon[score] = [animal]
    print(score_to_taxon)
    final_taxon = {}
    for sc,stt in score_to_taxon.items():
        for i in stt:
            taxon = all_taxonomy[stt[0]]
            for t in stt:
                taxon = set(taxon)&set(all_taxonomy[t])
            for t in all_taxonomy[stt[0]]:
                if t in taxon:
                    for x in stt:
                        final_taxon[x]=t[1]
                    break
    taxon = []
    for animal in df_output['animal'].values:
        taxon.append(final_taxon[animal])
    df_output['taxon']=taxon
    
    for taxon in set(df_output['taxon'].values):
        animals_in_taxon_out = set(df_output['animal'][df_output['taxon']==taxon].values)
        score = 0
        for a in animals_in_taxon_out:
            score += df_distance[a].values
        score = score/len(animals_in_taxon_out)
        ranking_taxon.append([taxon,score])
    ranking_taxon = sorted(ranking_taxon, key=lambda x:x[1])[::-1]
    name_file = 'ZooMSMS_results_'+file_name+'.csv'
    try:
        with open(path/'Output_Classicol'/sample_path/name_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['ZooMSMS analysis output (~: isoBLAST match, *: Unique PSM)'])
            writer.writerow(['Isoblast analysis revealed '+ str(len(set(df_output['found_match'].values)))+' unique peptide matches with the database'])
            for taxon,size in ranking_taxon:
                writer.writerow(' ')
                writer.writerow(['Taxonomic match: '+taxon+' Average Score= '+str(size)])
                animals = sorted(list(set(df_output['animal'][df_output['taxon']==taxon].values)))
                writer.writerow(['Containing these species:'+', '.join(animals)])
                for a in animals:
                    proteins = list(set(df_output['protein'][df_output['animal']==a].values))
                    for p in proteins:
                        writer.writerow([p])
                        peptides = df_output[['found_match','type','PTM','title']][df_output['protein']==p].values
                        done = []
                        for pep,types,ptm,t in peptides:
                            if pep+ptm in done:
                                continue
                            done.append(pep+ptm)
                            unique = ''.join(list(df_plot['uniqueness'][(df_plot['found_match']==pep) & (df_plot['taxon']==taxon)].values))
                            
                            if 'U_' in unique:
                                pep = '*'+pep
                            if types != 'Original':
                                pep = '~'+pep
                            writer.writerow([pep,ptm,t])
        print('Done saving')
        file_loc = 'ZooMSMS_results_'+file_name+'.csv'
    except:
        with open(path/'Output_Classicol'/sample_path/'ZooMSMS_results.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['ZooMSMS analysis output (~: isoBLAST match, *: Unique PSM)'])
            writer.writerow(['Isoblast analysis revealed '+ str(len(set(df_output['found_match'].values)))+' unique peptide matches with the database'])
            for taxon,size in ranking_taxon:
                writer.writerow(' ')
                writer.writerow(['Taxonomic match: '+taxon+' Average Score= '+str(size)])
                animals = sorted(list(set(df_output['animal'][df_output['taxon']==taxon].values)))
                writer.writerow(['Containing these species:'+', '.join(animals)])
                for a in animals:
                    proteins = list(set(df_output['protein'][df_output['animal']==a].values))
                    for p in proteins:
                        writer.writerow([p])
                        peptides = df_output[['found_match','type','PTM','title']][df_output['protein']==p].values
                        done = []
                        for pep,types,ptm,t in peptides:
                            if pep+ptm in done:
                                continue
                            done.append(pep+ptm)
                            unique = ''.join(list(df_plot['uniqueness'][(df_plot['found_match']==pep) & (df_plot['taxon']==taxon)].values))
                            
                            if 'U_' in unique:
                                pep = '*'+pep
                            if types != 'Original':
                                pep = '~'+pep
                            writer.writerow([pep,ptm,t])
        print('Done saving')
        file_loc = 'ZooMSMS_results.csv'
    return path/'Output_Classicol'/sample_path/file_loc

def make_connection_graph(df_output, file_name,final_output,found_animals,taxonomy):
    taxon = []
    all_taxonomy = {k:v for k,v in taxonomy.items() if k in found_animals}
    for i in df_output['animal'].values:
        taxon.append(final_output[i])
    df_output['taxon']=taxon

    df_plot = pd.DataFrame(columns = ['taxon','found_match','mascot_peptide','protein'])
    temp = [tuple(el) for el in df_output[['taxon','found_match','mascot_peptide','protein']].values]
    for i in set(temp):
        df_temp = pd.DataFrame(np.array(i).reshape(1,-1),columns = ['taxon','found_match','mascot_peptide','protein'])
        df_plot = pd.concat([df_plot,df_temp],ignore_index=True)

    pro = []
    temps = {}
    for i in df_plot['protein'].values:
        if i in temps:
            pro.append(temps[i])
        else:
            temp = i.split('[')[0]
            if '=' in temp:
                temp = i.split('=')[1]
                pro.append('_'.join(temp.split(' ')[0]))
                temps[i]='_'.join(temp.split(' ')[0])
            else:
                pro.append('_'.join(temp.split(' ')[1:]))
                temps[i]='_'.join(temp.split(' ')[1:])
    df_plot['protein']=pro
    
    unique = []
    temps = {}
    for i in df_plot['mascot_peptide'].values:
        if i in temps:
            unique.append(temps[i])
        else:
            temp = sorted(list(set(df_plot['taxon'][df_plot['mascot_peptide']==i])))
            if len(temp)==1:
                unique.append('U_'+' + '.join(temp))
                temps[i]='U_'+' + '.join(temp)
            else:
                if temp[0] in all_taxonomy:
                    start = all_taxonomy[temp[0]]
                    for k,v in all_taxonomy.items():
                        start = set(v)&set(start)
                    for t in all_taxonomy[temp[0]]:
                        if t in start:
                            unique.append('Shared_'+t[1])
                            temps[i]='Shared_'+t[1]
                            break
                else:
                    unique.append('other')
                    temps[i]='other'
    df_plot['uniqueness']=unique
    return df_plot, df_output
 
def rescore(path, file_loc,f_name,sample_path,taxonomy):
    file_name = file_loc.stem+'.csv'
    file_name_print = file_name.replace('.csv','')
    try:
        df = pd.read_csv(path/'Output_Classicol'/sample_path/file_name,names=['a','b','c'])
    except:
        df = pd.read_csv(path/'Output_Classicol'/sample_path/'ZooMSMS_result.csv',names=['a','b','c'])
    df = df.iloc[2:]
    
    add = []
    brackets = ''
    s=0
    score = []
    group = 0
    groups=[]
    pro = ''
    prot = []
    for i in df['a'].values:
        if '[' in i:
            keep = i
            i = i.split('[')[-1]
            i = i.split(']')[0]
            try:
                i = float(i)
                s = i
                brackets = ''
                group += 1
            except:
                pro=keep
                brackets = i
            score.append(s)
            prot.append('')
            add.append('')
            groups.append(group)
        elif 'OS=' in i:
            keep = i
            i = i.split('OS=')[-1]
            i = i.split(' OX')[0]
            score.append(s)
            pro=keep
            brackets = i
            prot.append('')
            add.append('')
            groups.append(group)
        elif 'Containing' in i:
            prot.append('')
            add.append('')
            score.append(0)
            groups.append(group)
        else:
            prot.append(pro)
            add.append(brackets)
            score.append(s)
            groups.append(group)
    amb_pep = [num if '~' in num else '' for num in df['a'].values]
    peps = [num.replace('~','') for num in df['a'].values]
    u_pep = [True if '*' in num else False for num in df['a'].values]
    peps = [num.replace('*','') for num in peps]
    df['a']=peps
    df['species']=add
    df['score']=score
    df['group']=groups
    df['protein']=prot
    df['ambig_pep']=amb_pep
    df['unique_peps']=u_pep
    df = df[df['species']!='']
    df.columns = ['peptides', 'PTM','title','species','score','group','protein','ambig_peps','unique_peps']

    
    columns = list(set(df['peptides'].values))
    z = []
    y = sorted(list(set(df['species'].values)))
    group_to_sp = {k:v for k,v in df[['species','group']].values}
    sp_to_group = {}
    for k,v in group_to_sp.items():
        if v in sp_to_group:
            sp_to_group[v]=sp_to_group[v]+ '+' +k
        else:
            sp_to_group[v]=k
            
    for a in y:
        temp = [1 if i in set(df['peptides'][df['species']==a].values) else 0 for i in columns]
        z.append(temp)
    top5 = {}
    top_sp = list(set(df['species'].values))
    for i in top_sp:
        sc = list(set(df['score'][df['species']==i].values))[0]
        top5[i]=sc
    top5sc = sorted(list(top5.values()))[::-1]
    top_x = 10
    top5sc = top5sc[0:top_x]
    top5 = {k:v for k,v in top5.items() if v in top5sc}
    score_to_cluster = {0:1}
    sp_to_cl = {el:0 for el in top5.keys()}
    y = top5.keys()
    for k,v in score_to_cluster.items():
        species = [sp for sp in y if sp_to_cl[sp]==k]
        if len(species)>1:
            temp = df[df['species'].isin(species)]
            common = []
            for i in set(temp['peptides'].values):
                if len(set(temp['species'][temp['peptides']==i].values))==len(species):
                    common.append(i)
            if len(common)==len(set(temp['peptides'])):
                common = []
            temp = temp[temp['peptides'].isin(common)==False]
            count = list(temp['peptides'].values)
            count = [count.count(num)/len(set(temp['species'][temp['peptides']==num])) for num in temp['peptides'].values]
            combine_species = []
            for i in temp['group'].values:
                sp = '+'.join(list(set(temp['species'][temp['group']==i])))
                combine_species.append(sp)
            temp['species']=combine_species
            temp['count']=count
            temp = temp.drop_duplicates()
            fig = px.bar(temp, x="peptides", y='count', color="species", title=file_name_print)
            fig.update_layout(height=600, width=1500, barmode = 'stack', xaxis={'categoryorder':'total descending'})
            fig.update_xaxes(showticklabels=False)
            name_file = 'Barplot_uniquePeptides_'+file_name_print+'.html'
            try:
                fig.write_html(path/'Output_Classicol'/sample_path/name_file)
            except:
                fig.write_html(path/'Output_Classicol'/sample_path/'Barplot_uniquePeptides.html')
            # plotly.offline.plot(fig)
            time.sleep(5)
            sp = list(set(list(temp['species'].values)))
            spy2 = []
            all_peptides = list(set(list(temp['peptides'].values)))
            array1 = [1 if pep in all_peptides else 0 for pep in all_peptides]
            seq_concat = '_'.join(all_peptides)
            uniques = [num for num in all_peptides if len(temp['species'][temp['peptides']==num].values)==1]
            weights_bc = [1/seq_concat.count(val) if (val not in uniques) or (seq_concat.count(val)>1) else 2 for val in all_peptides]
            for t in sp:
                temp_pep = temp['peptides'][temp['species']==t].values
                array2 = [1 if pep in temp_pep else 0 for pep in all_peptides]
                score = braycurtis(array1, array2,w=weights_bc)
                spy2.append(1-score)
            spy1 = [temp['score'][temp['species']==t].values[0] for t in sp]
            df_scores = pd.DataFrame()
            df_scores['species']=sp
            df_scores['score original']=spy1
            df_scores['Rescore']=spy2
            df_scores = df_scores.sort_values(by='score original', ascending=False)
            fig = go.Figure()
            fig.add_trace(
                go.Scatter(x=df_scores['species'], y=df_scores['score original'],name='Original score'))
            fig.add_trace(
                go.Scatter(x=df_scores['species'], y=df_scores['Rescore'],name='Rescore'))
            fig.update_layout(title = 'Original score VS Rescore of '+file_name_print,height=600, width=900)
            try:
                name_file = 'Rescore_'+file_name_print+'.html'
                fig.write_html(path/'Output_Classicol'/sample_path/name_file)
            except:
                fig.write_html(path/'Output_Classicol'/sample_path/'Rescore.html')
            # plotly.offline.plot(fig)
            time.sleep(3)
        else:
            df_scores = pd.DataFrame()
            df_scores['species']=['']
            df_scores['score original']=['']
            df_scores['Rescore']=['']
            fig = px.bar(df, x=['all peptides'], y=[1],color=species, title=file_name_print)
            try:
                name_file = 'Barplot_uniquePeptides_'+file_name_print+'.html'
                fig.write_html(path/'Output_Classicol'/sample_path/name_file)
            except:
                fig.write_html(path/'Output_Classicol'/sample_path/'Barplot_uniquePeptides.html')
            # plotly.offline.plot(fig)
            time.sleep(3)
    to_summ = make_output_file_after_rescoring(path,df,df_scores, f_name,sample_path,taxonomy)
    return to_summ

def find_child_nodes(taxonomy,taxon):
    #find all species in taxonomy under same LCA
    children = []
    for k,v in taxonomy.items():
        if taxon==v[1][1]:
            if k!=taxon and k.isdigit()==False and 'unclassified' not in k and 'sp.' not in k:
                children.append(k)
    
    return children

def go_to_species(c,species_list,taxonomy):
    for x in c:
        if x.lower()=='environmental samples' or 'sp.' in x:
            continue
        t_lin = taxonomy[x]
        if 'species' in t_lin[0][0]:
            species_list.append(x)
            try:
                gtsp = find_child_nodes(taxonomy,x)
                species_list = species_list +gtsp
            except:
                continue
        else:
            try:
                gtsp = find_child_nodes(taxonomy,x)
                gts = []
                for g in gtsp:
                    g_no_numbers = ''.join([num for num in g if num.isdigit()==False])
                    if g_no_numbers==g:
                        gts.append(g)
                species_list = species_list + go_to_species(gts,[],taxonomy)
            except:
                continue
    return species_list

def make_output_file_after_rescoring(path,df_og,df_rescore, file_name,sample_path,taxons):
    print('generating final output file')
    name_file = 'Taxonomic_results_after_rescoring_'+file_name+'.csv'
    try:
        with open(path/'Output_Classicol'/sample_path/name_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['Species','Species_Rank','Original_score','Rescored','Unique_peptide_amongst_candidates',
                             'isoBLAST_generated_Peptide','Protein','Peptide','PTMs','Representative_PSM_title'])
            for p,ptm,sp,t,sc,gr,amb,pro,u in df_og[['peptides','PTM','species','title','score','group','ambig_peps','protein','unique_peps']].values:
                if len(amb)>0:
                    amb = True
                else:
                    amb = False
                rescore = 'Not amongst top candidates used for rescoring'
                for spr,resc in df_rescore[['species','Rescore']].values:
                    if sp in spr:
                        rescore = resc
                        break
                writer.writerow([sp,gr,sc,rescore,u,amb,pro,p,ptm,t])
    except:
        with open(path/'Output_Classicol'/sample_path/'Taxonomic_results_after_rescoring.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['Species','Species_Rank','Original_score','Rescored','Unique_peptide_amongst_candidates',
                             'isoBLAST_generated_Peptide','Protein','Peptide','PTMs','Representative_PSM_title'])
            for p,ptm,sp,t,sc,gr,amb,pro,u in df_og[['peptides','PTM','species','title','score','group','ambig_peps','protein','unique_peps']].values:
                if len(amb)>0:
                    amb = True
                else:
                    amb = False
                rescore = 'Not amongst top candidates used for rescoring'
                for spr,resc in df_rescore[['species','Rescore']].values:
                    if sp in spr:
                        rescore = resc
                        break
                writer.writerow([sp,gr,sc,rescore,u,amb,pro,p,ptm,t])
    summary=[]
    for gr in range(1,4):
        temp = df_og[df_og['group']==gr]
        if len(temp)==0:
            summary = summary+['']*6
            continue
        temp_sp=sorted(list(set(temp['species'].values)))
        taxa=taxons[temp_sp[0]]
        for x in temp_sp:
            t = set(taxons[x])&set(taxa)
            new_t = []
            for i in t:
                if i in taxa:
                    new_t.append(i)
            taxa=new_t
        for t in taxons[temp_sp[0]]:
            if t in taxa:
                taxa = [t]
                break
        taxa = taxa[0][1]+' ('+', '.join(temp_sp)+')'
        rescore='Not rescored'
        for spr,resc in df_rescore[['species','Rescore']].values:
            if spr in temp_sp:
                rescore = resc
                break
        summary = summary + [taxa,
                             list(set(temp['score'][temp['species'].isin(temp_sp)].values))[0],
                             rescore,
                             len(set(temp['peptides'].values)),
                             len(set(temp['ambig_peps'])),
                             len(temp[temp['unique_peps']==True]),
                             ]
    return summary

def background_check(c,ip,taxonomy):
    out = {}
    for q in c:
        m = go_to_species([q],[],taxonomy)
        if len(set(m)&set(ip))==0:
            out[q] = m
    return out

def make_sunburst_with_missing2(dfs,ip_animals,found_animals,path,sample_path,taxonomy):
    print('Start making sunburst with missingness')
    df_missing = dfs[['mascot_peptide','found_match','animal','type']]
    #ip_animals are all animals that were searched against
    # additional = {}
    additional_animals = []
    
    overlap = set()
    all_taxonomy = {k:v for k,v in taxonomy.items() if k in found_animals}
    for k,v in all_taxonomy.items():
        if k in ip_animals:
            overlap=overlap|set([el[1] for el in v])
    overlap_all_input_animals=set()
    for k,v in taxonomy.items():
        if k in ip_animals:
            overlap_all_input_animals=overlap_all_input_animals|set([el[1] for el in v])
    overlap=list(overlap)
    overlap_all_input_animals=list(overlap_all_input_animals)
    print(overlap)
    print('{} taxa already assigned'.format(len(overlap)))
    taxonomy_temp = {k:v for k,v in taxonomy.items() if len(v)>2}
    for f_child in overlap:
        if taxonomy[f_child][0][0]=='species':
            children = find_child_nodes(taxonomy,f_child)
            children = [el for el in children if 'sp.' not in el and el.isdigit()==False]
            missing = background_check(children,ip_animals,taxonomy)
            for k,aa in missing.items():
                if 'unclassified' not in k:#len(aa)==0 and 
                    additional_animals.append(k)
        else:
            #take all taxons not in database but linked to found lineages. Check if higher up lineage also in results, because sometimes same name for different taxa
            additional_animals = list(set(additional_animals + [k for k,v in taxonomy_temp.items() if v[1][1]==f_child and v[2][1] in overlap and k not in overlap and k not in overlap_all_input_animals]))
    print('Looking for species related to outcome yielded {} missing taxa'.format(len(additional_animals)))
    additional_animals = [el for el in additional_animals if any(element[0]=='order' for element in taxonomy[el])]
    #take now highest taxon representing missing lineages from analysis
    #subspecies next taxon is in overlap, other missing species highest taxon needs to be in overlap

    all_missing_animals = found_animals+additional_animals
    for m in set(additional_animals):
        df_missing_add = np.array(['','',m,'Original'])
        df_missing_add = pd.DataFrame(df_missing_add.reshape(1,-1),columns=df_missing.columns)
        df_missing = pd.concat([df_missing,df_missing_add],ignore_index=True)

    print('done looking')
    labels = []
    values = []
    parents = []
    done = []
    values_iso = []
    braycurtis_dist = []
    
    all_peptides = list(set(list(df_missing['found_match'].values)))#mascot_peptide
    array1 = [1 if pep in all_peptides else 0 for pep in all_peptides]
    seq_concat = '_'.join(all_peptides)
    weights_bc = [1/seq_concat.count(val) for val in all_peptides]         
    print('start assigning branch scores missing')
    
    combo_taxon = []
    for key,val in taxonomy.items():
        if key in all_missing_animals:
            combo_taxon = combo_taxon+val
    combo_taxon=set(combo_taxon)
    
    for animal in all_missing_animals:#first the found animals are mapped, than the missing species are mapped to the existing tree
        
        a_taxon = taxonomy[animal]
        labels.append(animal)
        if animal in additional_animals:
            values.append('tbd')
            values_iso.append(0)
            braycurtis_dist.append('nan')
            for t in a_taxon:
                if t[1]==animal:
                    continue
                if t[1] in parents or t[1] in labels: #for subspecies it will match the label, for missing species it will match the parent nodes
                    parents.append(t[1])
                    break
            continue #we don't need to go in deeper, the parent node should be present already
        else:
            values.append(len(set(df_missing['found_match'][df_missing['animal']==animal].values)))
            values_iso.append(len(set(df_missing['found_match'][(df_missing['animal']==animal) & (df_missing['type'].str.contains('isobaric')==True)].values)))
            temp_peptides = set(df_missing['found_match'][df_missing['animal']==animal].values)#mascot_peptide
            array2 = [1 if pep in temp_peptides else 0 for pep in all_peptides]
            score = braycurtis(array1, array2,w=weights_bc)
            braycurtis_dist.append(1-score)#columns need to be the same
        
        ai = []
        stop = False
        
        previous_taxon = ''
        for t in a_taxon:
            if stop == True:
                break
            if previous_taxon == t[1]:
                continue
            else:
                previous_taxon = t[1]
            if (t in combo_taxon or t[1] in parents) and t[1]!=animal:
                animals_involved = [key for key,val in all_taxonomy.items() if t in val and key in df_missing['animal'].values]
                if sorted(animals_involved) == sorted(ai):
                    if len(animals_involved) == len(set(df_missing['animal'])):
                        parents.append('')
                        stop = True
                        continue
                ai = animals_involved
                add = len(set(df_missing['found_match'][df_missing['animal'].isin(animals_involved)].values))-1
                add_iso = len(set(df_missing['found_match'][(df_missing['animal'].isin(animals_involved)) & (df_missing['type'].str.contains('isobaric')==True)].values))
                if t[1] in all_missing_animals:
                    if t[1] in animal and t[0]=='species' and t[1]!=animal:
                        parents.append(t[1])
                        break
                if [t[1]] in done and 'species' not in t:
                    parents.append(t[1])
                    break
                done.append([t[1]])
                parents.append(t[1])
                labels.append(t[1])
                values.append(add) 
                values_iso.append(add_iso)
                temp_peptides = set(df_missing['found_match'][df_missing['animal'].isin(animals_involved)].values)#mascot_peptide
                #we want only the shared peptides/level. At species level the uniqueness does count
                unique_peptides = []
                animals_involved_og = [key for key,val in all_taxonomy.items() if t in val and key in dfs['animal'].values]
                for peptide_test in temp_peptides:
                    test_pep = set(dfs['animal'][(dfs['found_match']==peptide_test)&(dfs['animal'].isin(animals_involved_og))].values)#mascot_peptide
                    if len(test_pep)!= len(animals_involved_og):
                        unique_peptides.append(peptide_test)
                temp_peptides = temp_peptides^set(unique_peptides)
                array2 = [1 if pep in temp_peptides else 0 for pep in all_peptides]
                score = braycurtis(array1, array2,w=weights_bc)
                braycurtis_dist.append(1-score)#columns need to be the same
        if len(parents)<len(values):
            parents.append('')
    temporary_dataframe = pd.DataFrame()
    temporary_dataframe['values']=values
    temporary_dataframe['parents']=parents
    temporary_dataframe['labels']=labels
    temporary_dataframe['values_iso']=values_iso
    temporary_dataframe['braycurtis_dist']=braycurtis_dist
    vals = []
    print('adjusting for missing species')
    for k,v in temporary_dataframe[['parents','values']].values:
        if v == 'tbd':
            done = False
            for k2,v2 in temporary_dataframe[['labels','values']][temporary_dataframe['labels']==k].values:
                if k==k2 and v2!= 'tbd':
                    vals.append(v2)
                    done = True
                    break
            if done == False:
                print('Unable to assign score to {}'.format(k))
                vals.append(0)
        else:
            vals.append(v)
    temporary_dataframe['values']=vals
    values = []
    parents = []
    labels = []
    values_iso = []
    braycurtis_dist = []
    done = []
    for x in temporary_dataframe.values:
        if (x[1],x[2]) not in done and x[2] not in labels:
            values.append(x[0])
            parents.append(x[1])
            labels.append(x[2])
            values_iso.append(x[3])
            braycurtis_dist.append(x[4])
            done.append((x[1],x[2]))
    
    print('Start making figure')
    fig = go.Figure()
    fig.add_trace(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="remainder",
        meta = values_iso,
        marker=dict(
        colors=braycurtis_dist,
        coloraxis='coloraxis1',
        colorscale='Jet',
        showscale=True,
        cmid=0.5),
    
    hovertemplate='<b>%{label} </b> <br> Taxon score: %{color:.3f}<br> Total peptide count: %{value:.0f} <br> isoBLAST peptidoforms: %{meta:.0f}'))
    fig.update_layout(title = 'Score of output species to sample with species not in DB: '+file_name)
    fig.update_layout(autosize=True,margin = dict(t=30, l=0, r=0, b=0))
    fig.update_layout(coloraxis_colorbar_title='Score')
    try:
        name_file = 'sunburst_including_missing_species_'+file_name+'.html'
        fig.write_html(path/'Output_Classicol'/sample_path/name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'sunburst_including_missing_species.html')
    # plotly.offline.plot(fig)
    # time.sleep(2)
    return temporary_dataframe
    
def load_manual_files(path,test_file,AA_codes):
    #peprec file columns = ['peptide','ptm'] ptm= 0|N-term|1|Oxidation|-1|C-term
    print('open file {}'.format(test_file))
    df = pd.read_csv(test_file, header=0)
    df = df.fillna('')
    df.columns = ['pep_seq','pep_var_mod']
    file_name = str(test_file).split('/')[-1]
    file_name = file_name.split('.')[0]
    df['pep_scan_title']=[file_name+'_'+str(num) for num in range(0,len(df))]#add_unique_title
    df['prot_desc']=['unknown_protein']*len(df)
    #turn modification in the style of mascot
    #'pep_var_mod_pos' 0.001100.0
    #'pep_var_mod' Oxidation (P); Deamidated(NQ)
    new_mod = []
    new_mod_pos = []
    for p,pvm in df[['pep_seq','pep_var_mod']].values:
        if pvm == '':
            new_mod.append('')
            new_mod_pos.append('0.'+'0'*len(p)+'.0')
        else:
            pvm = pvm.split('|')
            new_pvm = {}
            for i in range(0,len(pvm),2):
                new_pvm[int(pvm[i])]=pvm[i+1]
            temp_mod={}
            adj_mod_pos=['0']*len(p)
            n_term='0.'
            c_term='.0'
            for k,v in new_pvm.items():
                if k==0:
                    n_term='1.'
                    v=v+' (N-term)'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
                elif k==-1:
                    c_term='.1'
                    v=v+' (C-term)'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
                else:
                    adj_mod_pos[k-1]='1'
                    v=v+' ('+p[k-1]+')'
                    if v in temp_mod.keys():
                        temp_mod[v]=temp_mod[v]+1
                    else:
                        temp_mod[v]=1
            adj_mod = '; '.join([str(v)+' '+k if v>1 else k for k,v in temp_mod.items()])
            new_mod.append(adj_mod)
            new_mod_pos.append(n_term+''.join(adj_mod_pos)+c_term)
    df['pep_var_mod']=new_mod
    df['pep_var_mod_pos']=new_mod_pos
    
    df_4_uni= df
    df2=df
    unimod_db, unimod = do_unimod(path,df_4_uni['pep_var_mod'].values)
    ids = {}
    for p, a,m in unimod_db[['PTM','AA','mass']].values:
        if a=='N-term':
            a='!'
        elif a=='C-term':
            a='*'
        add = '?'
        while add+a in AA_codes.keys():
            add+='?'
        if a=='!' or a=='*':
            AA_codes[add+a]=float(m)
        else:
            AA_codes[add+a]=AA_codes[a]+float(m)
        ids[add+a]=p
    
    return df2, unimod_db, ids,AA_codes


def do_alignment_threaded(all_seqs, seq_to_name, calculator, Z_distance_csv, sequence_db, cpu_count):
    with ProcessPoolExecutor(max_workers=cpu_count) as executor:
        done = []
        X = []
        all_seqs_name = [seq_to_name[t] for t in all_seqs if seq_to_name[t] in Z_distance_csv.columns]
        Z_temp_all = Z_distance_csv.loc[all_seqs_name]
        Z_temp_all = Z_temp_all[all_seqs_name]#inlcude only seqs that are relevant
        saving_later = False
        with tqdm(all_seqs, desc="Starting on alignment...") as pbar:
            for seq in pbar:
                pbar.set_description(f"Performing alignment on {seq_to_name[seq]}")
                temp = []
                align_seqs = []
                name_seq = seq_to_name[seq]
                if name_seq in all_seqs_name:
                    Z_temp = Z_temp_all.loc[[name_seq]]#make small subset for quick lookup
                    for t in all_seqs:
                        name_t = seq_to_name[t]
                        if name_t not in done and (name_seq not in Z_temp.columns
                                                or name_t not in Z_temp.columns):#if not considered this run and not in previous runs
                            align_seqs.append(t)
                        elif name_seq in Z_temp.columns and name_t in Z_temp.columns:#if considered previously
                            if list(Z_temp[name_t])[0]=='nothing':#if it was not matched previously
                                if name_t not in done:#if it has not been calculated this run already
                                    align_seqs.append(t)
                else:
                    align_seqs = [el for el in all_seqs if seq_to_name[el] not in done]
                if len(align_seqs)>0:
                    saving_later = True
                    results = []
                    futures = [executor.submit(thread_align, seq, t, calculator, sequence_db[t]) \
                            for t in set(align_seqs)]
                    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing", leave=False):
                        results.append(future.result())
                    results = {val: key for key, val in results}
                
                for t in all_seqs:
                    name_t = seq_to_name[t]
                    if name_t in done:#find the distance in the distance matrix that is being calculated
                        temp.append(float(X[all_seqs.index(t)][all_seqs.index(seq)]))
                    elif name_seq == name_t:#if same sequence, assign 0
                        temp.append(0)
                    elif t not in align_seqs and seq not in align_seqs:#if already calculated before
                        matching = list(Z_temp[name_t].values)[0]
                        if matching!='nothing':
                            temp.append(float(matching))
                        else:
                            temp.append(float(results[sequence_db[t]]))
                    else:
                        temp.append(float(results[sequence_db[t]]))
                done.append(name_seq)
                X.append(temp)
    
    return np.array(X),saving_later

def ClassiCOL_analysis(
        path: pathlib.Path,
        sample_path: str,
        file_name: str,
        sequence_db: dict[Seq, str],
        df2: pd.DataFrame,
        unimod_db, 
        ids,
        data_to_remember,
        AA_codes: dict[str, float],
        lim_tax: str, 
        demo_tf: bool,
        cpu_count: int,
    ):
    summary_output_file: list[typing.Any] = []
    ######Begin of contamination removal###########
    print('getting rid of Trypsin, Lys-C')
    # contamination = []
    # for i in df2['prot_desc'].values:
    #     contamination.append('keratin' in str(i).lower())
    # df2['contamination']=contamination
    # df2 = df2[df2['contamination']==False]
    
    df2 = df2.drop(['prot_desc'], axis=1)
    drop = []
    done = []
    print('Contaminants have been deleted')
    print('Deleting double peptides for faster classification')
    for p,v,vp in df2[['pep_seq','pep_var_mod','pep_var_mod_pos']].values:
        if (p,v,vp) not in done:
            done.append((p,v,vp))
            drop.append(False)
        else:
            drop.append(True)
    df2['double']=drop
    df2=df2[df2['double']==False]
    df2 = df2.drop_duplicates()
    sampleprep = ''.join([str(num) for num,el in sequence_db.items() if 'trypsin' in el.lower() or 'pseudomonas' in el.lower()])
    insampleprep = []
    for i in df2['pep_seq'].values:
        if i in sampleprep:
            insampleprep.append(True)
        else:
            insampleprep.append(False)
    df2['insp']=insampleprep
    df2 = df2[df2['insp']==False] #do not need to check these peptides
    ##################End of contamination removal########################
    ##################Start of batch search remembering###################
    print('Adjusting data for batch search')
    title_to_ptm_og = {}
    if len(data_to_remember) == 0:
        already = [False]*len(df2)
    else:
        already = []
        for k,v,v2,t in df2[['pep_seq','pep_var_mod_pos','pep_var_mod','pep_scan_title']].values:
            v = ''.join([el if (el == '0' or el=='.') else '1' for el in v])
            if any(mp==k and pl==v and ptm==v2 for mp,pl,ptm in data_to_remember[['mascot_peptide','PTM_loc','PTM_og']].values):
                already.append(True)
                title_to_ptm_og[t]=v2
            else:
                already.append(False)
    df2['already']=already
    df_already_found = df2[df2['already']==True]
    print('Reduced the dataframe by {} peptides because they were already matched'.format(len(df_already_found)))
    df2 = df2[df2['already']==False]
    # df2 = df2[df2['pep_scan_title'].isin(df_already_found['pep_scan_title'].values)==False]
    ######################End of remembering peptides from batchsearch##############################
    print('making matrix')
    mass_matrix = make_matrix(AA_codes)
    print('Binary matrix creation successful')
    
    print('starting isoBLAST search')
    columns=['animal','protein','mascot_peptide','found_match','switch','location','type','PTM','title']
    input_animals,skip_animals,taxonomy = animals_from_db_input(sequence_db, lim_tax,demo_tf,path)
    input_animals = sorted(input_animals)
    print(f'Using {cpu_count} CPUs')
    total_peptide_begin = len(df2)
    check_all_psms = False
    
    unimod_masses = [float(0)]+[-t for t in unimod_db['mass'].values] #ptms added is lower backbone mass
    unimod_masses = unimod_masses+[num - t for num in unimod_masses for t in unimod_db['mass'].values] #2 ptms added is lower backbone mass 
    unimod_masses = np.array(sorted(list(set(unimod_masses))), dtype=np.float64)

    unimod_lookup = {
        (row['PTM'], row['AA']): row['mass']
        for _, row in unimod_db.iterrows()
    }
 
    seq_db = {str(key) : val for key, val in sequence_db.items()}
    df_out = find_isoblast(
        input_animals,
        cpu_count,
        df2,
        skip_animals,
        seq_db,
        unimod_db,
        mass_matrix,
        ids,
        check_all_psms,
        AA_codes,
        columns,
        unimod_masses,
        unimod_lookup)
    ####################ISOblast is now finished#############################
    print('Adding position information')
    temp = {(a1,a2):''.join([el if (el == '0' or el=='.') else '1' for el in a3]) for a1,a2,a3 in df2[['pep_seq','pep_scan_title','pep_var_mod_pos']].values}
    add = [temp[(a1,a2)] for a1,a2 in df_out[['mascot_peptide','title']].values]
    df_out['PTM_loc']=add
    ##############Positional inforamtion has been added#########################
    print('Removing spectra with isoBLAST linked to trypsin and Lys-C')
    linked_to_contaminants = []
    for p,t in df_out[['protein','title']].values:
        if 'Pseudomonas' in p or 'trypsin' in p.lower():
            linked_to_contaminants.append(t)
    dfs = df_out[df_out['title'].isin(linked_to_contaminants) == False]
    #################Isobaric peptides to contaminants removed#############################
    title_already_found = df_already_found['pep_scan_title'].values
    retreived = pd.DataFrame(columns=dfs.columns)
    for t in set(title_already_found):
        temp = df_already_found[df_already_found['pep_scan_title']==t]
        done = []
        for k,v,v2 in temp[['pep_seq','pep_var_mod_pos','pep_var_mod']].values:
            if (k,v,v2) in done:
                continue
            done.append((k,v,v2))
            v = ''.join([el if (el == '0' or el=='.') else '1' for el in v])
            temp_found = data_to_remember[(data_to_remember['mascot_peptide']==k) &
                                          (data_to_remember['PTM_loc']==v) &
                                          (data_to_remember['PTM_og']==v2)]
            temp_found['title']=[t]*len(temp_found)
            retreived = pd.concat([retreived,temp_found],ignore_index=True)
    print('Adding {} peptides from batch search results'.format(len(retreived)))
    dfs = pd.concat([dfs,retreived],ignore_index=True)
    dfs['PTM_og']=[list(df2['pep_var_mod'][df2['pep_scan_title']==el].values)[0] if el in df2['pep_scan_title'].values else title_to_ptm_og[el] for el in dfs['title'].values]
    ##################Added the peptides from the batchsearch####################################
    print('Reducing the data to 1 peptide sequence/peptide')
    inds = dfs[['protein','found_match','type']].drop_duplicates().index
    dfs = dfs[dfs.index.isin(inds)]
    print('Removing species that have less than 10% of total peptides')#would be done otherwise after the classification, so will go quicker now
    remove_species = []
    all_peps = set(dfs['found_match'].values)
    for el in set(dfs['animal']):
        if len(set(dfs['found_match'][dfs['animal']==el].values))<len(all_peps)*0.1:
            remove_species.append(el)
    dfs = dfs[dfs['animal'].isin(remove_species)==False]
    
    print('Removing single hit wonders')
    keep = []
    for pro in set(dfs['protein'].values):
        temp = set(dfs['found_match'][dfs['protein']==pro].values)
        if len(temp)<=10:
            drop_peps = [el for el in temp if any(el in t and el!=t for t in temp)]
            temp = temp^set(drop_peps) #if low amount of peptides remove the ones that have a ladder effect. Else there could be accidental single hit wonders
        if len(temp)>1: #>2 peptides per proteins
            keep.append(pro)
    dfs = dfs[dfs['protein'].isin(keep)]
    index, locs_dict, columns, z = calc_distance_collagen(sequence_db, dfs)
    
    #separate the uncertain ones from the analysis
    print('Removing uncertain peptides that are single hit wonders')
    nots = list(dfs['found_match'][dfs['type']!='With uncertainty'].values)
    real_uncertain = [el for el in dfs['found_match'][dfs['type']=='With uncertainty'].values if el not in nots]
    keep_uncertain = dfs[(dfs['type']=='With uncertainty') & (dfs['found_match'].isin(real_uncertain))]
    dfs = dfs[dfs['found_match'].isin(real_uncertain)==False]
    if len(set(title_already_found))>0:
        df_adding = dfs[data_to_remember.columns][dfs['title'].isin(title_already_found)==False]#keep outcome, but only the new ones because we already have the rest stored
    else:
        df_adding = dfs[data_to_remember.columns]
    ################Left uncertain ones out of analysis#######################
    ################Starting on the protein tree##############################
    if len(set(dfs['protein'].values))>1 and len(set(dfs['animal'].values))>1:
        print('Starting on the heatmap and collagen sequence multiple alignment')
        if not (path / 'MISC' / 'collagen_distance.csv').exists():
            print('No distance file found ... \n creating csv file')
            (path / 'MISC').mkdir(parents=True, exist_ok=True)
            with open(path/'MISC/collagen_distance.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                writer.writerow(['index','collagen_seq1','collagen_seq2','distance'])
                writer.writerow(['test',0,0,0])
                writer.writerow(['test2',0,0,0])
                writer.writerow(['test3',0,0,0])
            
        ###############File has been created and is ready to be read########################
        Z_distance_csv = pd.read_csv(path/'MISC/collagen_distance.csv', header=0)      
        Z_distance_csv.index = Z_distance_csv['index'].values
        Z_distance_csv = Z_distance_csv[list(Z_distance_csv.columns)[1:]]
        
        #save each output so this step can go faster in the future
        # matrix = substitution_matrices.load('BLOSUM90')
        print('Calculating distances with {} CPUs'.format(cpu_count))
        names = index
        temp_db = {val:key for key,val in sequence_db.items()}
        all_seqs= [temp_db[key] for key in names]
        seq_to_name = {v:k for k,v in temp_db.items()}
        calculator = DistanceCalculator('blosum90')
        done = []
        distance,saving_out = do_alignment_threaded(all_seqs, seq_to_name, calculator, Z_distance_csv, sequence_db, cpu_count)
        labels = names
        print('end calculating distances')
        
        ############Distances are calculated############################
        all_animals = {el:el for el in input_animals}
        taxon_distance, df_heatmap_values, heat_data,names,animals = make_plots_coverage_per_animal(dfs, sequence_db, all_animals, distance, labels, index, locs_dict, columns, z, names, file_name,sample_path,taxonomy)
        print('find outliers')
        output = find_animals(distance,names,taxon_distance,animals,df_heatmap_values)#eiwit niveau, adapt to filter out bullshit
        #find locations in tree based on protein
        ################Save distance matrix######################
        print('Saving distances ...')
        distance =pd.DataFrame(distance,columns=names,index=names)
        if saving_out ==True:
            print('Start retrieving calculations ...')
            in_dist_csv = list(Z_distance_csv.columns)
            print('Adding new columns')
            new_to_add=[]
            for col in distance.columns:
                if col not in Z_distance_csv.columns:
                    Z_distances: list[float|str] = []
                    new_to_add.append(col)
                    for col2 in in_dist_csv:
                        if col2 in names:
                            Z_distances.append(float(distance[col][col2]))#freshly calculated
                        else:
                            Z_distances.append('nothing')#not calculated this time
                    #adding new calculations to dataframe
                    Z_distance_csv[col]=Z_distances
            if len(new_to_add)!=0:
                print('Adding new rows')
                adapt = [['']*len(Z_distance_csv.columns) for i in range(0,len(new_to_add))]
                Z_distance_csv = pd.concat([Z_distance_csv,pd.DataFrame(np.array(adapt),columns=Z_distance_csv.columns)],ignore_index=True)
                Z_distance_csv.index = Z_distance_csv.columns
                for nums, col in enumerate(new_to_add):
                    if (nums/len(new_to_add))%10==0:
                        print('{} of {} done'.format(nums,len(new_to_add)))
                    Z_distance_csv.loc[col]=Z_distance_csv[col].values
                print('Adding new block')
                for col in new_to_add:
                    for row in new_to_add:
                        Z_distance_csv[col][row]=float(distance[col][row])
            print('Checking if new combinations were made')
            save_cols = Z_distance_csv.columns
            for nums,col in enumerate(save_cols):
                if (nums)%1000==0:
                    print('{} of {} done'.format(nums,len(save_cols)))
                if col not in distance.columns or col in new_to_add:#already added or not calculated this round
                    continue
                Z_distances = Z_distance_csv[col].values
                if 'nothing' in Z_distances:
                    change = []
                    for num,el in enumerate(Z_distances):
                        if el == 'nothing':
                            if col in names and save_cols[num] in names:
                                change.append(float(distance[col][save_cols[num]]))#freshly calculated
                                continue
                        change.append(el)
                    Z_distance_csv[col]=change
                
            ##############Writing the output csv file#####################
            Z_distance_csv = Z_distance_csv.drop(['collagen_seq1','collagen_seq2','distance'],axis=0,errors='ignore')
            Z_distance_csv = Z_distance_csv.drop(['collagen_seq1','collagen_seq2','distance'],axis=1,errors='ignore')
            print('Save as csv file')
            with open(path/'MISC/collagen_distance.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                writer.writerow(['index']+list(Z_distance_csv.columns))
                for q in Z_distance_csv.columns:
                    writer.writerow([q]+list(Z_distance_csv[q].values))
            Z_distance_csv = 'saving memory'
        ######################Done saving matrix##############################
        print('Done saving')
        taxon_distance = pd.DataFrame(np.array(taxon_distance),columns=animals,index=animals)
        print('Starting on the walk down the taxonomic tree')
        
        output_final, recall, dfs = find_animals2_1(distance,names,taxon_distance,animals,df_heatmap_values,output,dfs)#find animal at protein level
        
        print('End of the walk')
        #find uniqueness for each location in tree
        # df_output = pd.DataFrame(columns=dfs.columns)
        
        #rescore the isoblasts to the df_output
        # df_isobaric = df_output[df_output['type'].str.contains('isobaric')]
        # df_other = df_output[df_output['type'].str.contains('isobaric')==False]
        # delete = []
        # for i in df_isobaric['mascot_peptide'].values:
        #     if i in df_other['mascot_peptide'].values:
        #         delete.append(True)
        #     else:
        #         delete.append(False)
        # df_isobaric['del']=delete
        # df_isobaric = df_isobaric[df_isobaric['del']==False]
        # df_isobaric = df_isobaric.drop(['del'],axis=1)
        # df_output = pd.concat([df_other,df_isobaric],ignore_index=True)
        incl_animals = []
        for t in output_final:
            for i in t:
                i = tuple(i)
                if i not in recall:
                    continue
                incl_animals = incl_animals+list(i)
                # animal_pep = dfs[dfs['animal'].isin(i)]
                # df_output = pd.concat([df_output,animal_pep],ignore_index=True)
        df_output = dfs[dfs['animal'].isin(incl_animals)]     
        # delete_animals = []
        # for i in set(df_output['animal'].values):
        #     if len(set(df_output['found_match'][df_output['animal']==i]))==1:
        #         delete_animals.append(i)
        # df_output=df_output[df_output['animal'].isin(delete_animals)==False]
        
        #set the taxonomy
        found_animals = []
        for t in output_final:
            for i in t:
                found_animals = found_animals+list(i)
    
        if len(found_animals)>0:
            final_output = {} 
            for ts in output_final:
                for i in ts:
                    taxon = taxonomy[i[0]]
                    for t in i:
                        taxon = set(taxon)&set(taxonomy[t])
                    for t in taxonomy[i[0]]:
                        if t in taxon:
                            for x in i:
                                final_output[x]=t[1]
                            break
            #If an animal is found that had an uncertain peptide match, than add this peptide to the output. All other uncertain peptides can be discarded
            keep_uncertain = keep_uncertain[keep_uncertain['animal'].isin(set(df_output['animal'].values))]
            keep_uncertain = keep_uncertain[df_output.columns]
            df_output = pd.concat([df_output,keep_uncertain],ignore_index=True)
            print('Determining uniqueness and ambiguity')
            df_plot,df_output = make_connection_graph(df_output, file_name,final_output,found_animals,taxonomy)
            print('Output contains {} species.'.format(len(set(df_output['animal'].values))))
            bc, labels,file_out_final = make_sunburst(dfs,all_animals,df_output,file_name,input_animals,df_plot,path,sample_path,taxonomy)
            ################################################################################
            df_adding['PTM_loc'] = [''.join([el if (el == '0' or el=='.') else '1' for el in element]) for element in df_adding['PTM_loc'].values]
            df_adding['switch']=[', '.join(num) for num in df_adding['switch'].values]
            df_adding['location']=[str(num) for num in df_adding['location'].values]
            df_adding = pd.DataFrame(df_adding, columns=data_to_remember.columns)
            df_adding = df_adding.drop_duplicates()
            df_adding['location']=[eval(num) for num in df_adding['location'].values]
            data_to_remember = pd.concat([data_to_remember,df_adding],ignore_index=True)
            ################################################################################
            to_summary=rescore(path,file_out_final,file_name,sample_path,taxonomy)
            summary_output_file=[lim_tax,total_peptide_begin]+to_summary
        else:
            name_file = 'ZooMSMS_results_'+file_name+'.csv'
            with open(path/'Output_Classicol'/sample_path/name_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                writer.writerow(['ZooMSMS analysis found nothing'])
                summary_output_file=[lim_tax,total_peptide_begin]+['nothing_found']+['']*17
    else:
        name_file = 'ZooMSMS_results_'+file_name+'.csv'
        with open(path/'Output_Classicol'/sample_path/name_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
            writer.writerow(['ZooMSMS analysis found nothing'])
            summary_output_file=[lim_tax,total_peptide_begin]+['nothing_found']+['']*17

    return summary_output_file, data_to_remember, taxonomy

def find_isoblast(
    input_animals,
    cpu_count,
    df2,
    skip_animals,
    seq_db,
    unimod_db,
    mass_matrix,
    ids,
    check_all_psms,
    AA_codes,
    columns,
    unimod_masses,
    unimod_lookup):
    uncertain = retreive_uncertain()
    df_out = pd.DataFrame(columns=columns)
    b = len(df2)
    with ProcessPoolExecutor(max_workers=cpu_count) as process_executor:
        # Map futures to animals for tracking
        future_to_animal = {
            process_executor.submit(
                species_check, a, sequence_db, df2, unimod_db, mass_matrix,
                ids, check_all_psms, AA_codes, uncertain,
                unimod_masses, unimod_lookup, skip_animals, columns
            ): a
            for a in input_animals
        }
        with tqdm(total=len(future_to_animal), desc="Starting isoBLASTing...") as animal_bar:
            for future in as_completed(future_to_animal):
                a = future_to_animal[future]
                animal_bar.set_description(f"Finished isoBLASTing on {a}, {b} peptides were checked")
                try:
                    df_add = future.result()
                    if not df_add.empty:
                        df_out = pd.concat([df_out, df_add], ignore_index=True)
                except Exception as e:
                    animal_bar.write(f"Error processing {a}: {e}")
                finally:
                    animal_bar.update(1)
    return df_out

def species_check(a, sequence_db, df2, unimod_db, mass_matrix, ids, check_all_psms,
                   AA_codes, uncertain, unimod_masses,
                   unimod_lookup, skip_animals, columns):
    df_out = pd.DataFrame(columns=columns)
    
    results = [
        thread_worker(sequence, name, df2, unimod_db, mass_matrix, ids,
                      check_all_psms, AA_codes, uncertain,
                      unimod_masses, unimod_lookup)
        for sequence, name in sequence_db.items()
        if not include_species(name, a, skip_animals)
    ]
    for data in results:
        if False in data:
            continue
        df_add = pd.DataFrame(data=data, columns=columns)
        df_out = pd.concat([df_out, df_add], ignore_index=True)
    return df_out

def isinput_species(name):
    
    if 'OS=' in name:
        anim = name.split('OS=')[-1]
        anim = anim.split(' OX=')[0]
    elif '[' in name:
        anim = name.split('[')[1]
        anim = anim.split(']')[0]
    elif '|' in name:
        anim = name.split('|')[1]
    else:
        anim = name
    return anim

def include_species(name, animal, skip_animals):
    anim = isinput_species(name)
    skip = anim in skip_animals or anim != animal
    return skip

def find_LCA(taxonomy,animals,adjust=False):
    animals = list(animals)
    taxa = {k:taxonomy[k] for k in animals}
    lca = set(taxa[animals[0]])
    for k,v in taxa.items():
        lca = set(v)&lca
    lca0 = taxa[animals[0]]
    for i in lca0:
        if i in lca:
            if adjust == True and i[1] in animals:
                continue
            return i[1]
    return 'root'

def relatedness_finder(considered_species,taxonomy):
    td = {}
    for c in considered_species:
        related = []
        for a in considered_species:
            if a==c:
                continue
            l = find_LCA(taxonomy,[c,a])
            for loc,i in enumerate(taxonomy[c]):
                if l in i:
                    related.append([a,loc])
                    break
        td[c]=[el[0] for el in sorted(related, key = lambda x:x[1])]
    return td

def reiter_grouping(taxonomic_groups,taxonomy):
    tg = {}
    to_check =[]
    related = []
    enough = set()
    for order,species in taxonomic_groups.items():
        if len(species)<4:#go a bit higher to combine, but not too high > still needs to make some sense
            to_check.append(order)
            related = related + species
        else:
            for sp in species:
                enough = enough|set(taxonomy[sp])#restrictions already assigned
    related = relatedness_finder(related, taxonomy)
    change = {}
    if len(to_check)>1:
        for order in to_check:
            for contenders in [[element for element in related[el] if element not in taxonomic_groups[order]][0] for el in taxonomic_groups[order]]:
                flca = find_LCA(taxonomy,taxonomic_groups[order]+[contenders])
                line = [el[0] for el in taxonomy[contenders]]
                if 'class' not in line:
                    inds= 0
                else:
                    inds = line.index('class')
                lines = [el[1] for el in taxonomy[contenders]][inds:]
                if any(flca in el for el in enough)==False and flca not in lines:
                    if flca not in change:
                        change[flca]=taxonomic_groups[order]+[contenders]
                    else:
                        change[flca]=list(set(change[flca]+taxonomic_groups[order]+[contenders]))
                else:
                    break
    species_assigned_new = []
    for k,v in change.items():
        if any(len(set(el)&set(v))>0 for el in change.values() if v!=el):
            matches = [el for el,c in change.items() if len(set(v)&set(c))>0]
            match = []
            for a in matches:
                match = list(set(match + change[a]))
            related = relatedness_finder(match, taxonomy)
            overlapping = set(change[matches[0]])
            for a in matches:
                overlapping = overlapping&set(change[a])
            for m in overlapping:
                closest = related[m]
                for ms,vs in change.items():
                    if closest[0] in vs:
                        tg[ms]=vs
                        species_assigned_new = species_assigned_new + vs
                        break
        else:
            tg[k]=v
            species_assigned_new = species_assigned_new +v
    for k,v in taxonomic_groups.items():
        if any(el in species_assigned_new for el in v)==False:
            species_assigned_new = species_assigned_new +v
            tg[k]=v

    return tg

def parse_input(path_to_output_file,file_extinct,mixture,sequences,taxonomy):
    species_in = animals_from_db_input_mix(sequences)
    try:
        df = pd.read_csv(path_to_output_file/file_extinct, header=0).fillna(0)
    except:
        df = pd.read_csv(path_to_output_file/'Taxonomic_results_after_rescoring.csv', header=0).fillna(0)
    taxonomy_missing = {k:v for k,v in taxonomy.items() if k not in species_in+list(set(df.Species))}
    taxonomy = {k:v for k,v in taxonomy.items() if k in species_in+list(set(df.Species))}
    lca_limit = 'order'#split mammals,fish,birds,etc
    
    if mixture == False:
        print('Single bone')
        a_df = []
        b_df = []
        df = df.sort_values(by=['Original_score'], ascending=False)
        max_score = max(df['Original_score'].values)
        for a,b,c in df[['Species','Species_Rank','Original_score']].values:
            if max_score>0.65 and b<=15:
                if 0.5<c and len(a_df)<6:
                    if b not in b_df and a not in a_df:
                        a_df.append(a)
                        b_df.append(b)
                elif c<0.5:
                    continue
                if b in b_df and a in a_df:#1 candidate with same peptides is good enough
                    continue
                a_df.append(a)
                b_df.append(b)
            else:
                if b in b_df or b>15 or c<max(df['Original_score'].values)-0.1:#1 candidate with same peptides is good enough
                    continue
                a_df.append(a)
                b_df.append(b)
        df=df[df['Species'].isin(a_df)]#no doubles that have the same peptide content
    else:
        print('Potential mixture')
        a_df = list(set(df.Species))
    taxonomic_groups={}
    print('Top species share {} as LCA'.format(lca_limit))
    for sp in set(df.Species):
        taxon = taxonomy[sp]
        if any(lca_limit in el for el in taxon)==False:#maybe adapt this so no random species not considered
            continue
        limit = [el[1] for el in taxon if lca_limit in el][0]
        if limit in taxonomic_groups:
            taxonomic_groups[limit] = taxonomic_groups[limit] + [sp]
        else:
            taxonomic_groups[limit] = [sp]
    taxonomic_groups = reiter_grouping(taxonomic_groups,taxonomy)

    rank_lca = {loc:[find_LCA(taxonomy,a_df[:loc+1]),a_df[:loc+1]] for loc in range(0,len(a_df))}
    
    miss = []#higher taxons need to be included at this stage
    for v in taxonomy.values():
        for i in v:
            if i[1] not in taxonomy:
                miss.append(i[1])
        miss = list(set(miss))
    for i in miss:
        if i in taxonomy_missing:
            taxonomy[i]=taxonomy_missing[i]
            taxonomy_missing.pop(i)
    #iterate within the groups, because otherwise too many 'X' residues which are not completely relevant
    return df, taxonomy, taxonomic_groups,species_in,rank_lca,lca_limit, taxonomy_missing

def find_taxonomic_tree(anim,mammals):
    tree = []
    for i in anim:
        start = mammals[i]
        temp = []
        for q in anim:
            test = mammals[q]
            temp_start = start
            start_taxon = set([el[0] for el in test])&set([ele[0] for ele in start])
            for a in start:
                if a[0] in start_taxon:
                    start_taxon = a[0]
                    s = []
                    do = False
                    for st in start:
                        if start_taxon in st:
                            do=True
                        if do==True:
                            s.append(st)
                    temp_start = s
                    do = False
                    t = []
                    for st in test:
                        if start_taxon in st:
                            do=True
                        if do==True:
                            t.append(st)
                    test=t
            diff = set(temp_start)^set(test)
            temp.append(len(set(start)&diff))
        tree.append(temp)
    branches = {}
    X_test = linkage(tree,'ward')
    i = 0
    one_cluster = False
    while one_cluster == False:
        cluster = list(fcluster(X_test,t=i,criterion='distance'))
        if len(set(cluster))==1:
            one_cluster = True
        if cluster not in branches.values():
            branches[i] = cluster
        i += 0.1
    tree = {}
    for i in sorted(list(branches.keys()))[::-1]:
        cluster_tree = {t:[] for t in branches[i]}
        for l,t in enumerate(branches[i]):
            cluster_tree[t] = sorted(cluster_tree[t]+[anim[l]])
        tree[i]=cluster_tree
    return tree

def aligning(seq1,seq2,matrix):
    
    seq1 = Seq(str(seq1).replace('&',''))
    seq2 = Seq(str(seq2).replace('&',''))
    try:
        aligner = Align.PairwiseAligner(gap_score=-10)
        aligner.substitution_matrix = matrix
        align = aligner.align(seq1,seq2)
        align = next(align)

        return align
    except:
        return ['failed']

def multiple_alignment(consensus,next_seq,calculator,sub_matrix):
    sub_matrix = substitution_matrices.load('BLOSUM90')
    sub_m_X = 5+int(consensus.count('X')/10)
    sub_matrix[22]=[-1]*22+[sub_m_X]+[-1]#X residue
    sub_matrix[:,22]=[-1]*22+[sub_m_X]+[-1]
    sub_matrix[21]=[-10]*21+[sub_m_X*5]+[-10,-10]#Z residue
    sub_matrix[:,21]=[-10]*21+[sub_m_X*5]+[-10,-10]
    keep_c=consensus
    next_seq=str(next_seq)[::-1]
    consensus=str(consensus)[::-1]
    c = []
    temp = ''
    for el in consensus:
        if el=='K' or el=='R':
            if 25>len(temp)>5 and next_seq.count(temp+el)==1 and consensus.count(temp+el)==1 and temp.count('X')<=len(temp)/5:
                c.append(temp+el)
            temp=''
            continue
        temp+=el
    c.append(temp)
    f = 'Z'*5
    ns = 'Z'*5
    for i in c:
        if i in next_seq and i in consensus and next_seq.count(i)==1 and consensus.count(i)==1:
            csplit = consensus.split(i)[0]
            nssplit = next_seq.split(i)[0]
            if len(csplit)>0 or len(nssplit)>0:
                if len(nssplit)==0:
                    nssplit += 'X'*len(csplit)
                if len(csplit)==0:
                    csplit += 'X'*len(nssplit)
            alignment = aligning(str('Z'*3+csplit+i+'Z'*3),str('Z'*3+nssplit+i+'Z'*3),sub_matrix)
            if alignment[0]=='failed':
                print('alignment FAILED')
                return keep_c,'Failed','Failed' #means that the new one is not good
            ns+=alignment[1].replace('Z','')
            f+=alignment[0].replace('Z','')
            next_seq=''.join(next_seq.split(i)[-1])
            consensus=''.join(consensus.split(i)[-1])
            
            while len(f)>len(ns):
                ns+='X'
            while len(f)<len(ns):
                f+='X'
            
                
    if len(next_seq)>0 or len(consensus)>0:
        if len(next_seq)==0:
            next_seq+='X'*len(consensus)

        if len(consensus)==0:
            consensus+='X'*len(next_seq)

        alignment = aligning(str('Z'*3+consensus+'Z'*3),str('Z'*3+next_seq+'Z'*3),sub_matrix)
        if alignment[0]=='failed':
            print('alignment FAILED')
            return keep_c,'Failed','Failed' #means that the new one is not good
        ns+=alignment[1].replace('Z','')
        f+=alignment[0].replace('Z','')
            
        while len(f)>len(ns):
            ns+='X'
        while len(f)<len(ns):
            f+='X'

    next_seq=ns+'Z'*5
    consensus=f+'Z'*5
    
    a1 = consensus.replace('Z','')[::-1]
    a2 = next_seq.replace('Z','')[::-1]
    
    consensus_seq = ''
    old_consensus = ''
    old_seq = ''
    for i in range(0,len(a1)):
        if a1[i]==a2[i]:
            if a1[i]=='-':
                continue
            consensus_seq += a1[i]
            old_consensus += a1[i]
            old_seq += a2[i]
        else:
            consensus_seq += 'X'
            if a1[i]=='-':
                old_consensus += 'X'
            else:
                old_consensus += a1[i]
            if a2[i]=='-':
                old_seq += 'X'
            else:
                old_seq += a2[i]
    return consensus_seq, old_consensus,old_seq

def mapping(seqs_4_align,consensus_df,df):
    #map the peptides and create a theoretical sequence
    done = []
    insilico = [['X']]*len(consensus_df.columns)
    for sp in set(df['Species'].values):
        print('Mapping peptides to {}'.format(sp))
        temp = list(set(df['Peptides'][(df['Species']==sp) & (df['Peptides'].isin(done)==False)].values))
        og_seq = seqs_4_align[sp]
        consensus_aligned_og = ''.join(list(consensus_df.loc[sp].values))
        
        for peptide in temp:
            if peptide not in og_seq:#means that it belongs to a different protein
                continue
            if peptide in consensus_aligned_og:
                done.append(peptide)#not 2 times same 
                positions = [(match.start(),match.end()) for match in re.finditer(peptide, consensus_aligned_og)]
                for i in positions:
                    for x in range(i[0],i[1]):
                        if peptide[x-i[0]] in insilico[x]:
                            continue
                        insilico[x] = insilico[x]+[peptide[x-i[0]]]

            else:#means that it will be aligned but with insertions or deletions
                positions = [(match.start(),match.end()) for match in re.finditer(peptide, str(og_seq))]
                for nr_position,i in enumerate(positions):
                    add_insertions = 0
                    while i[0]+add_insertions <len(consensus_aligned_og) and consensus_aligned_og[:i[0]+add_insertions].replace('-','').count(peptide)!=nr_position+1:#locate AA position1
                        add_insertions+=1
                    insert_adjust = False
                    for minus in range(1,i[0]+add_insertions):
                        if consensus_aligned_og[:i[0]+add_insertions][-minus:].replace('-','') == peptide:
                            add_insertions -= minus
                            insert_adjust = True
                            break
                    if insert_adjust==False:
                        print('Could not locate insert for {}'.format(peptide))
                        continue
                    done.append(peptide)#not 2 times same 
                    temp_insert = add_insertions
                    for x in range(i[0],i[1]):
                        if x+add_insertions>=len(consensus_aligned_og):
                            continue
                        while add_insertions+1 != i[1]+add_insertions and consensus_aligned_og[x+add_insertions]=='-' and x+add_insertions+1<len(consensus_aligned_og):
                            add_insertions+=1
                        if peptide[x-i[0]] in insilico[x+add_insertions]:
                            continue
                        insilico[x+add_insertions] = insilico[x+add_insertions]+[peptide[x-i[0]]]
                    print('Insertion found: {} became {}'.format(peptide,consensus_aligned_og[temp_insert+i[0]:add_insertions+i[1]]))
    return insilico

def find_search_space(seq_to_anim,animals_to_include,seqs,anim,mammals,df,taxa_name,protein_name,taxonomy,path,sample_path):
    calculator = DistanceCalculator('blosum90')
    sub_matrix = substitution_matrices.load('BLOSUM90')#substitution_matrices.load()
    sub_matrix[22]=[-1]*22+[5]+[-1]#X residue
    sub_matrix[:,22]=[-1]*22+[5]+[-1]
    sub_matrix[21]=[-10]*21+[10]+[-10,-10]#Z residue
    sub_matrix[:,21]=[-10]*21+[10]+[-10,-10]
    
    sequences = {v:k for k,v in seqs.items()}
    seqs_4_align = {val:sequences[key] for key,val in seq_to_anim.items()}
    t_tree = find_taxonomic_tree(anim,mammals)
    consensus_tree = {}

    for i in sorted(t_tree.keys()):
        level = t_tree[i]
        for group,species in level.items():
            species = sorted(species)
            if len(species)<2 or tuple(species) in consensus_tree:
                continue #no alignment possible
            in_consensus_tree = [el for el in consensus_tree.keys() if len(set(el)&set(species))!=0]
            not_in_consensus_tree = list(set([el for el in species if any(el in k for k in in_consensus_tree)==False]))
            
            longest = []#FIND MAXIMAL LEVEL FOR EACH OF THE SPECIES
            for l in sorted(in_consensus_tree, key=lambda x:len(x)):
                for x in longest:
                    if set(x).issubset(set(l))==True:
                        longest = [el for el in longest if el!=x]
                longest.append(l)
            in_consensus_tree = longest
            if len(in_consensus_tree) != 0:
                #Yes consensus = add to existing consensus, even add consensus to consensus
                consensus = consensus_tree[in_consensus_tree[0]]
                in_consensus_tree = [el for el in in_consensus_tree if el != in_consensus_tree[0]]
                while len(in_consensus_tree)>0:#first align consensus seqs
                    next_seq = consensus_tree[in_consensus_tree[0]]
                    consensus,old_c,old_s = multiple_alignment(consensus,next_seq,calculator,sub_matrix)
                    consensus_tree[in_consensus_tree[0]]=old_s
                    in_consensus_tree = [el for el in in_consensus_tree if el != in_consensus_tree[0]]
                    
                
                while len(not_in_consensus_tree)>0:#align the new ones
                    next_seq = seqs_4_align[not_in_consensus_tree[0]]
                    consensus,old_c,old_s = multiple_alignment(consensus,next_seq,calculator,sub_matrix)
                    seqs_4_align[not_in_consensus_tree[0]]=old_s
                    not_in_consensus_tree = [el for el in not_in_consensus_tree if el != not_in_consensus_tree[0]]
                    
            else:
                #No consensus = make consensus and save
                #1. take the longest sequence
                longest = ['',0]
                for s in species:
                    if longest[1]<len(seqs_4_align[s]):
                        longest = [s,len(seqs_4_align[s])]
                #2. take next longest sequence and align to the longest one
                consensus = seqs_4_align[longest[0]]
                not_in_consensus_tree = [el for el in not_in_consensus_tree]# if el!=longest[0]]
                while len(not_in_consensus_tree)>0:
                    next_seq = seqs_4_align[not_in_consensus_tree[0]]
                    
                    consensus,old_c,old_s = multiple_alignment(consensus,next_seq,calculator,sub_matrix)
                    seqs_4_align[not_in_consensus_tree[0]]=old_s
                    # seqs_4_align[longest[0]]=old_c
                    not_in_consensus_tree = [el for el in not_in_consensus_tree if el != not_in_consensus_tree[0]]
            #save the new consensus sequence
            consensus_tree[tuple(species)]=str(consensus)

    cols = consensus_tree[sorted(list(consensus_tree.keys()), key=lambda x:len(x))[-1]]
    ind = []
    matrix = []
    unique_taxa = []
    tree_matrix = sub_matrix
    tree_matrix[22]=[0]*22+[5*cols.count('X')/20]+[0]#X residue
    tree_matrix[:,22]=[0]*22+[5*cols.count('X')/20]+[0]
    for k,v in consensus_tree.items():
        LCA = find_LCA(taxonomy,k)
        unique_taxa.append(LCA)
        ind.append(LCA+'_'+str(len(k)))
        if len(v)!=len(cols):
            alignment = aligning(cols,v,tree_matrix)
            add = alignment[1]
        else:
            add=v
        matrix.append(list(add))
    for k,v in seqs_4_align.items():
        smallest_cols = sorted([el for el in consensus_tree.keys() if k in el],key=lambda x:len(x))[0]
        LCA = find_LCA(taxonomy,smallest_cols)
        LCA=LCA+'_'+str(len(smallest_cols))
        LCA=ind.index(LCA)
        ind.append(k)
        s_cols = ''.join(matrix[LCA]).replace('-','X')
        alignment = aligning(s_cols,v,tree_matrix)
        add = alignment[1]
        if len(add)<len(cols):
            add+='-'*(len(cols)-len(add))
        while len(add)>len(cols):
            add = add[:-1]
        matrix.append(list(add))
    remove = []
    for i in set(unique_taxa):
        counter = []
        for q in ind:
            if '_' in q and i in q:
                counter.append(int(q.split('_')[-1]))
        if len(counter)>0:
            remove = remove + [i+'_'+str(el) for el in counter if el!= max(counter)]
    matrix = [el for num,el in enumerate(matrix) if ind[num] not in remove]
    ind = [el.split('_')[0] for el in ind if el not in remove]
    consensus_df = pd.DataFrame(matrix,columns=list(cols),index=ind)
    
    consensus_df = consensus_df[~consensus_df.index.duplicated(keep='first')]
    keep = []
    for i in range(0,len(consensus_df.columns)):
        if list(consensus_df.iloc[:,i].values).count('-')<len(consensus_df.iloc[:,i].values)*0.1 and set(consensus_df.iloc[:,i].values)!={'-','X'}:
            keep.append(i)
    consensus_df=consensus_df.iloc[:,keep]
    
    
    higher_taxa = [el for el in ind if el not in seqs_4_align.keys()]
    belongs = {}
    for i in higher_taxa:
        belongs[i]=[]
        for x in seqs_4_align.keys():
            if any(i in el for el in mammals[x])==True:
                belongs[i]=belongs[i]+[x]
    print('Turning it into numbers')
    consensus_in_numbers = {}
    for i in higher_taxa:
        consensus = ''.join(list(consensus_df.loc[i].values))
        number_seq = ''
        for num,el in enumerate(consensus):
            if el == 'X':
                count = []
                for x in belongs[i]:
                    temp = ''.join(list(consensus_df.loc[x].values))
                    count.append(temp[num])
                number_seq += str(len(set(count)))
            else:
                number_seq += '1'
        consensus_in_numbers[i]=number_seq
    
    #import the results
    df=df[df['Species'].isin(consensus_df.index)]
    mapped = mapping(seqs_4_align,consensus_df,df)
             
    #make plot to visualize mapping
    fig = make_subplots(rows=len(higher_taxa), cols=1,
                        subplot_titles=higher_taxa,shared_xaxes=True,shared_yaxes=True)
    count = 1
    for k in higher_taxa:
        v = consensus_in_numbers[k]
        # overlap = [1 if el=='X' else max([int(el) for el in v])+1 for el in insilico[key]]
        fig.add_trace(go.Scatter(x=[el for el in range(0,len(v))], y=[int(el) for el in v],
                            mode='lines',fill=None,
                            name=k),row=count,col=1)
        count += 1
    if 'PREDICTED' in protein_name:
        protein_name=protein_name.split('PREDICTED')[0]
    protein_name = protein_name.replace(':','')
    fig.update_layout(title_text='Mutations found within '+taxa_name+' matching to proteins grouping with '+protein_name)
    
    try:
        name_file = taxa_name+'_multiple_alignment_'+protein_name+'.html'
        fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'aligned' /name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'aligned' / 'multiple_alignment.html')
    
    return consensus_df,mapped

def do_consensus(df, sequences, taxonomy, taxonomic_groups,restrict,path,sample_path,animals_in_input):
    total_consensus_df = {}
    total_insilico = {}
    seperate_outcome = []
    all_sequences_considered = []
    reverse_sequences={v:k for k,v in sequences.items()}
    considered_peptides_after_filtering = set()
    #from here split the proteins accorind to the input 
    for taxa in taxonomic_groups.keys():
        mammals = {k:v for k,v in taxonomy.items() if any(taxa in el for el in v) and k in animals_in_input}
        animals_to_include = [el for el in mammals.keys()]
        #group per protein
        temp = df[['Peptides','Protein','Species','Original_score']][df['Species'].isin(animals_to_include)]
        max_score = max(list(temp['Original_score'].values))
        max_score_species = list(set(temp['Species'][temp['Original_score']==max_score].values))
        if len(max_score_species)>1:#1 representative
            max_score_species = sorted([(el,list(set(temp['Protein'][temp['Species']==el].values))) for el in max_score_species], key=lambda x:x[1])[::-1][0][0]
        else:
            max_score_species =max_score_species[0]
        proteins = list(set(temp['Protein'][temp['Species']==max_score_species].values))
        protein_to_peptides={}
        peptides_included = set()
        peptides_included_highest = set()
        proteins = sorted([(el,len(reverse_sequences[el])) for el in proteins], key=lambda x:x[1])[::-1]
        proteins = [el[0] for el in proteins]
        for x in proteins:
            peptides=set(temp['Peptides'][temp['Protein']==x].values)
            isoform = False
            for already_added_pr, already_added_vals in protein_to_peptides.items():
                if len(peptides^set(already_added_vals))<=0.1*len(peptides):
                    print(f'!!Removing {x}, because it was found to be an isoform from another that will be used for multiple-alignment based on peptide content!!')
                    isoform = True
                    break
            if isoform==True:
                continue
            coverage = [0]*len(reverse_sequences[x])
            seq = str(reverse_sequences[x])
            for p in peptides:
                if p in seq:
                    ind = seq.index(p)
                    for r in range(0,len(p)):
                        coverage[ind+r]=1
            if len([el for el in coverage if el==1])<=len(seq)*0.25: #we need at least 25% coverage
                print(f'Coverage of {x} is not enough for mixture analysis, minimum of 25% required!')
                continue
            peptides_diff = peptides_included^peptides
            peptides_included = peptides_included|peptides
            peptides = peptides&peptides_diff
            peptides_included_highest = peptides_included_highest|set(peptides)
            protein_to_peptides[x]=peptides
        print('Doing {}'.format(taxa))
        
        proteins_all = list(set(temp['Protein'].values))
        for x in proteins_all:
            peptides=set(temp['Peptides'][temp['Protein']==x].values)
            coverage = [0]*len(reverse_sequences[x])
            seq = str(reverse_sequences[x])
            for p in peptides:
                if p in seq:
                    ind = seq.index(p)
                    for r in range(0,len(p)):
                        coverage[ind+r]=1
            if len([el for el in coverage if el==1])<=len(seq)*0.2 or len(set(peptides)&set(peptides_included_highest))<len(peptides)/2: #we need at least 20% coverage
                continue
            considered_peptides_after_filtering = considered_peptides_after_filtering|set(peptides)
        
        representatives = []
        for a in set(df['Species'].values):
            lin=taxonomy[a]
            lca = []
            for b in lin:
                if b[1]==taxa or b[0]==restrict:
                    break
                lca.append(b[1])
                if 'order' in b[0]:
                    break#we restrict here
            representatives = representatives+lca
        temp_df = {}
        anims = []
        for protein,peps in protein_to_peptides.items():
            length_protein = len(''.join([str(key) for key,val in sequences.items() if val==protein]))
            # if len(peps)<=0.02*len(reverse_sequences[protein]):
            #     continue
            print('Start on taxon: {} matching protein groups of {}'.format(taxa,protein))
            seqs = {}
            anim = []
            linked_side_lin = []
            for a in animals_to_include:
                if any(el[1] in representatives for el in taxonomy[a])==False:#if for example only pecora species, add suidae and camelidae too
                    if len([el[1] for el in taxonomy[a] if el[1] in linked_side_lin])>2:#3 reps per side lineage if possible to have some more potential mutational zones
                        continue
                    side = []
                    for el in taxonomy[a]:
                        if 'order' in el:
                            break
                        side.append(el[1])
                    linked_side_lin = linked_side_lin + side
                temp = []
                for k,v in sequences.items():
                    sp=animals_from_db_input_mix({k:v})
                    if sp[0]!=a:
                        continue
                    temp.append([v,len([el for el in peps if el in str(k)])])
                if len(temp)==0:
                    continue
                
                temp = sorted(temp,key=lambda x:x[1])
                largest = temp[-1][1]
                temp = sorted([el[0] for el in temp if el[1]==largest])
                if largest<=5:
                    continue
                if largest<=len(peps)/4 and len(peps)>100 and a not in set(df['Species'].values):#If the true sequence is missing
                    continue
                if len(str(reverse_sequences[temp[0]]))<length_protein*0.9:#means protein too trunctated to consider or a totally different one
                    continue                
                anim.append(a)
                seqs[str(reverse_sequences[temp[0]])]=temp[0]
            temp_df[protein]=seqs
            anims.append(anim)
        if any(len(el)<2 for el in anims) or len(anims)==0:
            continue
        anim = set(anims[0])
        for an in anims:
            anim = anim&set(an)
        anim=list(anim)
        print('Performing multiple alignment on {}'.format(anim))
        for protein,seqs in temp_df.items():
            animals_to_include_now = anim
            seq_to_anim = {}
            for i in seqs.values():
                for a in anim:
                    if a in i:
                        seq_to_anim[i]=a
                        all_sequences_considered.append(i)
                        break
            animals_to_include_now = [el for el in set(seq_to_anim.values())]
            anim = [el for el in set(seq_to_anim.values())]
            if len(anim)<2 or len(animals_to_include_now)<2:
                seperate_outcome.append(anim[0])
                continue
            consensus_df,insilico = find_search_space(seq_to_anim,animals_to_include_now,seqs,anim,mammals,df,taxa,protein,taxonomy,path,sample_path)
            total_consensus_df[(taxa,protein)] = consensus_df
            total_insilico[(taxa,protein)] = insilico
            time.sleep(2)
    return total_consensus_df,total_insilico,seperate_outcome,all_sequences_considered,considered_peptides_after_filtering

def delete_groups(taxonomic_groups, df,path,sample_path,file_extinct,considered_peptides_after_filtering):
    removing =True
    removed = []
    colour_add = []
    name_add = []
    u_add = []
    u_per_sp = []
    df = df[df['Peptides'].isin(considered_peptides_after_filtering)]
    include_order = 0
    include_order_per_species = 0
    while removing ==True:
        order_uniqueness = []
        names_groups = []
        all_in_groups = []
        order_unique_per_species = []
        for key,n in taxonomic_groups.items():
            if key in removed:
                continue
            all_in_groups = all_in_groups+n
        df = df[df['Species'].isin(all_in_groups)] 
        for group,species in taxonomic_groups.items():
            if group in removed:
                continue
            
            names_groups.append(group)
            order_peps = set(df['Peptides'][df['Species'].isin(species)].values)
            other_peps = set(df['Peptides'][df['Species'].isin(species)==False].values)
            diff = set(order_peps)^set(other_peps)
            order_uniqueness.append(len(set(order_peps)&diff))
            order_unique_per_species.append(len(set(order_peps)&diff)/len(species))
        if True not in [False if el>=max(2,0.1*max(order_uniqueness)) and order_unique_per_species[loc]>=1.5 else True for loc,el in enumerate(order_uniqueness)] or len(order_uniqueness)<0:
            removing = False
            include_order=min(order_uniqueness)
            include_order_per_species=min(order_unique_per_species)
        elif len(order_uniqueness)>0:
            # r = [(names_groups[loc],el) for loc,el in enumerate(order_uniqueness) if el<max(2,0.1*max(order_uniqueness)) or order_unique_per_species[loc]<2]
            r = [(names_groups[loc],el,order_uniqueness[loc]) for loc,el in enumerate(order_unique_per_species) if el<1.5 or order_uniqueness[loc]<max(2,0.1*max(order_uniqueness))]
            if len(r)==0 and min(order_uniqueness)>2:
                removing = False
                include_order=min(order_uniqueness)
                break
            r_num = sorted(r,key=lambda x:x[2])
            r_num = r_num[:3]#decide which one of the lowest 3
            r_num = sorted(r_num,key=lambda x:x[1])
            r_num_order =r_num[0][2]
            r_num = r_num[0][1]
            
            add = [el[0] for el in r if el[1]==r_num and el[2]==r_num_order]
            colour_add = colour_add + ['crimson']*len(add)
            name_add = name_add + add
            u_add = u_add + [r_num_order]*len(add)
            u_per_sp = u_per_sp + [r_num]*len(add)
            removed = removed + add
            include_order=min(order_uniqueness)
            include_order_per_species=min(order_unique_per_species)
        
        
    colours = ['forestgreen' if el>=max(2,include_order) and order_unique_per_species[loc]>=max(1.5,include_order_per_species)  else 'crimson' for loc,el in enumerate(order_uniqueness)]+colour_add[::-1]
    fig = go.Figure()
    fig.add_trace(go.Bar(x=names_groups+name_add[::-1], y=order_uniqueness+u_add[::-1],
                    marker_color=colours,hovertext=[f'{element:.2f} unique peptides per species at order level' for element in order_unique_per_species+u_per_sp[::-1]],
                    name='Order uniqueness relative to number of candidate species'))
    fig.update_layout(title=dict(text="Sequential unique peptides at order level {} relative to species".format(file_extinct)),
                      yaxis=dict(
                            title=dict(
                                text="Unique peptide count",
                                font=dict(
                                    size=16
                                )
                            ),
                        ),
                      xaxis=dict(
                            title=dict(
                                text="Order (red=removed from possibilities)",
                                font=dict(
                                    size=16
                                )
                            ),
                            categoryarray=list(np.array(names_groups)[np.argsort(np.array(order_uniqueness))[::-1]]),
                        ),
                      )
    name_file = 'Before_start_order_uniqueness_'+file_extinct+'.html'
    
    try:
        fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'taxonomic_output' /name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'taxonomic_output' / 'Before_start_uniqueness.html')    
    
    return [el for loc,el in enumerate(names_groups) if order_uniqueness[loc]>=max(2,include_order) and order_unique_per_species[loc]>=max(1.5,include_order_per_species)]#At least 2 unique peptides to be considered. More is there is a high abundance of uniqueness (high detectability, low chance of having missed the unique ones of other species less abundant in the sample)


def score_directionallity(num,df_c,AA,lineage1,lca,lineage2,df,db_species,taxonomy):
    #adapt here so that if subspecies on both sides the analysis includes subspecies, else species is enough to lower the bias
    l_r_path = lineage1+[lca]+lineage2[::-1]
    slice_df = df_c.iloc[:,num]
    species1 = lineage1[0]
    species2 = lineage2[0]
    #Uniqueness_score, probability that the residue is unique under the LCA
    species_under_taxon = {el:list(background_check_mix([el],[],taxonomy).values())[0] for el in l_r_path}
                           # if ' ' not in el}#only higher levels
    species_linked_residue = list(slice_df.index[slice_df==AA])
    sut = species_under_taxon[lca]
    AA_in_db = len([el for el in species_linked_residue if el in db_species and el in sut])
    species_in_db = len([el for el in sut if el in db_species])
    if species_in_db==0:
        species_in_db += 1
    missing = len([el for el in sut if el not in db_species])
    species_under_lca = len(sut)
    if species_under_lca==0:
        species_under_lca += 1
    P_unique = 1-((AA_in_db/species_in_db)+(1-missing/species_under_lca))
    
    #mutation location in the tree   
    #determine lineage where the mutation is found
    add_to_calc = []
    if species1 in species_linked_residue:
        track = species1
    else:
        track=species2
    if 'theoretical' in species1:
        add_to_calc.append(species1)
    if 'theoretical' in species2:
        add_to_calc.append(species2)
    #if both have the same amino acide residue than the score should be the same
    t = taxonomy[track]
    x=0
    Precent = 0
    while lca not in t[x]:
        back = background_check_mix([t[x][1]],[],taxonomy)[t[x][1]]+add_to_calc
        sp_in_db_with_AA = len([el for el in back if el in species_linked_residue])
        back_in_db = len([el for el in back if el in slice_df.index])+1
        missing = len([el for el in back if el not in slice_df.index])
        if back_in_db==0:
            back_in_db=1
        mutationalChance = (1-(sp_in_db_with_AA/back_in_db))*(1-(missing/(len(back)+1)))
        Precent += mutationalChance
        x+=1
    if x==0:
        x=1
    residue_score = (P_unique+(Precent/x))
    return residue_score

def look_for_complementarity(compare,df_c,df,species_in,protein_name,taxonomy,insilico):
    ancestor = find_LCA(taxonomy,list(compare))
    c1 = taxonomy[compare[0]]
    twodim = []
    for x in c1:
        if ancestor not in x:
            twodim.append(x[1])
        else:
            break
    twodimreverse = []
    c2 = taxonomy[compare[1]]
    for x in c2:
        if ancestor not in x:
            twodimreverse.append(x[1])
        else:
            break
    if ancestor==compare[0]:
        itself=True
        ancestor = [el[1] for el in taxonomy[compare[0]] if el[1] in df_c.index and el[1]!= compare[0]][0]
    else:
        itself=False
    split_lineage = twodim+[ancestor]+twodimreverse[::-1]
    plot_data = {"line_x": [], "line_y": [], "Baseline": [], "Taxon_direction": [],"weight":[],'residue':[]}
    pep_locs = np.array([el+'_'+str(num) for num,el in enumerate(df_c.columns)])#LCA = the X list, otherwise at large scale way too slow
    pep_locs_plot = list(pep_locs)
    count_complement = 0
    temp_lca = ancestor
    while temp_lca not in df_c.index:
        temp_lca = find_LCA(taxonomy,[temp_lca,compare[0]],True)
    temp_compare = df_c.loc[temp_lca].values
    compare_VS1 = df_c.loc[compare[0]].values
    compare_VS2 = df_c.loc[compare[1]].values
    temp_VS = df_c[df_c.index.isin(split_lineage)]
    reduce_complete_peptides = 0#due to missingness in the sequences sometimes a complete peptide can be matched, however this will bias towards genetic mixes or towards the more complete sequence, we reduce these sequences to 1 representative residue
    for num,pl in enumerate(pep_locs):
        if 'X_' in pl and len(insilico[num])!=1 and temp_compare[num]=='X' and compare_VS1[num]!=compare_VS2[num]:#check if the location is found and mutational
            #find matching taxon   
            temp = temp_VS.iloc[:,num]
            for c,AA in enumerate(insilico[num]):
                if AA == 'X':
                    continue
                including = list(temp.index[temp==AA])
                if len(including)==0:
                    direction=ancestor
                    direction_weight = 0
                elif num==reduce_complete_peptides+1:
                    direction=ancestor
                    direction_weight = 1
                else:
                    direction=find_LCA(taxonomy,including)
                    direction_weight = score_directionallity(num,df_c,AA,twodim, ancestor,twodimreverse,df,species_in,taxonomy)
                if c!=1:
                    pep_locs_plot = pep_locs_plot[:num+count_complement]+[pl]+pep_locs_plot[num+count_complement:]
                    count_complement += 1
                plot_data["line_y"].extend([pl,pl, None])#for arrows
                plot_data["line_x"].extend([ancestor,direction,None])#change pl to taxon score#for arrows
                plot_data["Baseline"].extend([ancestor])#sequence location in plot
                plot_data['Taxon_direction'].extend([direction])#taxon direction score
                plot_data["weight"].extend([direction_weight])#add a weight to mutational residues
                plot_data["residue"].extend([AA])
            reduce_complete_peptides = num
        else:#non-relevant locations
            plot_data["line_y"].extend([pl,pl, None])#for arrows
            plot_data["line_x"].extend([ancestor,None,None])#change pl to taxon score#for arrows
            plot_data["Baseline"].extend([ancestor])#sequence location in plot
            plot_data['Taxon_direction'].extend([ancestor])#taxon direction score
            if 'X_' in pl and len(insilico[num])!=1 and temp_compare[num]=='X' and itself==True:
                plot_data["weight"].extend([1])#if detected and mutational space but we keep it apart for next round. The location might influence other outcomes
            else:
                plot_data["weight"].extend([0])#No weight needed for non-unique residues
            plot_data["residue"].extend([''])
    plot_data['weight']=[el/max(list(plot_data['weight'])) if max(list(plot_data['weight'])) != 0 
                                                                  else el for el in plot_data['weight']] #min-max normalize
    plot_data['pep_locs_plot']=pep_locs_plot
    
    return [compare[0] +'_VS_'+compare[1],plot_data]

def initial_iteration(total_consensus_df, total_insilico,df,taxonomy,cpu_count,species_in):
    pep_locs_all = {}
    mixtures_all = {}
    with ProcessPoolExecutor(max_workers=cpu_count) as executor2:
        with tqdm(list(total_consensus_df), desc="Starting on comparison...") as pbar:
            for protein_al in pbar:
                pbar.set_description(f"Performing 1-VS-1 comparison on {protein_al}")
                df_c = total_consensus_df[protein_al]
                
                #map all peptides found to the consensus
                #link residues to species, only 'x' residues
                insilico= total_insilico[protein_al]
                output_species = [el for el in df_c.index if el in set(df['Species'].values)]
                out_species = []
                if len(output_species)!=1:
                    for el in output_species:
                        for num in output_species:
                            if el!=num and sorted([el,num]) not in out_species:
                                out_species.append(sorted([el,num]))
                else:
                    out_species = [[el,el] for el in output_species]
                found = []
                test = [el for el in out_species]
                adjust_theoretical = []#No need to combine sequences if a lower taxonomic combination can be made
                while len(found) != len(output_species) and len(found+output_species)>1:     
                    at = []
                    for combos in test:
                        at.append(find_LCA(taxonomy,combos))
                    for combos in test:
                        add_combo = True
                        c1 = taxonomy[combos[0]]
                        c2 = taxonomy[combos[1]]
                        for x in c1:
                            if x[1] in at:
                                most_recent = x[1]
                                break
                        for x in c2:
                            if x[1] in at:
                                if x[1]!= most_recent:
                                    adjust_theoretical.append(combos)
                                    add_combo = False
                                break
                        if add_combo == True:
                            found = list(set(found + combos))
                    test = [el for el in test if el[0] not in found and el[1] not in found]
                    singles = []
                    for el in test:
                        singles = singles+el
                    singles = list(set(singles))
                    # adjust_theoretical = [el for el in adjust_theoretical if el not in test]
                    for el in singles:
                        found.append(el)
                        out_species.append([el,el])
                    test= []
                    
                    # if len(test)==0 and len(output_species)!= len(found):
                    # for x in adjust_theoretical:
                    #     if x[0] not in found or x[1] not in found:
                    #         adjust_theoretical = [el for el in adjust_theoretical if el != x]
                    #         found = list(set(found + x))
                    #         break
                    for x in adjust_theoretical:
                        if x[0] not in found or x[1] not in found:
                            # adjust_theoretical = [el for el in adjust_theoretical if el != x]
                            for t in x:
                                if any(t in element for element in out_species if element not in adjust_theoretical)==False:
                                    out_species.append([t,t])
                            found = list(set(found + x))
                            break
                out_species = [el for el in out_species if el not in adjust_theoretical]
                
                try:
                    results = []
                    futures = [executor2.submit(look_for_complementarity,compare,df_c,df,species_in,protein_al,taxonomy,insilico) \
                            for compare in out_species]
                    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing", leave=False):
                        results.append(future.result())
                    mixtures = {key:val for key, val in results}
                    pep_locs_all[protein_al] = np.array([el+'_'+str(num) for num,el in enumerate(df_c.columns) if 'X' in el])
                    mixtures_all[protein_al]=mixtures
                except Exception as e:
                    pbar.write(f"Error processing {protein_al}: {e}")
                finally:
                    pbar.update(1)
    return pep_locs_all, mixtures_all

def find_trace(t,tb):
    if 'theoretical' not in t:
        return [t]
    trace = []
    for x in tb[int(t.split('_')[-1])]:
        trace = trace + find_trace(x,tb)
    return trace

def find_similar_resiudes(search_space,data,sp_data,tax,theoretical_trace_back,taxonomy,acomp):#search_space,unique_left_not_in_others,VSleft,tax,theoretical_trace_back
    tax = [el[1] for el in tax]
    tax_search = max([loc for loc,el in enumerate(tax) if el in search_space])
    tax = tax[:tax_search+1]
    tax_trace = []
    for anim in find_trace(tax[0],theoretical_trace_back):#bring taxonomy together from theoretical
        tax = tax + [el[1] for el in taxonomy[anim] if el in search_space+tax[:tax_search+1]]
        max_match = [loc for loc,el in enumerate(taxonomy[anim]) if el[1] in tax]
        if len(max_match)>0:
            max_match = max(max_match)
            tax_trace = tax_trace + [el[1] for loc,el in enumerate(taxonomy[anim]) if loc <=max_match]
    #if combination is made with the species it may not account for uniqueness, because it was donated by the original
    tax_trace=tax_trace+['theoretical_'+str(key) for key,val in theoretical_trace_back.items() if any(tax[0] in q for q in find_trace(val[0],theoretical_trace_back)+find_trace(val[1],theoretical_trace_back))]
    for a in acomp:
        if tax[0] in a and len(set(a))>1:
            parallel = find_trace([el for el in a if el!=tax[0]][0],theoretical_trace_back)
            for p in parallel:
                if p not in taxonomy:
                    continue
                p = [par[1] for par in taxonomy[p] if par[1] in search_space]
                tax_trace = tax_trace + p#if compared and survived, the residues can belong to either so they remain deterministic in the comparison
            
    for a in acomp:#compared in parallel, meaning that if it survived the individuals than it must survive the combined ones
        if tax[0] not in a:
            tr = find_trace(a[0],theoretical_trace_back)+find_trace(a[1],theoretical_trace_back)
            if any(element not in tax_trace for element in tr)==False:
                parallel = ['theoretical_'+str(k) for k,v in theoretical_trace_back.items() if sorted(v)==sorted(a)]
                if len(parallel)>0:
                    for p in parallel:
                        if p not in taxonomy:
                            continue
                        p = [par[1] for par in taxonomy[p] if par[1] in search_space]
                        tax_trace = tax_trace+p
    tax_trace=list(set(tax_trace))
    keep_res = []
    sp_data=sp_data.T
    sp_data = sp_data[sp_data['taxon'].isin(tax)].T
    for q in set(data.columns):
        temp_d = data.T
        temp_sp = sp_data[q]
        temp_sp=temp_sp.loc['residue']
        if len(temp_sp)>1:
            temp_sp=temp_sp.values
            temp_sp=list(temp_sp)[0]
        temp = temp_d[(temp_d['taxon'].isin(tax_trace)==False) & (temp_d['residue']==temp_sp)].T
        if q not in temp:
            keep_res.append(q)
            continue
        temp = temp[q]
        temp = temp.loc['residue']
        if len(temp)>1:
            temp=temp.values
        if temp_sp not in temp:
            keep_res.append(q)
    keep_res = [el for el in set(data.columns) if el not in keep_res] #keep the ones with another explanation
    return keep_res

def find_potential_missing_species(search_space,mix,peps,taxonomy,theoretical_trace_back,acomp):#search_space,mixtures_all,pep_locs_all,taxonomy,theoretical_trace_back
    out = {}
    combos = []
    for val in mix.values():
        combos = combos + list(val.keys())
    combos = set(combos)
    
    weights = []
    taxon_dir=[]
    pep_locs= []
    filter_pep = []
    residue = []
    acomp = [el for el in acomp if any(element.split('_VS_') == el for element in combos)==False]
    for i in combos:
        print('including {}'.format(i))
        for num,protein in enumerate(list(mix.keys())):
            temp = mix[protein]
            pep_locs = pep_locs+[el+'_'+str(num) for el in temp[i]['pep_locs_plot']]
            filter_pep = filter_pep+[el+'_'+str(num) for el in peps[protein]]
            weights = weights + temp[i]['weight']
            taxon_dir = taxon_dir + temp[i]['Taxon_direction']
            residue = residue + temp[i]['residue']
    
    double_df_all = pd.DataFrame([weights,taxon_dir,residue],
                          columns=pep_locs,index=['weight','taxon','residue'])
    double_df_all = double_df_all[[el for el in set(double_df_all.columns) if 'X' in el]]
    double_df_all = double_df_all.T
    double_df_all = double_df_all[double_df_all['weight']>0]
    double_df_all['weight'] = [0]*len(double_df_all)
    double_df_all['ind']=double_df_all.index
    double_df_all.index = [r for r in range(0,len(double_df_all))]
    double_df_all = double_df_all.drop_duplicates()
    double_df_all.index = double_df_all['ind']
    double_df_all=double_df_all[['weight','taxon','residue']]
    double_df_all = double_df_all.T
    double_df_all = double_df_all[[el for el in set(double_df_all.columns) if list(double_df_all.columns).count(el)>1]]

    for i in combos:
        left = i.split('_VS_')[0]
        right = i.split('_VS_')[-1]
        lca = find_LCA(taxonomy,[left,right])
        left = [el[1] for el in taxonomy[left]]
        left = left[:left.index(lca)]
        right = [el[1] for el in taxonomy[right]]
        right = right[:right.index(lca)]
        weights = []
        taxon_dir=[]
        pep_locs= []
        filter_pep = []
        residue = []
        for num,protein in enumerate(list(mix.keys())):
            temp = mix[protein]
            pep_locs = pep_locs+[el+'_'+str(num) for el in temp[i]['pep_locs_plot']]
            filter_pep = filter_pep+[el+'_'+str(num) for el in peps[protein]]
            weights = weights + temp[i]['weight']
            taxon_dir = taxon_dir + temp[i]['Taxon_direction']
            residue = residue + temp[i]['residue']
        oneVSone=pd.DataFrame([weights,taxon_dir,residue],
                              columns=pep_locs,index=['weight','taxon','residue'])

        oneVSone = oneVSone[filter_pep]
        oneVSone.loc['weight']=oneVSone.loc['weight'].apply(lambda x:float(x))
        oneVSone = oneVSone.T
        oneVSone = oneVSone[oneVSone['weight']>0].T
        
        oneVSone_noLCA = oneVSone.loc[:, ~(oneVSone == lca).any()]#remove LCA in 1VS1 comparison
        
        physical_double = [el for el in set(oneVSone.columns) if list(oneVSone.columns).count(el)>1]#overlap, both side have evidence
        complete_mix = [el for el in set(oneVSone_noLCA.columns) if list(oneVSone.columns).count(el)==1]#complement, unique for 1 side
        
        if len(physical_double)>0:
            VSleft = oneVSone_noLCA.T
            VSleft = VSleft[(VSleft['taxon'].isin(left))]
            VSleft = VSleft.T
            
            VSr = oneVSone_noLCA.T
            VSr = VSr[(VSr['taxon'].isin(right))]
            VSr = VSr.T
            complete_mix_left = [el for el in set(VSleft.columns) if list(oneVSone.columns).count(el)==1]
            complete_mix_right = [el for el in set(VSr.columns) if list(oneVSone.columns).count(el)==1]
            tm= False
            
            #Can the complementary locational residue be found in other candidates
            unique_left_not_in_others = double_df_all[[element for element in complete_mix_left if element in double_df_all.columns]].T
            unique_left_not_in_others = unique_left_not_in_others[unique_left_not_in_others['taxon'].isin(left+[lca]+right)==False].T
            
            tax = taxonomy[left[0]]
            del_residue = find_similar_resiudes(search_space,unique_left_not_in_others,VSleft,tax,theoretical_trace_back,taxonomy,acomp)
            keep_residue_cml = [el for el in complete_mix_left if el not in del_residue]
            
            tax=taxonomy[right[0]]
            unique_right_not_in_others = double_df_all[[element for element in complete_mix_right if element in double_df_all.columns]].T
            unique_right_not_in_others = unique_right_not_in_others[unique_right_not_in_others['taxon'].isin(left+[lca]+right)==False].T
            del_residue = find_similar_resiudes(search_space,unique_right_not_in_others,VSr,tax,theoretical_trace_back,taxonomy,acomp)
            keep_residue_cmr = [el for el in complete_mix_right if el not in del_residue]
            
            #Can te overlap locations be found in other candidates. If yes than non informative location
            uol = [el for el in set(VSleft.columns) if list(oneVSone.columns).count(el)>1]
            uor = [el for el in set(VSr.columns) if list(oneVSone.columns).count(el)>1]
            unique_overlap_left = double_df_all[[element for element in uol if element in double_df_all.columns]].T
            unique_overlap_left = unique_overlap_left[unique_overlap_left['taxon'].isin(left+[lca]+right)==False].T
            tax=taxonomy[left[0]]
            del_residue = find_similar_resiudes(search_space,unique_overlap_left,VSleft,tax,theoretical_trace_back,taxonomy,acomp)
            keep_residue_uol = [el for el in uol if el not in del_residue]
            
            
            unique_overlap_right = double_df_all[[element for element in uor if element in double_df_all.columns]].T
            unique_overlap_right = unique_overlap_right[unique_overlap_right['taxon'].isin(left+[lca]+right)==False].T
            tax=taxonomy[right[0]]
            del_residue = find_similar_resiudes(search_space,unique_overlap_right,VSr,tax,theoretical_trace_back,taxonomy,acomp)
            keep_residue_uor = [el for el in uor if el not in del_residue]
            
            #is there still uniqueness on the left side
            det_pos_1 = uol+complete_mix_left
            common_pos_1 = keep_residue_cml+keep_residue_uol
            det_pos_2 = uor+complete_mix_right
            common_pos_2 = keep_residue_cmr+keep_residue_uor
            #need for multiple doubles and they have to be true doubles and in total overlap needs to be higher than 10%
            if len(set(det_pos_1)&set(common_pos_1))>1 and len(physical_double)>1 and (len(set(det_pos_1)&set(common_pos_1))/(len(set(det_pos_1))+len(set(det_pos_2)))>0.25 or len(set(det_pos_1))/len(set(det_pos_2))>0.25):
                i = i+'!TOTAL_MIX$1'
                tm = True
            
            if len(set(det_pos_2)&set(common_pos_2))>1 and len(physical_double)>1 and (len(set(det_pos_2)&set(common_pos_2))/(len(set(det_pos_1))+len(set(det_pos_2)))>0.25  or len(set(det_pos_2))/len(set(det_pos_1))>0.25):
                i = i+'!TOTAL_MIX$2'
                tm = True
            if tm==False:
                i= i+'!MIX'
            
        if len(physical_double)>0: #real, 1 or 2 can be accidents
            left_df = oneVSone.T
            left_df = left_df[(left_df['taxon'].isin(left)) | (left_df.index.isin(physical_double)==False)]
            left_df = left_df.T
            # oneVSone = oneVSone.loc[:, ~(oneVSone == lca).any()]#remove LCA in 1VS1 comparison
            score = []
            total_average = np.sum(left_df.loc['weight'])
            for k in set(left_df.loc['taxon'].values):
                temp = left_df.loc[:, (left_df == k).any()]
                temp = temp.loc['weight'].values
                score.append([k,np.sum(temp)/total_average])
            #reiter the physical mixture analysis, keep mixture if sign that one side complete
            temp_i = i
            if any(element[0] not in left and element[0]!=lca for element in score)==False and 'TOTAL_MIX$1' not in i and len(left_df.loc[:, (left_df.isin(left)).any()].columns)>2:
                if '!MIX' in i:
                    temp_i = temp_i.replace('!MIX','')
                if 'TOTAL_MIX$2' in i:
                    temp_i = temp_i.replace('TOTAL_MIX$2','TOTAL_MIX$1!TOTAL_MIX$2')
                else:
                    temp_i = temp_i+'!TOTAL_MIX$1'
            out[temp_i+'!1']=score
            #############RIGHT SIDE##############
            r_df = oneVSone.T
            r_df = r_df[(r_df['taxon'].isin(right)) | (r_df.index.isin(physical_double)==False)]
            r_df = r_df.T
            
            score = []
            total_average = np.sum(r_df.loc['weight'])
            for k in set(r_df.loc['taxon'].values):
                temp = r_df.loc[:, (r_df == k).any()]
                temp = temp.loc['weight'].values
                score.append([k,np.sum(temp)/total_average])
            temp_i = i
            if any(element[0] not in right and element[0]!=lca for element in score)==False and 'TOTAL_MIX$2' not in i and len(r_df.loc[:, (r_df.isin(right)).any()].columns)>2:
                if '!MIX' in i:
                    temp_i = temp_i.replace('!MIX','')
                temp_i = temp_i+'!TOTAL_MIX$2'
            out[temp_i+'!2']=score
            
        else:
            # oneVSone = oneVSone.loc[:, ~(oneVSone == lca).any()]#remove LCA in 1VS1 comparison
            score = []
            total_average = np.sum(oneVSone.loc['weight'])
            for k in set(oneVSone.loc['taxon'].values):
                temp = oneVSone.loc[:, (oneVSone == k).any()]
                temp = temp.loc['weight'].values
                score.append([k,np.sum(temp)/total_average])
            out[i+'!0']=score
    return out

def checking_physical_mix(k,theoretical_trace_back,physical_mix_trace_back):
    removing = False
    if 'theoretical' in k:
        k = k.split('_VS_')
        left = k[0]
        right = k[1]
        if 'theoretical' in left and 'theoretical' in right:
            removing = [left]
        #If theoretical not in both, we have true subsetters, because a theoretical cannot be compared with species it is comprised of
        elif 'theoretical' in left:
            #means that a combo of species equals a species we have, all under the left need to be removed
            removing = [left]
        else:
            #means that a combo of species equals a species we have, all under the left need to be removed
            removing = [right]
    return removing

def path_to_lca(taxonomy,ancestor,s):
    compare = s.split('_VS_')
    c1 = taxonomy[compare[0]]
    twodim = []
    for x in c1:
        if ancestor not in x:
            twodim.append(x[1])
        else:
            break
    twodimreverse = []
    c2 = taxonomy[compare[1]]
    for x in c2:
        if ancestor not in x:
            twodimreverse.append(x[1])
        else:
            break
    return twodim,twodimreverse

def add_to_consensus(total_consensus_df,new_theoretical,total_insilico,physical_mix):
    change_consensus_df = [(key,val) for key,val in total_consensus_df.items()]
    doubled_sequences = {k:[] for k in new_theoretical.keys()}
    for protein_al,df_c in change_consensus_df: 
        insilico = total_insilico[protein_al]

        for k,v in new_theoretical.items():
            temp1 = df_c.loc[v[0]]
            temp2 = df_c.loc[v[1]]
            pm = physical_mix['theoretical_'+str(k)]
            to_add = []
            insilico_match = []
            for x in range(0,len(df_c.columns)):
                if len(insilico[x])>1:
                    insilico_match.append(1)
                else:
                    insilico_match.append(0)
                if temp1[x]==temp2[x]:
                    to_add.append(temp1[x])
                elif temp1[x]=='-':
                    to_add.append(temp2[x])
                elif '-'==temp2[x]:
                    to_add.append(temp1[x])
                elif pm == 1 or pm==0:#left side mix
                    if temp1[x] in insilico[x] and temp1[x]!='X':#first check for presence of left sided residue
                        to_add.append(temp1[x])
                    elif temp2[x] in insilico[x] and temp2[x]!='X':
                        to_add.append(temp2[x])
                    else:
                        to_add.append('X')
                elif pm == 2: #right side mix
                    if temp2[x] in insilico[x] and temp2[x]!='X':#first check for presence of right sided residue
                        to_add.append(temp2[x])
                    elif temp1[x] in insilico[x] and temp1[x]!='X':
                        to_add.append(temp1[x])
                    else:
                        to_add.append('X')
            seq_add = ''.join([el for locs,el in enumerate(to_add) if insilico_match[locs]==1])
            for seq in df_c.index:
                left_over = ''.join([el for locs, el in enumerate(df_c.loc[seq].values) if insilico_match[locs]==1])
                if left_over ==seq_add:
                    doubled_sequences[k].extend([seq])
            df_c.loc['theoretical_'+str(k)]=to_add
        total_consensus_df[protein_al]=df_c
    ds = {}
    for k,v in doubled_sequences.items():
        v = [el for el in doubled_sequences if v.count(el)>=len(change_consensus_df)]#all proteins considered need to be the same
        if len(v)>0:
            ds[k]=['theoretical_'+str(k)]#we do not want the newest identical match
    return total_consensus_df, ds

def find_most_likely_combo(cs,td,ac,taxonomy,lca,theoretical_trace_back):
    combo_likely = {}
    data_temp = []
    for k in cs:
        temp = []
        td_k = td[k]
        td_k = [el for el in td_k if sorted([el,k]) not in ac]
        for k2 in cs:
            trace1 = find_trace(k,theoretical_trace_back)
            trace2 = find_trace(k2,theoretical_trace_back)
            if k==k2 or sorted([k,k2]) in ac or k2 not in td_k:
                temp.append(1000) 
            elif len(set(trace1)&set(trace2))>0:
                temp.append(1000) 
            else:
                temp.append(td_k.index(k2))
        data_temp.append(temp)
    data_df = pd.DataFrame(data_temp,columns=cs, index=cs)
    summed = []
    for k in cs:
        for k2 in cs:
            if any(element[0] == sorted([k,k2]) for element in summed):
                continue
            summed.append([sorted([k,k2]),data_df[k][k2]+data_df[k2][k]])
    summed = sorted(summed, key=lambda x:x[1])
    done = []
    for element in summed:
        if any(el in element[0] for el in done) or find_LCA(taxonomy,element[0])==lca:
            combo_likely[tuple(element[0])]=False
        else:
            combo_likely[tuple(element[0])]=True
            done = done + element[0]
    ranked_input = []
    others = []
    for k,v in combo_likely.items():
        if k[0]!=k[1] and v==True:
            ranked_input = ranked_input + list(k)
        elif v==True:
            others.append(k[0])
    return ranked_input + others
    

def find_combinations_of_species(considered_species,theoretical_trace_back,already_compared,taxonomy,final_output_genetic_mixtures):
    #make closest relatives first
    #higher combos not possible if subsetter still in the mix
    #maximize amount of species that can be combined, now 4 species toghether
    out_species = []
    temp_combined = []
    grouped_species=[]
    traces = {}

    considered_species = [el for el in considered_species if 'theoretical' in el]+[el for el in considered_species if 'theoretical' not in el]
    td_considered_species = relatedness_finder(considered_species,taxonomy)
    LCA = find_LCA(taxonomy, considered_species)
    done = []
    most_likely_combo = find_most_likely_combo(considered_species,td_considered_species,already_compared,taxonomy,LCA,theoretical_trace_back)
    #group very close related species first
    for x in most_likely_combo + [el for el in considered_species if el not in most_likely_combo]:
        traces[x]=find_trace(x,theoretical_trace_back)
    for x in most_likely_combo + [el for el in considered_species if el not in most_likely_combo]:
        if any(x in element for element in out_species):
            continue
        trace_x = find_trace(x,theoretical_trace_back)
        traces[x]=trace_x
        for x2 in td_considered_species[x]:
            if x==x2 or any(x2 in element for element in out_species):#Do not compare sample to itself
                continue
            if x in grouped_species or x2 in grouped_species or find_LCA(taxonomy,[x,x2])==LCA:
                continue
            trace_x2 = find_trace(x2,theoretical_trace_back)
            if sorted([x,x2]) in already_compared or sorted(trace_x)==sorted(trace_x2):#already compared #if x in taxonomy and x2 in taxonomy
                continue
            #When the combination has already been compared before with more animals included + same for combinations that will be checked
            if any(len(set(find_trace(el[0],theoretical_trace_back)+find_trace(el[1],theoretical_trace_back))&set(trace_x+trace_x2))==len(set(trace_x+trace_x2)) and len(set(find_trace(el[0],theoretical_trace_back)+find_trace(el[1],theoretical_trace_back)))>len(set(trace_x+trace_x2)) for el in already_compared+out_species):
                continue
            #no need in comparing 2 same species because result =0
            if sorted([x,x2]) not in out_species:
                if sorted(trace_x+trace_x2) in temp_combined or len(set(trace_x)&set(trace_x2))>0 or len(set(trace_x + trace_x2))>8:
                    continue #either next round or not needed
                if find_LCA(taxonomy,[x,x2])==LCA and len([find_LCA(taxonomy,[x,element])!=LCA for ele in considered_species for element in td_considered_species[ele] if ele !=x and ele not in done and sorted([x,element]) not in already_compared])>0:
                    #if a better match is possible lower than LCA level, this needs to be taken
                    continue
                temp_combined = temp_combined + [sorted(trace_x+trace_x2)]
                out_species.append(sorted([x,x2]))
                done = done + [x,x2]
                grouped_species.append(x)
                grouped_species.append(x2)
                break#1 comaprison per iteration, otherwise too many for quick search
    #group further species, LCA != total LCA
    already_compared_but_mixed = [el for el in considered_species if any(el in num for num in out_species)==False]
    if len(already_compared_but_mixed)==len(considered_species):
        LCA=[]
    done = []
    if (len(out_species)>0 or LCA==[]) and len(already_compared_but_mixed)>0:
        for x in already_compared_but_mixed:
            if any(x in element for element in out_species):
                continue
            closest = td_considered_species[x]  
            if x in grouped_species:#allocated with another one that was lca distant
                continue
            for x2 in closest:
                if any(x2 in element for element in out_species):
                    continue
                if sorted([x,x2]) in out_species or sorted([x,x2]) in already_compared or len(set(traces[x])&set(traces[x2]))>0:
                    continue
                if x2 in grouped_species and x2 in already_compared_but_mixed:#this means that it will be included already, otherwise to much senseless 1VS1 comaprisons
                    continue
                if find_LCA(taxonomy,[x,x2])==LCA and len([find_LCA(taxonomy,[x,element])!=LCA for ele in already_compared_but_mixed for element in td_considered_species[ele] if ele !=x and ele not in done and sorted([x,element]) not in already_compared])>0:
                    #this means that there is another option better than this one. 
                    continue
                elif find_LCA(taxonomy,[x,x2])!=LCA:
                    out_species.append(sorted([x,x2]))
                    done = done + [x,x2]
                    grouped_species.append(x)
                else:
                    if find_LCA(taxonomy,[x2,td_considered_species[x2][0]])==LCA:#if both LCA distance to closest than they can be combined
                        out_species.append(sorted([x,x2]))
                        done = done + [x,x2]
                        grouped_species.append(x)
                break
    already_compared_but_mixed = [el for el in considered_species if any(el in num for num in out_species)==False]
    r = []
    if len(already_compared_but_mixed)==len(considered_species):#final iteration is passed, we stop here
        for x in already_compared_but_mixed:
            r.append(x)
            if 'theoretical' in x:
                x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
            else:
                x = (x,x,0)
            final_output_genetic_mixtures[x]=1
    elif len(already_compared_but_mixed)>0:#keep them for the next round
        for x in already_compared_but_mixed:
            out_species.append([x,x])
            grouped_species.append(x)
    
    return out_species,r,final_output_genetic_mixtures

def further_iteration(cpu_count,total_consensus_df,total_insilico,out_species,species_in,taxonomy,df):
    with ProcessPoolExecutor(max_workers=cpu_count) as executor3:
        mixtures_all = {}
        pep_locs_all = {}
        with tqdm(list(total_consensus_df), desc="Starting on comparison...") as pbar:
            for protein_al in pbar:
                pbar.set_description(f"Performing 1-VS-1 comparison on {protein_al}")
                insilico= total_insilico[protein_al]
                df_c = total_consensus_df[protein_al]
        
                if len(out_species)!=0:#make the consensus sequence here already so we don't have to do it later after convergence
                    try:
                        results = []
                        futures = [executor3.submit(look_for_complementarity,compare,df_c,df,species_in,protein_al,taxonomy,insilico) \
                                for compare in out_species]
                        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing", leave=False):
                            results.append(future.result())
                        mixtures = {key:val for key, val in results}
                        pep_locs_all[protein_al] = np.array([el+'_'+str(num) for num,el in enumerate(df_c.columns) if 'X' in el])
                        mixtures_all[protein_al]=mixtures
                    except Exception as e:
                        pbar.write(f"Error processing {protein_al}: {e}")
                    finally:
                        pbar.update(1)
    return mixtures_all,pep_locs_all

def plot_new_tree(dist_matrix,name_to_coverage,file_extinct,path, sample_path,taxonomy):
    print('Plotting taxonomic tree')
    fig = ff.create_dendrogram(dist_matrix, 
                               labels = [el+'_BC-score= '+name_to_coverage[el] if el in name_to_coverage else el for el in dist_matrix.index],
                               )
    fig.update_layout(autosize=True,title=file_extinct)
    file_extinct = file_extinct.replace('Taxonomic_results_','')
    
    # labs = fig['layout']['xaxis']['ticktext']#get labels
    # x_locs = fig['layout']['xaxis']['tickvals']#get x-axis location labels
    # data = fig['data']
    # lca = find_LCA(taxonomy,[el for el in labs if 'BC-score' not in el])
    # clusters = {}
    # species_line = {}
    # y_dist = []
    # for loc,species in enumerate(labs):
    #     if 'BC-score' in species or 'Theoretical' in species or 'Database' in species:#this is what we are trying to locate
    #         continue
    #     aver = x_locs[loc]#get X loc of tree
    #     species_line[species]=[aver]
    #     added = True
    #     former_y = 0
    #     while added==True:#lca_now != lca:
    #         added = False
    #         for q in data:
    #             if any(aver-0.001<=el<=aver+0.001 for el in q['x']) and max(q['y'])>=former_y:#new_node
    #                 former_y = max(q['y'])
    #                 aver = np.average(q['x'])
    #                 y_dist = y_dist + list(q['y'])
    #                 added = True
    #                 if aver in species_line[species]:
    #                     added=False
    #                 species_line[species]=species_line[species]+[aver]
    #                 break
    # added_already = []
    # x_dist = (x_locs[0]+x_locs[1])/4
    # y_dist = max(y_dist)/(x_dist*4)
    # zeros = []
    # for loc,species in enumerate(labs):
    #     if 'BC-score' in species or 'Theoretical' in species or 'Database' in species:#this is what we are trying to locate
    #         continue
        
    #     taxa = taxonomy[species]#get lineage
    #     aver = x_locs[loc]#get X loc of tree
    #     added_to_lin = [species]#to recover node LCA
    #     added = True
    #     former_lca = species
    #     while added==True:#lca_now != lca:
    #         added = False
    #         for q in data:
    #             if aver in q['x']:#new_node
    #                 aver = np.average(q['x'])
    #                 added = True
    #                 x0=aver-x_dist
    #                 x1=aver+x_dist
                    
    #                 y0=max(q['y'])-y_dist/5
    #                 y1=max(q['y'])+y_dist/5
    #                 if min(y0,y1)<=0:
    #                     y0 =0
    #                 include = []
    #                 for sp,av in species_line.items():
    #                     if aver in av:
    #                         include.append(sp)
    #                 if len(include)>0:
    #                     lca_now= find_LCA(taxonomy,include+[species])
    #                 else:
    #                     lca_now = lca
    #                 level = [el[0] for el in taxa if el[1]==lca_now][0]
    #                 if lca_now == former_lca:
    #                     continue
    #                 former_lca = lca_now
    #                     # lca_now = [taxa[num+1][1] for num,el in enumerate(taxa) if el[0]==level][0]
    #                 if max(q['y'])==0 and lca_now in zeros:
    #                     continue#no point in having these genusses plotted again
    #                 elif  max(q['y'])==0:
    #                     zeros.append(lca_now)#we plot 1 genus collection
    #                 if (x0,x1,y0,y1) in added_already:
    #                     continue
    #                 added_already.append((x0,x1,y0,y1))
    #                 adding = [dict(type="rect",
    #                         xref="x", yref="y",
    #                         x0=x0, y0=y0,
    #                         x1=x1, y1=y1,
    #                         opacity=0.5,
    #                         fillcolor="cornflowerblue",
    #                         line_color="black",
    #                         layer="above",
    #                         line_width=0,
    #                         label = dict(text=lca_now),)]
                    
    #                 if level in clusters:
    #                     if adding not in clusters[level]:
    #                         clusters[level]=clusters[level]+adding
    #                 else:
    #                     clusters[level]=adding
    #                 break               

    
    # menus = [dict(label="None",
    #      method="relayout",
    #      args=["shapes", []]),
    #          ]
    # all_args = []
    # for cl,m in clusters.items():
    #     add = dict(label=cl,method = "relayout",args=['shapes',m])
    #     all_args= all_args+m
    #     menus.append(add)
    # menus.append(dict(label='All',method = "relayout",args=['shapes',all_args]))
    # fig.update_layout(
    # updatemenus=[
    #     dict(buttons=list(menus),
    #         )],
    # paper_bgcolor='rgba(0,0,0,0)',
    # plot_bgcolor='rgba(0,0,0,0)')
    # fig.update_yaxes(showticklabels=False,ticks=None, title=None)
    # fig.update_xaxes(ticklabelposition='outside bottom',ticks='outside')
    
    
    name_file=file_extinct+'.html'
    try:
        fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'taxonomic_output' /name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'taxonomic_output' / 'measured_tree_output.html') 
    
    return

def transform_taxons(transform_taxa,taxonomy,sp_in,golca,lca_all):#{},taxonomy,species_in,golca,lca_all
    for n in sp_in:
        t = taxonomy[n]
        adjusted_taxa = {}
        loc = 0
        while any('genus' in el for el in t) and 'genus' not in t[loc]:
            loc += 1
        if loc==0:#genus not in lineage or first in lineage
            adjusted_taxa[loc] = golca[0]
        else:
            adjusted_taxa[0]=0
            adding = golca[0]/(loc)
            for i in range(1,loc+1):
                adjusted_taxa[i] = max(adjusted_taxa.values())+adding
        remember_loc = loc
        while any('order' in el for el in t) and 'order' not in t[loc]:
            loc += 1
        if loc==0:#genus not in lineage or first in lineage
            adjusted_taxa[loc] = golca[0]
        else:
            if loc-remember_loc>0:
                adding = (golca[1]-max(adjusted_taxa.values()))/(loc-remember_loc)
            else:
                adding = (golca[1]-max(adjusted_taxa.values()))
            for i in range(remember_loc+1,loc+1):
                adjusted_taxa[i] = max(adjusted_taxa.values())+adding
                
        remember_loc = loc
        if any(lca_all in el for el in t[:loc])==False:
            while any(lca_all in el for el in t) and lca_all not in t[loc]:
                loc += 1
                if loc>len(t):
                    break
            if loc==0:#genus not in lineage or first in lineage
                adjusted_taxa[loc] = golca[0]
            else:
                if loc-remember_loc>0:
                    adding = (golca[2]-max(adjusted_taxa.values()))/(loc-remember_loc)
                else:
                    adding = (golca[2]-max(adjusted_taxa.values()))
                for i in range(remember_loc+1,loc+1):
                    adjusted_taxa[i] = max(adjusted_taxa.values())+adding
        
        transform_taxa[n]=adjusted_taxa
    return transform_taxa

def find_species_intervals(in_names,tree_matrix, transform_taxa,dist_matrix):#added
    species_intervals = {}
    for i in in_names:
        intervals = {}
        if len(sorted(list(set(tree_matrix[i].values)))[1:])>0:
            inters = sorted(list(set(tree_matrix[i].values)))[1:]
        else:
            inters = [sorted(list(set(tree_matrix[i].values)))[0]]
        for x in inters:#do not take 0 as first, because it is start value
            #value linked to species that are linked to next taxon
            if len(intervals)==0:
                previous_value = 0
            else:
                previous_value = max(intervals[max(intervals.keys())])
            linked = tree_matrix[i][tree_matrix[i]==x].index
            next_value = [el for el in dist_matrix[i][linked].values if el>previous_value]#others will be further away in the tree but taxonomically at the same distance
            if len(next_value)>0:
                if x!=max(sorted(list(set(tree_matrix[i].values)))[1:]):
                    intervals[x]=[previous_value,np.average(np.array(next_value))]
                else:
                    intervals[x]=[previous_value,max(next_value)]
            else:
                intervals[x]=[previous_value,previous_value]
        intervals[max(intervals.keys())+1]=[max(intervals[max(intervals.keys())]),1]#in case something at the order level is even further away from measured
        species_intervals[i]=intervals
    return species_intervals

def main_side_lineages(ranked, found,y_labs,species_to_taxon,taxonomy,lca_loc,golca,o_n):
    temp_data = {'line_x':[], 'line_y':[],'x':[],'x_adj':[],'y':[]}
    main_lineage = []
    side_lineages = {}
    lower_scored = []
    sides = []
    for loc,el in enumerate([el for el in ranked[found].values][:lca_loc]):
        taxon = y_labs[loc]
        tax_lin = [element[1] for element in taxonomy[taxon]]
        if ranked.index[loc] != found:
            same_scoring = [y_labs[locs] for locs,element in enumerate([el for el in ranked[found].values][:lca_loc]) if element==el]#taxa same distance
            same_scoring_lca = find_LCA(taxonomy,same_scoring)#LCA of same distance taxa
            temp_data['x'].extend([el])#distance
            if species_to_taxon[ranked.index[loc]]==0 or el==0:#species match, perfect taxon match
                temp_data['x_adj'].extend([el])#no arrow needed
                temp_data['line_x'].extend([el,el,None])#no arrow needed
                #addition to the main or side lineages
                if len(main_lineage)==0 and any(element in y_labs for element in tax_lin):#ADJUST FOR MUTLIPLE EQUAL SCORING OVIS ammon-aries-...!!!!!!!!!!!!!!!!!!!!
                    main_lineage = tax_lin
                elif any(element in main_lineage for element in tax_lin):
                    side_lineages[taxon]=[element for element in tax_lin if element not in main_lineage]+[[element for element in tax_lin if element in main_lineage][0]]
                    sides= sides+[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
            elif (len(main_lineage)==0 or taxon in main_lineage or any(element in main_lineage for element in same_scoring+[same_scoring_lca])) and (any(element in y_labs for element in tax_lin) or taxon in main_lineage):
                #check if the taxon is part of the main or side lineage
                if any(element[1] in lower_scored for element in taxonomy[taxon])==False or taxon in main_lineage:#part of the main branch
                    adj = el-species_to_taxon[ranked.index[loc]]
                    if adj<0:#error because of database deficiency at this lowest level, sub-species can be possible too, which will be below 0 too.
                        adj=abs(adj)/2
                    temp_data['x_adj'].extend([adj])
                    temp_data['line_x'].extend([el,adj,None])
                    main_lineage = main_lineage + [element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]#first addition
                else:
                    adj = el+abs(el-species_to_taxon[ranked.index[loc]])
                    temp_data['x_adj'].extend([adj])
                    temp_data['line_x'].extend([el,adj,None])
                    side_lineages[taxon]=[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
                    sides= sides+[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
            else:#side lineages
                if any(element[1] in main_lineage for element in taxonomy[taxon]):
                    side_lineages[taxon]=[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
                    sides= sides+[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
                    adj = el+abs(el-species_to_taxon[ranked.index[loc]])
                    temp_data['x_adj'].extend([adj])
                    temp_data['line_x'].extend([el,adj,None])
                else: #if taxon in sides:
                    adj = el+abs(el-species_to_taxon[ranked.index[loc]])
                    temp_data['x_adj'].extend([adj])
                    temp_data['line_x'].extend([el,adj,None])
                    side_lineages[taxon]=[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]
                    sides= sides+[element[1] for element in taxonomy[taxon] if element[1] not in main_lineage]+[[element[1] for element in taxonomy[taxon] if element[1] in main_lineage][0]]   
            if taxon!=o_n:
                lower_scored.append(taxon)
            temp_data['line_y'].extend([taxon,y_labs[loc],None])
            temp_data['y'].extend([taxon])
    #combine_side_lineages
    colour_list = ['olive','teal','navy','purple','orange','yellowgreen','goldenrod','gold','seagreen']*10
    colours = {k:colour_list[loc] for loc,k in enumerate(main_lineage) if k in temp_data['y']}
    side_lin = {k:[] for k in main_lineage if k in temp_data['y']}
    side_lin['Unlinked']=[]
    for k,v in side_lineages.items():
        if v[-1] in colours:  
            colours[k]=colours[v[-1]]
            side_lin[v[-1]].extend([k])
        else:
            colours[k]='brown'
            side_lin['Unlinked'].extend([k])
    colours_out = {}
    for k,v in colours.items():
        if k not in main_lineage:
            colours_out[k]=v
        else:
            colours_out[k]='black' 
    colours = [colours_out[element] for element in temp_data['y']]
            
    return temp_data, colours,main_lineage,side_lin

def find_most_likely(taxonomy,species,dots,transform_taxa,order_l):
    points = []
    points_species = []
    link = [(species[el],dots[el]) for el in range(0,len(species))]
    link = sorted(link, key=lambda x:x[1])
    link_dict = {k:v for k,v in link}
    done = []
    most_relevant_species = [0]
    assigned = []
    range_assigned = [0]
    for el in link:
        lin = [element[1] for element in taxonomy[el[0]]]
        if el[1]>order_l:#not relevant
            continue
        if any(ele[1]-el[1]<2 for ele in link if ele[0] in lin and ele[0]!=el[0])==False and el[0].count(' ')>0:#non-sense outcome, taxa should be in a close ladder
            continue
        if el[0].count(' ')>0 and len(most_relevant_species)==1:#adjust most relevant species outcome, other species are thus less relevant and more distant 
            most_relevant_species.append(el[1])
        if transform_taxa[el[0]][0]>el[1] and el[1]>1 and (el[1]!=min(dots) and min(dots)>3):#not fitting, and not the minimal
            continue
        if el[0].count(' ')>0 and (el[1]>3 or max(most_relevant_species)<el[1]):#species matches in database need to be distance <2 to be somewhat relevant
            continue
        track = [el[1]]
        for q in lin:
            if q not in link_dict:#not in options so cannot be checked anyways
                continue
            if (q==el[0] or q in done) and link_dict[q]!=min(dots):#taxonomy not in result
                continue
            elif link_dict[q]==min(dots) and len([ele for ele in lin if ele in assigned or ele in done])<=1:#Same distance as the minimum so as likely
                track.append(link_dict[q])
            elif link_dict[q]<=max(range_assigned) and len([ele for ele in lin if ele in assigned or ele in done])<=1 and el[0].count(' ')==0: #higher taxon is lower than this one
                track.append(link_dict[q])
            elif link_dict[q]>=track[-1] and (len(assigned)<3 or (el[0].count(' ')>0 and len([ele for ele in lin if ele in assigned or ele in done])<=1)):#within lineage of most likely, and within top3 star
                track.append(link_dict[q])
            else:
                break
        if len(track)>1 or track[0]<=max(range_assigned):
            # print(el[0])
            points.append(el[1])
            points_species.append(el[0])
            if (len(assigned) <3 or el[1] in range_assigned) and el[0].count(' ')==0:#higher taxa top 3
                assigned.append(el[0])
                range_assigned.append(el[1])
        done.append(el[0])
    return points,points_species

def find_distance_species_taxon(transform_taxa,names,additional,taxonomy,golca,lca_all):
    sp_to_tax = {}
    for i in names:
        if i in additional and len(additional[i])>0:
            candidate = additional[i][0]
            c_og = candidate
            lin = taxonomy[candidate]
            while 'species' not in lin[0][0]:
                try:#in case weird things happen
                    children = find_child_nodes(taxonomy,candidate)
                    children = [el for el in children if 'sp.' not in el and el.isdigit()==False]
                    children = [el for el in children if c_og in [ta[1] for ta in taxonomy[el]]]
                except:
                    break
                if len(children)==0:#means that some stupid fock gave the same name to something not mammalian
                    lin = [v for k,v in taxonomy.items() if any(c_og in el for el in v)][0]
                    break
                candidate = children[0]
                lin = taxonomy[candidate]
            taxonomy[candidate]=lin
            transform_taxa = transform_taxons(transform_taxa,taxonomy,[candidate],golca,lca_all)
            lin = [el[1] for el in lin]
            loc = lin.index(i)
        else:
            candidate = i
            try:
                lin = [el[1] for el in taxonomy[candidate]]
            except: 
                children = find_child_nodes(taxonomy,candidate)
                
                children = [el for el in children if 'sp.' not in el and el.isdigit()==False]
                children = [el for el in children if candidate in [ta[1] for ta in taxonomy[el]]]
                candidate = children[0]
                taxonomy[candidate]=taxonomy[candidate]
                transform_taxa = transform_taxons(transform_taxa,taxonomy,[candidate],golca,lca_all)
                lin = [el[1] for el in taxonomy[candidate]]
            loc = lin.index(candidate)
        sp_to_tax[i]=transform_taxa[candidate][loc]
    return sp_to_tax, transform_taxa

def recover_sequences(dist_matrix,final_combined_multiple_sequence_alignment,sequences,genetic_mixtures,species_in,TF,subs,TF_all):
    seqs = {}
    split_data = set([el.split('_')[-1] for el in final_combined_multiple_sequence_alignment.columns])
    species_in = [el for el in species_in if ' ' in el]
    for sp in genetic_mixtures:
        if sp.replace('Database_match_','') not in final_combined_multiple_sequence_alignment.index:
            continue #will be there for the other order group
        temp = [el for el in dist_matrix.columns if (el in species_in and el in final_combined_multiple_sequence_alignment.index) or sp==el]
        temp = dist_matrix[temp]
        temp = temp.loc[temp.columns]
        closest = sorted([(el,list(temp.columns)[loc]) for loc,el in enumerate(list(temp[sp].values))], key=lambda x:x[0])
        closest = [el for el in closest if el[1] in species_in and el[1] in dist_matrix.columns][0][1]
        for number in split_data:
            temp = [ind for ind,el in enumerate(final_combined_multiple_sequence_alignment.columns) if '_'+number in el]
            found = [TF[ind] for ind,el in enumerate(final_combined_multiple_sequence_alignment.columns) if '_'+number in el]
            found_residue = [TF_all[ind] for ind,el in enumerate(final_combined_multiple_sequence_alignment.columns) if '_'+number in el]
            temp = final_combined_multiple_sequence_alignment.iloc[:,temp]
            th_sequence = ''.join(temp.loc[sp.replace('Database_match_','')])
            cl_sequence = ''.join(temp.loc[closest])
            cl_all = []
            for k,v in sequences.items():
                if closest not in v:
                    continue
                matches = [el for el in cl_sequence.split('K') if el in str(k)]
                cl_all.append((v,len(matches)))
            matches = sorted(cl_all,key=lambda x:x[1])[-1][0]
            correct = []
            for loc,i in enumerate(th_sequence):
                if i in found_residue[loc]:
                    correct.append(1)
                else:
                    correct.append(0)
            found = [el if correct[loc]==1 else 0 for loc,el in enumerate(found)]
            seqs[sp+'_'+number]=[th_sequence,cl_sequence,matches,found]
    return seqs

def thread_align_mix(i,seq,m,calculator):
    if i==seq:
        return 0
    i = Seq(str(i).replace('&',''))
    seq = Seq(str(seq).replace('&',''))
    try:
        aligner = Align.PairwiseAligner()
        align = aligner.align(i,seq)
        align = next(align)
        # align = list(align)[0]
        a=SeqRecord(align[0].replace('-','Z'),id='a')
        b=SeqRecord(align[1].replace('-','B'),id='b')
        align = MultipleSeqAlignment([a, b])
        dm = calculator.get_distance(align)
        return dm.matrix[1][0]
    except:
        return 10 

def go_to_species_mix(c,species_list,taxonomy):
    for x in c:
        if x.lower()=='environmental samples' or 'sp.' in x:
            continue
        t_lin = taxonomy[x]
        if 'species' in t_lin[0][0]:
            species_list.append(x)#we are not interested in subspecies, because highly variable 
            try:
                gtsp = [el for el in find_child_nodes[x] if 'sp.' not in el and 'unclassified' not in el and el.isdigit()==False]
                species_list = species_list +gtsp
            except:
                continue
        else:
            try:
                gtsp = [el for el in find_child_nodes[x] if 'sp.' not in el and 'unclassified' not in el and el.isdigit()==False]
                gts = []
                for g in gtsp:
                    g_no_numbers = ''.join([num for num in g if num.isdigit()==False])
                    if g_no_numbers==g:
                        gts.append(g)
                species_list = species_list + go_to_species_mix(gts,[],taxonomy)
            except:
                continue
    return species_list

def background_check_mix(c,ip,taxonomy):
    out = {}
    for q in c:
        m = go_to_species_mix([q],[],taxonomy)
        if len(set(m)&set(ip))==0:
            out[q] = m
    return out


def plot_taxonomy(file_extinct,theoretical_trace_back,species_in,sequences,others,taxonomy,combined_file_output,all_starting_peptides,df,path,sample_path,taxonomy_missing):
    final_output_genetic_mixtures = {}
    for k,v in combined_file_output.items():
        final_output_genetic_mixtures = final_output_genetic_mixtures | v[0]
    genetic_mixtures = [el for el,val in final_output_genetic_mixtures.items() if val>0]    
    gm = []
    for num,el in theoretical_trace_back.items():
        if any(num==element[-1] for element in genetic_mixtures)==True:#find theoreticals
            gm.append('theoretical_'+str(num))
        elif any(element.count(el[0])>1 and element.count(el[1])>1 for element in genetic_mixtures)==True:#find species
            gm.append(num)
    genetic_mixtures=gm
    genetic_mixtures = list(set(genetic_mixtures))
    #locate the genetic mixture in the taxonomic tree
    fcmsa_per_group = {}
    final_combined_multiple_sequence_alignment_per_group = {}
    species_to_peptides = {el:[] for el in genetic_mixtures}
    all_species = []
    TF_per_group = {}
    total_insilico_all = {}
    for group,v in combined_file_output.items():
        total_consensus_df = v[1]
        total_insilico = v[2]
        final_combined_multiple_sequence_alignment = pd.DataFrame()
        fcmsa = pd.DataFrame()
        TF = []
        TF_all = []
        protein_nr = 0
        for protein_al,df_c in total_consensus_df.items():
            protein_nr += 1
            insilico = total_insilico[protein_al]
            binary = [0 if len(el)==1 else 1 for el in insilico]
            TF = TF+[0 if len(el)==1 else 1 for el in insilico]
            TF_all = TF_all + [el for el in insilico]
            df_add_f = df_c.iloc[:, [num for num,el in enumerate(binary)]]
            df_add_f.columns = [el+'_'+str(protein_nr) for el in df_add_f.columns]
            fcmsa = pd.concat([fcmsa, df_add_f], axis=1)
            
            df_add = df_c.iloc[:, [num for num,el in enumerate(binary) if el!=0]]
            df_add.columns = [el+'_'+str(protein_nr) for el in df_add.columns]
            final_combined_multiple_sequence_alignment = pd.concat([final_combined_multiple_sequence_alignment, df_add], axis=1)
            for i in species_to_peptides.keys():
                i_ori = i
                i = i.replace('Database_match_','')
                if i not in final_combined_multiple_sequence_alignment.index:
                    continue
                test_species = ''.join(final_combined_multiple_sequence_alignment.loc[i]).replace('-','')
                species_to_peptides[i_ori] = list(set(species_to_peptides[i_ori]+[el for el in all_starting_peptides if el in test_species]))#redistribution of the peptide contents
        for q in genetic_mixtures:
            if 'Database_match' in q:
                if q.replace('Database_match_','') not in final_combined_multiple_sequence_alignment.index:
                    continue
                final_combined_multiple_sequence_alignment.loc[q]=final_combined_multiple_sequence_alignment.loc[q.replace('Database_match_','')]
                fcmsa.loc[q]=fcmsa.loc[q.replace('Database_match_','')]
        fcmsa_per_group[group] = fcmsa
        TF_per_group[group]=TF
        total_insilico_all[group]=TF_all
        final_combined_multiple_sequence_alignment_per_group[group] = final_combined_multiple_sequence_alignment
        all_species = list(set(all_species + [el for el in final_combined_multiple_sequence_alignment.index if el in species_in]))
    #remove subsetters here
    name_to_coverage = {}
    subsetters = []
    compared = []
    all_peps = set()
   
    for k,v in species_to_peptides.items():
        all_peps = all_peps|set(v)
        compared.append(k)
        jv = '!'.join(v)
        for k2,v2 in species_to_peptides.items():
            jv2 = '!'.join(v2)
            if k2 in compared or k2 == k or k2 in subsetters or k in subsetters:#already done or already discarded
                continue
            if any(el not in jv for el in v2) and any(el not in jv2 for el in v):#both species have unique peptides
                continue#so we keep it
            if any(el not in jv2 for el in v)==False and any(el not in jv for el in v2)==False:#species 1 is not a subset
                if set(find_trace(k2,theoretical_trace_back)).issubset(find_trace(k,theoretical_trace_back)):
                    print('{} is in {} but more complementarity was found'.format(k2,k))
                    continue
                print('{} is in {}'.format(k2,k))
                subsetters.append(k2)
    for k in others:
        species_to_peptides['Database_match_'+k]=list(set(df['Peptides'][df['Species']==k]))
    #We only focus for the final Bray-Curtis score on the peptides that were considered. Others will drag the score down which returns a more ambiguise result.
    all_starting_peptides = []  
    for v in species_to_peptides.values():
        all_starting_peptides = all_starting_peptides + v
    all_starting_peptides = list(set(all_starting_peptides))
    array1 = [1]*len(all_starting_peptides)#starting_peptides
    seq_concat = '_'.join(list(set(all_starting_peptides)))
    weights_bc = [1/seq_concat.count(val) for val in all_starting_peptides]
    remove_order = True
    delete_order = []
    all_orders = [el for el in list(combined_file_output.keys())+others if el not in delete_order]
    delete_order_species = []
    while remove_order==True and len(delete_order)<len(all_orders):
        ord_to_sp = {g:[] for g in all_orders}
        order_uniqueness = []
        names_groups = []
        delete_species = []
        order_score = []
        del_ord = []
        for group in [el for el in all_orders if el not in delete_order]:
            order_peps = []
            other_peps = []
            species_gr = []
            sp_to_sp = {}
            for species,peps in species_to_peptides.items():
                if species in delete_order_species:
                    continue
                sp_to_sp[species.replace('Database_match_','')]=species
                species = species.replace('Database_match_','')
                if any(group in el for el in taxonomy[species]):
                    order_peps = order_peps+peps
                    species_gr.append(species)
                    ord_to_sp[group].extend([species])
                else:
                    other_peps = other_peps+peps
            diff = set(order_peps)^set(other_peps)
            if len(set(order_peps)&diff)<=1:
                delete_species = delete_species + species_gr
                del_ord.append(group)
            array2 = [1 if el in order_peps else 0 for el in all_starting_peptides]
            score = braycurtis(array1, array2,w=weights_bc)
            order_score.append(1-score)
            order_uniqueness.append(len(set(order_peps)&diff))
            names_groups.append(group)
        if len(delete_species)==0:
            remove_order=False
        else:
            if 0 in order_uniqueness:
                delete_order.append(sorted([[order_score[loc],names_groups[loc]] for loc,el in enumerate(names_groups) if el in del_ord and el not in delete_order and order_uniqueness[loc]==0])[0][1])
            else:
                delete_order.append(sorted([[order_score[loc],names_groups[loc]] for loc,el in enumerate(names_groups) if el in del_ord and el not in delete_order])[0][1])
            species_discard = [sp_to_sp[sps] for sps in ord_to_sp[delete_order[-1]]]
            print('Removing {} containing: {}'.format(delete_order[-1],species_discard))
            delete_order_species = delete_order_species + species_discard
            
    print('Plotting the uniqness and order BC scores')
    arr_sorted=np.argsort(np.array(order_uniqueness))[::-1],
    names_groups = np.array(names_groups)[arr_sorted]
    order_uniqueness = np.array(order_uniqueness)[arr_sorted]
    order_score = np.array(order_score)[arr_sorted]
    
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Bar(x=names_groups, y=order_uniqueness,
                    marker_color='skyblue',
                    name='Order uniqueness'), secondary_y=False,)
    fig.add_trace(
        go.Scatter(x=names_groups, y=order_score, name="Order BC-score"),
        secondary_y=True,
    )
    fig.update_layout(title=dict(text="Unique peptides at order level {}".format(file_extinct)),
                      xaxis=dict(title=dict(text="Order",font=dict(size=16)),),)
    fig.update_yaxes(title_text="Order level Bc-score", secondary_y=True)
    
    fig.update_yaxes(title_text="Unique peptide count", secondary_y=False)
    
    name_file='Order_uniqueness_'+file_extinct+'.html'
    try:
        fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'taxonomic_output' /name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'taxonomic_output' / 'Before_tree_analysis.html') 
    print('Deleting {} because there is no uniqueness at the order level'.format(delete_order))
    delete_species = delete_order_species
    genetic_mixtures = [el for el in genetic_mixtures if el not in delete_species]+['Database_match_' + el for el in others if any(el in element for element in delete_species)==False]
    species_to_peptides = {k:v for k,v in species_to_peptides.items() if k.replace('Database_match_','') not in delete_species}
    genetic_mixtures = list(set(genetic_mixtures))
    #find combinatorial subsetters
    combo_all = set()
    for k,el in species_to_peptides.items():
        if k not in ''.join(others):
            combo_all = combo_all|set(el)
    sub_combos = []
    for k in [el for el in others if el.replace('Database_match_','') not in delete_species and el.replace('Database_match_','') not in others]:#v=peptides,k=species
        v= species_to_peptides[k.replace('Database_match_','')]
        if len([el for el in set(v) if el not in combo_all])<=2:
            print('{} is a subset of assigned theorethicals'.format(k))
            print([el for el in set(v) if el not in combo_all])
            sub_combos.append(k)
    
    array1 = [1]*len(all_starting_peptides)#starting_peptides
    seq_concat = '_'.join(all_starting_peptides)
    weights_bc = [1/seq_concat.count(val) for val in all_starting_peptides]
    for k,v in species_to_peptides.items():
        k = k.replace('Database_match_','')
        og_k = k
        if k in subsetters+sub_combos:#no need to do it for the discarded ones
            continue
        if k in others or 'theoretical' not in og_k:
            k='Database_match_'+k
        array2 = [1 if el in v else 0 for el in all_starting_peptides]
        score = braycurtis(array1, array2,w=weights_bc)
        name_to_coverage[k]="{:.2f}".format(1-score)#columns need to be the same

    
    print('Deleting subsetters {}'.format(subsetters))
    genetic_mixtures = [el for el in genetic_mixtures if el not in subsetters and el not in sub_combos]
    
    print('Assinging {}'.format(genetic_mixtures))
    matrix = substitution_matrices.load('BLOSUM90')
    calculator = DistanceCalculator('blosum90')
    dist_matrix = []
    names = [ele for ele in genetic_mixtures if ele not in subsetters]+[el for el in all_species]
    names = list(set(names))
    others = list(set(others))
    others_og = others
    others = ['Database_match_' + el for el in others if el not in sub_combos and el not in delete_species and el not in subsetters]
    print('Making final figure, phylogenetic plot')
    print('Calculate distances')
    sps = []
    maximal = 0
    rows = []
    all_fcmsa = {}
    name_to_group = {}
    for group,final_combined_multiple_sequence_alignment in final_combined_multiple_sequence_alignment_per_group.items():
        filtering = []
        for i in final_combined_multiple_sequence_alignment.index:
            if i in subsetters:#filter_them_out_of_analysis
                continue
            if i in species_in:
                filtering.append(i)
            elif i in genetic_mixtures:
                filtering.append(i)
        final_combined_multiple_sequence_alignment = final_combined_multiple_sequence_alignment.loc[filtering]
        final_combined_multiple_sequence_alignment = final_combined_multiple_sequence_alignment.loc[:, ~(final_combined_multiple_sequence_alignment.isin(['-'])).any()]#'X' was removed -> we loose too much important regions otherwise
        final_combined_multiple_sequence_alignment = final_combined_multiple_sequence_alignment[~(final_combined_multiple_sequence_alignment.isna()).any(axis=1)]
        sps = sps + list(final_combined_multiple_sequence_alignment.index)
        all_fcmsa[tuple(final_combined_multiple_sequence_alignment.index)] = final_combined_multiple_sequence_alignment
        for x in tuple(final_combined_multiple_sequence_alignment.index):
            name_to_group[x]=tuple(final_combined_multiple_sequence_alignment.index)
    for loc_1,seq1 in enumerate(names):
        seq1_name = seq1
        temp = []
        if seq1_name not in others and seq1_name in name_to_group:
            seq1 = all_fcmsa[name_to_group[seq1_name]].loc[seq1]
            seq1 = ''.join(list(seq1))
        if 'theoretical' in seq1_name or 'Database_' in seq1_name:# or 'Database_match_'+str(seq1_name) in names+others:
            found_loc = [0]*len(seq1)
            for p in all_starting_peptides:
                if p in seq1:
                    ind = seq1.index(p)
                    for r in range(0,len(p)):
                        found_loc[ind+r]=1
            seq1_adapted = ''.join([el for loc,el in enumerate(seq1) if found_loc[loc]==1])
        for loc_2,seq2 in enumerate(names):
            if seq2==seq1_name:
                temp.append(0)
                continue
            if 'theoretical' not in seq2 and 'theoretical' not in seq1_name and loc_2<loc_1:
                temp.append(dist_matrix[loc_2][loc_1])
                continue
            if seq2 in name_to_group and seq1_name in name_to_group and seq2 not in others and seq1_name not in others and name_to_group[seq1_name]==name_to_group[seq2]:
                seq2 = all_fcmsa[name_to_group[seq2]].loc[seq2]
                seq2 = ''.join(list(seq2))
                if 'theoretical' in seq1_name or 'Database_' in seq1_name:# or 'Database_match_'+str(seq1_name) in names+others:
                    seq2 = ''.join([el for loc,el in enumerate(seq2) if found_loc[loc]==1])
                    dist = thread_align_mix(seq1_adapted,seq2,matrix,calculator)
                else:
                    dist = thread_align_mix(seq1,seq2,matrix,calculator)
                if dist>maximal:
                    maximal = dist
                temp.append(dist)
            else:
                temp.append(1)#'distant'
        rows.append(seq1_name)
        dist_matrix.append(temp)
        
    print('end of calculations')
    
    
    dist_matrix = pd.DataFrame(dist_matrix,columns=names,index=rows)
    
    dist_matrix = dist_matrix.sort_index()
    dist_matrix = dist_matrix[sorted(list(dist_matrix.columns))]
    
    #we do not just do .T because we need only to change the theoreticals in relation to reals. it is basically 2 different matrices in 1 big matrix
    save_element = {}
    for element in dist_matrix.columns:
        save_element[element]=dist_matrix.loc[element].values
    for element,column in save_element.items():
        if 'theoretical' in element or 'Database_' in element:
            dist_matrix[element]=column
    
    for n in others:
        n_sp = n.replace('Database_match_','')
        dist_matrix[n_sp]=dist_matrix[n]
        dist_matrix.loc[n_sp]=list(dist_matrix[n_sp].values)+[0]
    
    
    plot_new_tree(dist_matrix,name_to_coverage,file_extinct,path,sample_path,taxonomy)
    #calculatre distance from species to species via the LCA = evolutionary distance
    #extrapolate the distances, so order genus/order/... levels are at the same distance, otherwise bias towards less documented lineages.
    #find_most common location of order
    c_order = {}
    for n in species_in:
        t = taxonomy[n]
        count = 0
        for i in t:
            if 'order' in i:
                if count in c_order:
                    c_order[count] =  c_order[count]+1
                else:
                    c_order[count] =1
                break
            count += 1
    #find_most_common location of genus
    c_genus = {}
    for n in species_in:
        t = taxonomy[n]
        count = 0
        for i in t:
            if 'genus' in i:
                if count in c_genus:
                    c_genus[count] =  c_genus[count]+1
                else:
                    c_genus[count] =1
                break
            count += 1
    #find_most common location of LCA
    in_names = [el for el in names if el in species_in]+[el.replace('Database_match_','') for el in others]
    lca_all = find_LCA(taxonomy,in_names)
    c_lca = {}
    for n in species_in:
        t = taxonomy[n]
        count = 0
        for i in t:
            if lca_all in i:
                if count in c_lca:
                    c_lca[count] =  c_lca[count]+1
                else:
                    c_lca[count] =1
                break
            count += 1
    golca = [[k for k,v in c_genus.items() if v==max(c_genus.values())][0],
             [k for k,v in c_order.items() if v==max(c_order.values())][0],
             [k for k,v in c_lca.items() if v==max(c_lca.values())][0]]
    
    transform_taxa = transform_taxons({},taxonomy,species_in,golca,lca_all)
    
    tree_matrix = []
    for n in in_names:
        begin = taxonomy[n]
        tree_row = []
        for n2 in in_names:
            if n==n2:
                tree_row.append(0)
                continue
            if dist_matrix[n][n2]==0:#if the measured distance equals to 0 than the taxonomy should be put to 0 too, otherwis the distance calculations might be biased for this theoretical species
                tree_row.append(0)
                continue
            lca = find_LCA(taxonomy,[n,n2])
            t = taxonomy[n2]
            
            dist_n_lca=transform_taxa[n][[el[1] for el in begin].index(lca)]
            
            dist_n2_lca = transform_taxa[n2][[el[1] for el in t].index(lca)]
            tree_row.append(dist_n_lca+dist_n2_lca)
        tree_matrix.append(tree_row)
    
    #normalize the tree with measured data
    #It can be that some distant species have the same detected peptides. This means that these need to be equalized.
    #So no difference possible to detect for this sample with the given database
    tree_matrix = pd.DataFrame(np.array(tree_matrix), columns=in_names, index=in_names)
    
    for i in tree_matrix.index:
        for q in tree_matrix.columns:
            if dist_matrix[q].loc[i]==0:
               tree_matrix[q].loc[i]=0
    
    #We can calculate the the distance intervals from each species in the dataset towards the lca with another species
    species_intervals = find_species_intervals(in_names,tree_matrix, transform_taxa,dist_matrix)
    
    
    #now calculate the location and the distance of new sequences with the database.
    #the theoretical is between taxonX+extrapolated_taxonX+1
    
    for i in genetic_mixtures:
        th_dist = dist_matrix[i][dist_matrix[i]<1].sort_values()
        th_dist = {el:th_dist.values[num]for num,el in enumerate(list(th_dist.index))}
        new_tree = []
        closest_in_tree = [el for el in th_dist.keys() if el in in_names]
        if len(closest_in_tree)>0:
            closest_in_tree = closest_in_tree[0]
        else:
            print('{} is very distant from everything'.format(i))
            th_dist = dist_matrix[i][dist_matrix[i]].sort_values()
            th_dist = {el:th_dist.values[num]for num,el in enumerate(list(th_dist.index))}
            closest_in_tree = [el for el in th_dist.keys() if el in in_names][0]
        all_names = []
        for n in in_names:
            if n not in th_dist:#take over distance when species are beyond order level and adjust for better results
            
            #look for closest in dist_matrix 
                distance = th_dist[closest_in_tree]
                extrapolate = species_intervals[closest_in_tree]
                for k,v in extrapolate.items():
                    if v[0]<=distance<=v[1]:
                        d_close_to_n = tree_matrix[n][closest_in_tree]
                        closest2i2 = [el for el in th_dist.keys() if el in in_names]
                        loc=0
                        if len(closest2i2)>1:
                            while loc < len(closest2i2)-1 and th_dist[closest2i2[loc]]==th_dist[closest2i2[0]]:
                                loc+=1
                        closest2_in_tree = closest2i2[loc]
                        
                        distance2 = th_dist[closest2_in_tree]
                        extrapolate2 = species_intervals[closest2_in_tree]
                        ks2 = [0]+list(extrapolate2.keys())
                        d_close2_th = [k2-((k2-ks2[loc2-1])*((extrapolate2[k2][1]-distance2)/(extrapolate2[k2][1]-extrapolate2[k2][0]))) for loc2,k2 in enumerate(ks2) if k2 != 0 and extrapolate2[k2][0]<=distance2<=extrapolate2[k2][1]][0]
                        d_close2_to_close = tree_matrix[closest2_in_tree][closest_in_tree]
                        #No we can determine the distance between theoretical and very distant
                        
                        new_tree.append(d_close_to_n+d_close2_th-d_close2_to_close)
                        all_names.append(n)
                        break
                    elif distance==0:
                        new_tree.append(distance)
                        all_names.append(n)
                        break
                    
            else:
                distance = th_dist[n]
                extrapolate = species_intervals[closest_in_tree]
                ks = [0]+list(extrapolate.keys())
                for loc,k in enumerate(ks):
                   
                    if k==0:
                        continue
                    v=extrapolate[k]
                    if v[0]<=distance<=v[1]:
                        all_names.append(n)
                        if k!=0:
                            new_tree.append(k-((k-ks[loc-1])*((v[1]-distance)/(v[1]-v[0]))))
                        else:
                            new_tree.append(k)
                        break
                    elif distance==0:
                        all_names.append(n)
                        new_tree.append(distance)
                        break
        tree_matrix[i]=new_tree
    
    #now we find the distances between the new sequences
    adding =  []
    for i in genetic_mixtures:
        temp = list(tree_matrix[i].values)            
        for i2 in genetic_mixtures:
            if i==i2:
                temp.append(0)  
            else:
                closest2i = [el for el in dist_matrix[i].sort_values().index if el in tree_matrix.index][0]
                closest2i2 = [el for el in dist_matrix[i2].sort_values().index if el in tree_matrix.index][0]
                temp.append(tree_matrix[i][closest2i2]+tree_matrix[i2][closest2i]-tree_matrix[closest2i][closest2i2])
        adding.append(temp)
    adding = pd.DataFrame(np.array(adding),columns=tree_matrix.columns, index=genetic_mixtures)
    tree_matrix=pd.concat([tree_matrix,adding],ignore_index=True)
    tree_matrix.index=tree_matrix.columns
    
    # New locations can be linked to missing taxa
    #do not to adjust for location of the theoreticals
    taxons = []
    for a in in_names:
        for b in in_names:
            if a==b:
                continue
            taxons.append(find_LCA(taxonomy,[a,b]))
    taxons=list(set(taxons))
    additional = {}
    already_found = []  
    
    all_orders = []
    for i_th in genetic_mixtures:#for all the species we found from database
        for i in [el for el in tree_matrix[i_th].sort_values().index if el in in_names and el not in already_found]:
            already_found.append(i)
            t_lin = taxonomy[i]#get taxonomic lineage
            #find children nodes for each considered lineage uptil the order
            order_loc = [el[0] for el in t_lin]
            if 'order' in order_loc:
                order_loc = order_loc.index('order')
            else:
                order_loc = 1
            for find_child in t_lin[:order_loc+1]:#Include order for example giraffa
                if find_child in already_found:
                    continue
                if 'order' in find_child:
                    all_orders.append(find_child[1])
                f_child = find_child[1]
                children = find_child_nodes(taxonomy_missing,f_child)
                children = [el for el in children if 'sp.' not in el and el.isdigit()==False]
                children = [el for el in children if f_child in [ta[1] for ta in taxonomy_missing[el]]]
                missing = background_check_mix(children,in_names,taxonomy_missing)
                additional_add = [el for el in missing.keys() if 'unclassified' not in el]
                additional[f_child]=additional_add#only keep 1 level lower
                already_found.append(find_child)
    
    not_missing_taxa = []
    for sp1 in tree_matrix.columns:
        for sp2 in tree_matrix.columns:
            if sp1 in in_names and sp2 in in_names and sp1!=sp2:
                lca = find_LCA(taxonomy, [sp1,sp2])
                if lca in tree_matrix.columns:#apparently there are cases like 'subspecies', 'Pezoporus flaviventris' that has as species 'Pezoporus wallicus' and this one is also in our database. 
                    continue
                t = taxonomy[sp1]
                t_order = [el[0] for el in t]
                if any('order' in element for element in t):
                    t_loc = t_order.index('order')
                else:
                    t_loc = 1
                if any(lca in element for element in t[:t_loc+1]):
                    if lca in additional:
                        if len(additional[lca])==0:
                            additional[lca] = [sp1,sp2]
                            not_missing_taxa.append(lca)
                    else:
                        additional[lca]=[sp1,sp2]
                        not_missing_taxa.append(lca)
    
    
    all_missing = []
    for m,mv in additional.items():
        mv = [el for el in mv if el.count(' ')<=1]#if ' ' not in el
        additional[m]=mv
        if m in all_orders:
            continue
        if len(mv)!=0:
            all_missing.append(m)
            for t in additional[m]:
                if t.count(' ')>1:#missing genera included, subspecies not relevant for this
                    continue
                try:
                    tax = taxonomy[t]
                except:
                    continue
                if any(lca_all in el for el in tax)==False:#double taxa so will never match as it should
                    continue
                all_missing.append(t)
                taxonomy[t] = tax
            if t in taxonomy:
                tree_matrix[m]=[100]*len(list(tree_matrix.index))
                add = pd.DataFrame([[100]*len(tree_matrix.columns)], index=[m], columns=tree_matrix.columns)
                tree_matrix = pd.concat([tree_matrix,add])
            mv = [el for el in mv if el in all_missing]
            additional[m]=mv
        else:
            continue
    for element in all_orders:
        print(f'{element} is order level and will not be considered for the distance calculations')
        additional.pop(element)
        
    for el in genetic_mixtures:
        if 'Database_match' in el:
            taxonomy[el]=taxonomy[el.replace('Database_match_','')]
    
    normalize_missing = {tuple(v):k for k,v in additional.items()}
    transform_taxa = transform_taxons(transform_taxa,taxonomy,all_missing+all_orders,golca,lca_all)
    for t,t_lower in additional.items():
        if len(t_lower) == 0:
            continue
        adding_missing = []
        t_og = t
        t = t_lower[0]
        for i in tree_matrix.columns:
            if i in in_names:#we have taxonomy#or 'Database_match_' in i
                # i = i.replace('Database_match_','')
                t_lin = taxonomy[i]
                t_lin = [el[1] for el in t_lin]
                t_lin2 = taxonomy[t]
                for x in t_lin2:
                    if x[1] in t_lin:
                        l1 = transform_taxa[i][t_lin.index(x[1])]
                        l2 = transform_taxa[t][t_lin2.index(x)]
                        t_lin2 = [el[1] for el in t_lin2]
                        walk_back_missing = transform_taxa[t][t_lin2.index(t_og)]-transform_taxa[t][t_lin2.index(t)]#transform_taxa[t][t_lin2.index(t_og)]#we do not have the species so distance different levels need to be taken into account
                        adding_missing.append(l1+l2-walk_back_missing)
                        break
            elif i not in all_missing and i not in additional.keys():#matching theoretical, no walk back because it is incuded in extrapolation
                closest2i = [el for el in tree_matrix[i].sort_values().index if el in in_names][0]
                closest2i2 = [el for el in tree_matrix[i].sort_values().index if el in in_names]
                loc=0
                if len(closest2i2)>1:
                    loc = 1
                    while loc < len(closest2i2)-1 and tree_matrix[closest2i2[loc]][i]==tree_matrix[closest2i2[0]][i]:
                        loc+=1
                closest2i2 = closest2i2[loc]
                t_lin = taxonomy[closest2i]
                t_lin = [el[1] for el in t_lin]

                t_lin2 = taxonomy[t]
                for x in t_lin2:
                    if x[1] in t_lin:
                        missing2i = transform_taxa[closest2i][t_lin.index(x[1])]+transform_taxa[t][t_lin2.index(x)]
                        break
                t_lin = taxonomy[closest2i2]
                t_lin = [el[1] for el in t_lin]
                
                for x in t_lin2:
                    if x[1] in t_lin:
                        missing2i2 = transform_taxa[closest2i2][t_lin.index(x[1])]+transform_taxa[t][t_lin2.index(x)]
                        break
                t_lin2 = [el[1] for el in t_lin2]
                walk_back = transform_taxa[t][t_lin2.index(t_og)]-transform_taxa[t][t_lin2.index(t)]
                
                if missing2i2==tree_matrix[closest2i][closest2i2]:#same dstance as the 1st hit
                    sc = abs((0.5*tree_matrix[closest2i][i])-tree_matrix[closest2i][closest2i2]+(0.5*tree_matrix[closest2i2][i])+(0.5*missing2i)+(0.5*missing2i2))
                else:
                    sc = abs(tree_matrix[closest2i2][i]-tree_matrix[closest2i][closest2i2]+missing2i)-walk_back
                adding_missing.append(abs(sc))
            
            else:#missing VS missing
                try:
                    i = additional[i][0]
                except:
                    print(f'Could not assign {i}')
                    adding_missing.append(0)
                    continue
                if t==i or [value for key,value in normalize_missing.items() if t in key]==[value for key,value in normalize_missing.items() if i in key]:
                    adding_missing.append(0)
                    continue
                
                t_lin1 = taxonomy[t]
                t_lin2 = taxonomy[i]
                for x in t_lin1:
                    if x in t_lin2:
                        l1 = transform_taxa[t][t_lin1.index(x)]
                        t_lin1 = [el[1] for el in t_lin1]
                        walk_back_missing1 = 0#transform_taxa[t][t_lin1.index(t)]
                        l2 = transform_taxa[i][t_lin2.index(x)]
                        t_lin2 = [el[1] for el in t_lin2]
                        walk_back_missing2 = 0#transform_taxa[i][t_lin2.index(i)]
                        adding_missing.append(l1+l2+walk_back_missing1+walk_back_missing2)#minus 1 because difference at lca not below
                        break
        
        tree_matrix[t_og]=adding_missing
        tree_matrix.loc[t_og]=adding_missing
        
    all_found = [(el,name_to_coverage[el]) for el in tree_matrix.columns if 'theoretical' in el or 'Database_' in el]
    all_found = sorted(all_found, key=lambda x:x[1])[::-1]#we want the top hits first
    all_found = [el[0] for el in all_found]
    sub_titles=['Taxonomic closeness of: '+el +' BC-score = {}'.format(name_to_coverage[el]) for el in all_found]
    
    if len(all_found)>1:
        minimal = 1/(len(all_found)-1)
    else:
        minimal = 0.025#we do not need spacing when only 1 plot
    fig = make_subplots(rows=len(all_found), cols=1, 
                        subplot_titles=sub_titles,shared_yaxes=False,
                        row_heights=[len(all_found)/100]*len(all_found))#,vertical_spacing = max(0.02, min(minimal,0.02))
    plot_nr=0
    species_to_taxon,transform_taxa = find_distance_species_taxon(transform_taxa,[el for el in tree_matrix.columns if el not in all_found],additional,taxonomy,golca,lca_all)
    for found in all_found:
        plot_nr +=1
        temp_tree_matrix = tree_matrix[tree_matrix.index.isin(all_found)==False]
        ranked = temp_tree_matrix.sort_values(by=[found])
        sp_close = [el for el in ranked.index if 'theoretical' not in el and 'Database_' not in el and ' ' in el][0]
        linkage = [el[0] for el in taxonomy[sp_close]]
        order_name = [el for el in linkage if 'order' == el]+[el for el in linkage if 'order' in el]
        if len(order_name)>0:
            order_name = order_name[0]
        else:
            order_name = linkage[-1]
        taxon_name_order = [el[1] for el in taxonomy[sp_close] if el[0]==order_name][0]

        order_level = transform_taxa[sp_close][linkage.index(order_name)]+0.5
        if len(ranked[ranked[found]<order_level])>25:
            ranked = ranked[ranked[found]<order_level]#not relevant to plot anyways
        else:
            ranked = ranked.loc[ranked.index[:30]]
        lca_loc = find_LCA(taxonomy, [el for el in ranked.index[:10]])
        group_best_match = [el[1] for el in taxonomy[ranked.index[0]] if el[1] in combined_file_output.keys()]
        if lca_loc in ranked.index:
            lca_loc = max(30,list(ranked.index).index(lca_loc))
        else:
            lca_loc = 30
        y_labs = [el for el in ranked.index[:lca_loc] if el != found]#Top 20 for each of the outcomes 
        
        temp_data,colours,main,side = main_side_lineages(ranked, found,y_labs,species_to_taxon,taxonomy,lca_loc,golca,taxon_name_order)
        
        #assing main and side lineages
        cell = [[el[1] for el in taxonomy[main[0]] if el[1] in y_labs][::-1]]#LCA to species
        for k,v in side.items():
            add_cells = []
            if len(v)==0:
                continue
            if k=='Unliked':
                cell.append(v+['']*(len(cell[0])-len(v)))
                continue
            for i in cell[0]:
                if i==k:
                    add_cells.append(k)
                    break
                else:
                    add_cells.append(i)
            add_cells = add_cells + [el for el in v if el not in add_cells]
            add_cells = add_cells + ['']*(len(cell[0])-len(add_cells))
            cell.append(add_cells)
        c_data=pd.DataFrame([[', '.join(additional[el]) if el in additional and len(additional[el])>0 and el not in not_missing_taxa else 'No lower branch missing' for el in y_labs],
                             ['Main lineage' if el in main else 'Side lineage' for el in y_labs]]).T
        
        fig.add_trace(go.Scatter(mode='markers+lines',x=temp_data['line_x'], y=temp_data['line_y'],
                                  name=found,
                                  marker=dict(symbol='arrow',color='black',angleref='previous',standoff=8),
                                              ),
                      row=plot_nr,col=1)
        fig.add_trace(go.Scatter(mode='lines+markers',x=temp_data['x'], y=temp_data['y'],
                                 name=found,customdata=c_data,marker_color='red',
                                 hovertemplate =
                                             "<b>%{y}</b><br>" +
                                             "Distance: %{x:,.4f}<br>" +
                                             "Missing under taxon: %{customdata[0]}<br>"+
                                             "Lineage type: %{customdata[1]}<br>",
                                             ),
                      row=plot_nr,col=1)
        fig.add_trace(go.Scatter(mode='markers',x=temp_data['x_adj'], y=temp_data['y'],
                                 name=found,customdata=c_data,marker_color = 'black',
                                 hovertemplate =
                                             "<b>%{y}</b><br>" +
                                             "Distance: %{x:,.4f}<br>" +
                                             "Missing under taxon: %{customdata[0]}<br>"+
                                             "Lineage type: %{customdata[1]}<br>",
                                             ),
                      row=plot_nr,col=1)
        most_likely,most_likely_species = find_most_likely(taxonomy,temp_data['y'],temp_data['x_adj'],transform_taxa,order_level-0.5)
        fig.add_trace(go.Scatter(mode='markers',x=most_likely, y=most_likely_species,
                                 name=found,customdata=c_data,
                                 marker=dict(color='royalblue'),marker_symbol='hexagram-open',marker_size=12,
                                             ),
                      row=plot_nr,col=1)
        
        all_match = []
        scrambled=False
        in_ranked = list(ranked[ranked[found]<order_level-0.5].index)
        for location, element in enumerate(in_ranked):
            now = [el[1] for el in taxonomy[element] if el[1] in combined_file_output.keys()]
            if any(el in all_match for el in now)==False and len(all_match)!=0:
                if len(in_ranked)==location+1:
                    break
                if all_match[0] in in_ranked[location+1:]:
                    scrambled = True
            all_match = now
        if min(temp_data['x_adj'])<=3 and len(most_likely)>0 and scrambled==False:
            colorgr = 'lightgreen'
        else:
            colorgr = 'firebrick'
        orange_flagged = [temp_data['y'][loc] for loc,el in enumerate(temp_data['x']) if el<=1 and temp_data['x_adj'][loc]<=1 and temp_data['y'][loc] not in all_orders]
        if len(orange_flagged)>0:
            if find_LCA(taxonomy,orange_flagged) in all_orders:
                colorgr = 'darkorange'
        for i in range(0,int(max(temp_data['x_adj']))+1):
            fig.add_vrect(
                x0=i, x1=i+1,
                fillcolor=colorgr, opacity=max(0,0.5-i/10),
                layer="below", line_width=0,
                row=plot_nr,col=1
            ),
            
        fig.update_xaxes(title_text="Distance", row=plot_nr, col=1)
    fig.update_layout(title_text='Distance of potential candidates to Top20 taxa, including taxa not in database of sample :{}'.format(file_extinct),
                  height=800*len(all_found), showlegend=False, template='plotly_white',
                 )

    name_file='Ranked_'+file_extinct+'.html'
    try:
        fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'taxonomic_output' /name_file)
    except:
        fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'taxonomic_output' / 'Ranked_classification.html') 
    
    
    seqs = {}
    for g,fcmsa in fcmsa_per_group.items():
        gm = [el for el in genetic_mixtures if el in fcmsa.index and el in tree_matrix.columns]
        if len(gm)>0:
            seqs = seqs | recover_sequences(tree_matrix,fcmsa,sequences,gm,[el for el in tree_matrix.columns if el not in genetic_mixtures],TF_per_group[g],subsetters,total_insilico_all[g])
    return seqs, sps

def ClassiCOL_mixture_analysis(
    path: pathlib.Path,
    sample_path: str,
    file_name: str,
    sequences: dict[Seq, str],
    cpu_count: int,
    Mixture_analysis: str,
    taxonomy: dict[str,list]
    ):
    mixture = True if Mixture_analysis=='M' else False
    #Take the taxon result file generated before 
    name_file='Taxonomic_results_after_rescoring_'+file_name+'.csv'
    
    output_sequences = {}
    all_seqs = {}
    all_sp_in = []
    print('Analysing sample for mixtures and species not in the database.')
    path_to_output_file = path / 'Output_Classicol' / sample_path
    file_extinct = sample_path
    combined_file_output = {}
    missed_outcome = []
    print('Assigning taxonomic groups under each order level.')
    df, taxonomy, taxonomic_groups,species_in,rank_lca,restrict,taxonomy_missing = parse_input(path_to_output_file, name_file,mixture,sequences,taxonomy)
    df.columns = ['Peptides' if el=='Peptide' else el for el in df.columns]
    all_starting_peptides = set(df['Peptides'].values)
    print('Making consensus sequences for all order level candidates.')
    for k,v in sequences.items():
        all_seqs[k]=v
    #Either make new if new species to database otherwise load saved
    total_consensus_df_all,total_insilico_all,single_outcome,asc,cpaf = do_consensus(df, sequences, taxonomy, taxonomic_groups,restrict,path,sample_path,species_in)
    df = df[df['Protein'].isin(asc)]
    print('Consensus sequences are made and saved.')
    if len(total_insilico_all)==0:
        print('No mixture possible')
        rev_s = {v:k for k,v in sequences.items()}
        output_sequences[file_extinct]=[rev_s[protein] for protein in list(set(df['Protein'].values)) if list(df['Protein'].values).count(protein)>20]
        return output_sequences, missed_outcome, all_starting_peptides
    group_nr=0
    theoretical_trace_back = {}
    for x in single_outcome:
        theoretical_trace_back['Database_match_'+x]=(x,x)#for if another branch has 1 candidate, distant beyond order level
    keep_groups = delete_groups(taxonomic_groups, df, path, sample_path,file_extinct,cpaf)
    taxonomic_groups = {k:v for k,v in taxonomic_groups.items() if k in keep_groups}
    for groups in taxonomic_groups.keys():
        if any(groups in el for el in total_consensus_df_all.keys())==False:
            for i in taxonomic_groups[groups]:
                missed_outcome.append(i)
                print('Not enough evidence for {} to be included'.format(i))
            continue
        total_consensus_df = {key:val for key,val in total_consensus_df_all.items() if groups in key}
        total_insilico = {key:val for key,val in total_insilico_all.items() if groups in key}
    
        taxa_in_all = set(total_consensus_df[list(total_consensus_df.keys())[0]].index)
        for i in total_consensus_df.keys():
            taxa_in_all = taxa_in_all&set(total_consensus_df[i].index)
        for i in total_consensus_df.keys():
            total_consensus_df[i]=total_consensus_df[i].loc[list(taxa_in_all)]
        #Now we have the consensus sequence and all the species sequences mapped to that one
        taxonomy_missing = taxonomy_missing | {k:v for k,v in taxonomy.items() if k not in taxa_in_all and k in species_in}#All that have been filtered out become missing taxa
        taxonomy = taxonomy|taxonomy_missing
        pep_locs_all, mixtures_all = initial_iteration(total_consensus_df, total_insilico,df,taxonomy,cpu_count,species_in)
            
        #many VS many
        stop_making_mixtures=False
        final_output_genetic_mixtures = {}
        for el in set(df['Species'].values):
            theoretical_trace_back['Database_match_'+el]=(el,el)
        
        animal_trace_back = {}
        physical_mix_trace_back ={}
        already_compared = []
        while stop_making_mixtures==False:
            print('Start iterating ...')
            new_theoretical = {}
            physical_mix = {}
            search_space = list(total_consensus_df[list(total_consensus_df.keys())[0]].index)
            taxon_contribution_score = find_potential_missing_species(search_space,mixtures_all,pep_locs_all,taxonomy,theoretical_trace_back,already_compared)
            lcas = {el:''.join(el.split('!')[0]) for el in taxon_contribution_score.keys()}
            lcas = {el:find_LCA(taxonomy,lcas[el].split('_VS_')) for el in lcas.keys()}
            no_difference = [el for el,val in taxon_contribution_score.items() if (len(val)==0 or [lcas[el],1] in val or el.split('_VS_')[0] == el.split('!')[0].split('_VS_')[1]) and '$' not in el]
            already_compared = already_compared + [sorted(el.split('!')[0].split('_VS_')) for el in taxon_contribution_score.keys()]
            taxon_contribution_score = {k:v for k,v in taxon_contribution_score.items() 
                                        if k not in no_difference}
            removal = []
            one_side_keep = []
            
            
            ##Check if 2 species are exactly the same, and if so add 1 to the removal pile
            double_remove = []
            for k in no_difference:
                pm = k.split('!')[-1].replace('!MIX','')
                k = k.split('!')[0]
                print('No overlap or complementarity considering: {}'.format(k))
                if k in double_remove:
                    continue
                if k.split('_VS_')[0]!=k.split('_VS_')[1]:
                    final_output_genetic_mixtures[tuple(k.split('_VS_'))]=0
                rem_keep = sorted(k.split('_VS_'))
                if pm != '0':
                    removing = checking_physical_mix(k,theoretical_trace_back,physical_mix_trace_back)
                    if removing != False:
                        removal = removal+removing
                        if rem_keep[0] not in removing:
                            one_side_keep.append(rem_keep[0])
                        else:
                            one_side_keep.append(rem_keep[1])
                        double_remove.append(k)
                        continue
                
                #We keep 1 candidate is multiple are equal, otherwise theorethicals that are equal might be compared
                if rem_keep[0] in removal or rem_keep[0] in one_side_keep or rem_keep[1] in removal or rem_keep[1] in one_side_keep:
                    if rem_keep[0] not in one_side_keep:
                        removal.append(rem_keep[0])
                        double_remove.append(k)
                    if rem_keep[1] not in one_side_keep:
                        double_remove.append(k)
                        removal.append(rem_keep[1])
                    continue
                double_remove.append(k)
                removal.append(rem_keep[0])#only remove 1 side for computational purposes!!
                one_side_keep.append(rem_keep[1])
            #look for subsetters, and if any are found, remove
            remove_from_tcs=[]
            
            for k,v in taxon_contribution_score.items():
                if k.count('$')>0:#LCA almost always higher when full mixture, half mixture will be dealed with in next loop
                    continue
                pm = k.split('!')[-1].replace('!MIX','')
                k=k.split('!')[0]
                #check for chance of genetic mix
                lca = find_LCA(taxonomy,k.split('_VS_'))
                k1,k2 = path_to_lca(taxonomy,lca,k)
                
                weighted_taxa = [el[0] for el in v]
                if len(set(k1)&set(weighted_taxa))==0 or len(set(k2)&set(weighted_taxa))==0:
                    final_output_genetic_mixtures[tuple(k.split('_VS_'))]=0
                    print('Only uniqueness for 1 side considering: {}'.format(k))
                    remove_from_tcs.append(k)
                    if len(set(k1)&set(weighted_taxa))==0:
                        removal.append(k1[0])
                        one_side_keep=[el for el in one_side_keep if el!=k1[0]]
                        one_side_keep.append(k2[0])
                        one_side_keep=[el for el in one_side_keep if el!=k1[0]]
                    else:
                        removal.append(k2[0])
                        one_side_keep=[el for el in one_side_keep if el!=k2[0]]
                        one_side_keep.append(k1[0])
                    continue
            taxon_contribution_score = {k:v for k,v in taxon_contribution_score.items() if any(el in k for el in remove_from_tcs)==False} 
            mixture_one_side_keep = []
            
            
            #remove the physical mixtures from the convergence
            del_from_tcs = []
            all_total_mix = [el for el in taxon_contribution_score.keys() if '!TOTAL' in el]
            mixed_keep = []
            temp_dict = {}
            #Now figure out which ones are actually mixed
            for k in all_total_mix:
                remember=k
                m_keep = False
                v = taxon_contribution_score[k]
                side = [int(k.split('$')[-1].split('!')[-1])][0]-1
                k = k.split('!TOTAL_MIX')[0]
                x = k.split('_VS_')[side]
                del_from_tcs.append(remember)             
                print('Potential physical mix of {}, genetic mix will also be tested'.format(remember))
                t = remember.replace('TOTAL_MIX$','').split('!')
                lca = find_LCA(taxonomy, t[0].split('_VS_'))
                add = []
                for tk in k.split('_VS_'):
                    if 'theoretical' in tk:
                        add = add+animal_trace_back[tk]
                    else:
                        add.append(tk)
                if t[-1] in t[:-1] and sorted(t[0].split('_VS_')) not in theoretical_trace_back.values() and remember.count('TOTAL_')==1 and lca!=groups:
                    splits = t[0].split('_VS_')
                    print('Mixing {}'.format(splits))#mixing of true mixtures only allowed up to 6 species, else too random
                    t=t[0]+'!MIX!'+t[-1]
                    if x not in removal:
                        #keep the ones with enough evidence
                        one_side_keep.append(x)
                    if len(set(find_trace(splits[0],theoretical_trace_back))|set(find_trace(splits[1],theoretical_trace_back)))>5 or lca==groups:
                        print('Cancelled mixing {}, because the mix indicates something away from {}'.format(k.split('_VS_'),groups))
                    else:
                        temp_dict[t]=v
                elif t[-1] in t[:-1] and sorted(t[0].split('_VS_')) not in theoretical_trace_back.values() and remember.count('TOTAL_')==1:
                    if x not in removal:
                        print('Keeping the unmixable {}'.format(t[0].split('_VS_')))
                        #keep the ones with enough evidence
                        one_side_keep.append(x)
                elif remember.count('$')==2:#need for separation, and keep them
                    k1,k2 = path_to_lca(taxonomy,lca,k)
                    splits = t[0].split('_VS_')
                    vals = [np.sum(np.array([el[1] for el in v if el[0] in k1])), np.sum(np.array([el[1] for el in v if el[0] in k2]))]
                    if vals[side]>0.9 and lca!=groups and len(set(find_trace(splits[0],theoretical_trace_back))|set(find_trace(splits[1],theoretical_trace_back)))<=5:#means that it can be both, because a total mix will be maximum lower than 0.9
                        t=t[0]+'!MIX!'+t[-1]
                        temp_dict[t]=v
                        print('Close related species will be mixed')
                    else:
                        print('Close related species show too much difference, not mixable')
                    if vals[side]>0.2:#if the comparison is too low scoring for this side, than it is probably a false positive, so not include anymore
                        m_keep = True
                elif sorted(t[0].split('_VS_')) in theoretical_trace_back.values() and remember.count('$')==2:#already did it last time and they are still possible, so need to keep it individually, otherwise it will be discarded >false negative
                    print('Already mixed, keeping individual candidates {}'.format(t[0].split('_VS_')))
                    m_keep = True
                else:
                    print('Too much difference to mix {}'.format(t[0].split('_VS_')))
                    removal = removal + [x] 
                    one_side_keep = [el for el in one_side_keep if el!= x]
                    mixed_keep = [el for el in mixed_keep if el!= x]
                    print('Mix does not make sense')
                if m_keep ==True and x not in removal:
                    mixed_keep.append(x)
            taxon_contribution_score = {k:v for k,v in taxon_contribution_score.items() if k not in del_from_tcs}
            taxon_contribution_score = taxon_contribution_score|temp_dict
            considered_groups = list(set([el for el in removal]))+list(set(one_side_keep))+list(set(mixed_keep))
            ##Build new theoretical sequences and save them for traceback, note the no_difference candidates that were kept are added here too
            keep = []
            for k,v in taxon_contribution_score.items():
                remember = k
                pm = k.split('!')[-1].replace('!MIX','')
                k=k.split('!')[0]
                #if a genetic mix overlays with another, we have a physical mix
                #add VS to physical mix
                #remove the VS from the possibility list, so they are not combined in any way
                #take the k vals into account for this, low on one side can mean physical still possible
                lca = find_LCA(taxonomy,k.split('_VS_'))
                temp_lca = lca
                lca_comb = lca
                k1,k2 = path_to_lca(taxonomy,lca,k)
                k1vals = np.sum(np.array([el[1] for el in v if el[0] in k1]))
                k2vals = np.sum(np.array([el[1] for el in v if el[0] in k2]))
                
                key_tot = list(total_consensus_df.keys())[0]
                while temp_lca not in list(total_consensus_df[key_tot].index):# and 'theoretical' not in ''.join(k.split('_VS_')):
                    temp_lca = find_LCA(taxonomy,[temp_lca,k.split('_VS_')[0]],True)
                lca = temp_lca

                klca = np.sum(np.array([el[1] for el in v if el[0] in [lca,lca_comb]]))
                if k1vals+k2vals<klca and klca>0.8 and pm!='0':#only select when many potential candidates, the more candidates the more random the combinations so only the ones scoring higher than the LCA are good
                    final_output_genetic_mixtures[tuple(k.split('_VS_'))]=0
                    print('Combination not valid to make: {}'.format(k))
                    new_theoretical = {k:v for k,v in new_theoretical.items() if k!=remember}
                    continue
                elif (k1vals+k2vals<klca and klca<=0.6 and '$' in remember) or remember.count('$')>1:
                    if remember.count('$')>1:
                        print('Potential mixture of {}'.format(k))
                        keep = keep + [k1[0],k2[0]]
                    elif '$'+pm in remember:
                        keep = keep + [k.split('_VS_')[int(pm)-1]]
                    continue
                add = []
                for tk in k.split('_VS_'):
                    considered_groups.append(tk)
                    if 'theoretical' in tk:
                        add = add+animal_trace_back[tk]
                    else:
                        add.append(tk)
                add=list(set(add))
                #only make a new theorethical if less than 6 species represented. Otherwise too random
               
                    
                #mixture => mix can be considered 2 time, with different residues to contribute
                considered_temp = list(set(one_side_keep)) + list(set(mixed_keep))+list(set(keep))
                if (len(set(find_trace(k.split('_VS_')[0],theoretical_trace_back))|set(find_trace(k.split('_VS_')[1],theoretical_trace_back)))>6 or find_LCA(taxonomy,k.split('_VS_'))==groups or any((k.split('_VS_')[0] in element_v or k.split('_VS_')[1] in element_v) and 'theoretical_'+str(element_k) in considered_temp for element_k,element_v in theoretical_trace_back.items())) and int(pm)!=0 and (('theoretical' in k.split('_VS_')[0] or 'theoretical' in k.split('_VS_')[1]) and len(set(find_trace(k.split('_VS_')[0],theoretical_trace_back))|set(find_trace(k.split('_VS_')[1],theoretical_trace_back)))>3):
                    if any((k.split('_VS_')[0] in element_v or k.split('_VS_')[1] in element_v) and 'theoretical_'+str(element_k) in considered_temp for element_k,element_v in theoretical_trace_back.items()):
                        #if species already donated to something else, it is not allowed to be mixed untill that last one is dissolved. otherwise combos within physical mixtures will go too far up the tree
                        print('Any of {} already donated peptides to something that is still within the list of possibilities.\n So it is not mixable anymore.'.format(k.split('_VS_')))
                    elif find_LCA(taxonomy,k.split('_VS_'))==groups:
                        print('Mixing {}, cancelled because the mix indicates something away from {}'.format(k.split('_VS_'),groups))
                    else:
                        print('Mixing {}, cancelled because too much combos made'.format(k.split('_VS_')))
                    for element in k.split('_VS_'):
                        if element not in removal and ((k2vals>k1vals-0.1 and int(pm)==2) or (k1vals>k2vals-0.1 and int(pm)==1)):
                            keep.append(element)
                        elif element in removal or ((k2vals<=k1vals-0.1 and int(pm)==2) or (k1vals<=k2vals-0.1 and int(pm)==1)):
                            print('{} removed from possibilities, because it did not match criteria'.format(element))
                    continue
                
                pm_adj = int(pm)
                if pm_adj==0:
                    print('Mixing {}, because total complement found completely donating from both sides'.format(k.split('_VS_')))
                elif pm_adj==2:
                    if k2vals<k1vals-0.1 or k1vals == 0:
                        print('Scoring determined that mixing {} towards the right side does not make sense'.format(k.split('_VS_')))
                        if  k1vals == 0:
                            print('Potential mixture of {}'.format(k2[0]))
                            keep = keep + [k2[0]]
                        continue
                    print('Mixing {}, because total complement found towards the right side'.format(k.split('_VS_')))
                else:
                    if k1vals<k2vals-0.1 or k2vals == 0:
                        print('Scoring determined that mixing {} towards the left side does not make sense'.format(k.split('_VS_')))
                        if  k2vals == 0:
                            print('Potential mixture of {}'.format(k1[0]))
                            keep = keep + [k1[0]]
                        continue
                    print('Mixing {}, because total complement found towards the left side'.format(k.split('_VS_')))
                group_nr+=1
                new_theoretical[group_nr]=k.split('_VS_')
                removal = removal + k.split('_VS_')#combo is made, now we can remove these sequences because they don't matter anymore
                animal_trace_back['theoretical_'+str(group_nr)]=sorted(add)
                physical_mix['theoretical_'+str(group_nr)]=pm_adj
                if remember.count('$')>1:
                    print('Potential mixture of {}'.format(k))
                    keep = keep + [k1[0],k2[0]]
                elif '$'+pm in remember:
                    keep = keep + [k.split('_VS_')[int(pm)-1]]

            considered_groups = list(set(list(set([el for el in removal]))+considered_groups+keep+list(set(mixture_one_side_keep))+mixed_keep))
            if len(new_theoretical)==0 and len(keep+mixed_keep+one_side_keep)==0:#Nothing to check anymore
                print('Stopping because nothing is left')
                if len(set(removal)^set(considered_groups))>0:#IF 1 oucome, meaning 1 species/mixture, we need to output this
                    for x in set(removal)^set(considered_groups):
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                else:
                    for x in keep+one_side_keep+mixed_keep+mixture_one_side_keep:
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                stop_making_mixtures = True
                continue
            
            theoretical_trace_back = theoretical_trace_back|new_theoretical#remember for traceback
            physical_mix_trace_back = physical_mix_trace_back|physical_mix#remember for traceback
            
            considered_species = []
            for x in new_theoretical.values():#add all that remain
                considered_species=considered_species+x
            #add the theoreticals to the 1VS1 to get a 1VSmany and manyVSmany
            considered_species= list(set(one_side_keep)) +[
                'theoretical_'+str(k) for k in new_theoretical.keys()]+ list(set(mixed_keep))+list(set(keep))
            
            #add the theoretical taxonomy to the real_taxonomy
            for k,v in new_theoretical.items():
                k = 'theoretical_'+str(k)
                lca = find_LCA(taxonomy,v)
                taxonomy[k]=[('species',k)]+taxonomy[lca]
            #add now the new_theoretical_alignment to the df_c
            total_consensus_df,doubled_sequences = add_to_consensus(total_consensus_df,new_theoretical,total_insilico,physical_mix)
            identicals = []
            for idents in doubled_sequences.values():
                identicals = identicals + idents
            considered_species = [el for el in considered_species if el not in identicals] #we do not want to include sequences that are double or already analysed before  

            if len(identicals)>0:
                print('Removed the following theorethicals {}'.format(', '.join(list(set(identicals)))))
            considered_species = list(set(considered_species))
            out_species,r,final_output_genetic_mixtures = find_combinations_of_species(considered_species,theoretical_trace_back,already_compared,taxonomy,final_output_genetic_mixtures)
            removal = removal + r
            
            # remove = ['theoretical_'+str(k) for k in new_theoretical.keys() if [el for el in set(doubled_sequences[k]) if doubled_sequences[k].count(el)>=len(total_consensus_df)]]
            # out_species = [el for el in out_species if len(set(remove)&set(el))==0]
            remove_combo = []
            out_species_extra =[]
            lca_done = []
            added = []
            print('Checking for irrelevant combinations')
            for element in out_species:#species at lca level when only few to begin with are done in the initial iteration
                lca_temp = find_LCA(taxonomy,element)
                left_done = [find_LCA(taxonomy,el) for el in already_compared if element[0] in el and len(set(el))!=1]
                right_done = [find_LCA(taxonomy,el) for el in already_compared if element[1] in el and len(set(el))!=1]
                #if a species is compared multiple times with somehting with LCA shared. Than we need to check if it makes sense to analyse all of them.
                if lca_temp == groups and (len([ele for ele in set(left_done) if ele != groups])>1 or len([ele for ele in set(right_done) if ele != groups])>1 or abs(left_done.count(groups) - right_done.count(groups))>3) and left_done.count(groups)>3 and right_done.count(groups)>3:
                    remove_combo.append(element)
                    already_compared = already_compared+[element]
                    if [element[0],element[0]] not in out_species_extra and element[0] not in added and element[0] not in lca_done:
                        out_species_extra = out_species_extra+[[element[0],element[0]]]
                    if [element[1],element[1]] not in out_species_extra and element[1] not in added and element[1] not in lca_done:#if a species makes it to the LCA (order) than it is in there most likely, no need to check against all other species 
                        out_species_extra = out_species_extra+[[element[1],element[1]]]
                    print('Dissolving {}'.format(element))
                elif lca_temp == groups and (element[0] in lca_done or element[1] in lca_done) and len([ele for ele in set(left_done) if ele != groups])>1 and len([ele for ele in set(right_done) if ele != groups])>3:
                    remove_combo.append(element)
                    already_compared = already_compared+[element]
                    if [element[0],element[0]] not in out_species_extra and element[0] not in lca_done and element[0] not in added:
                        out_species_extra = out_species_extra+[[element[0],element[0]]]
                    if [element[1],element[1]] not in out_species_extra and element[1] not in lca_done and element[1] not in added:#if a species makes it to the LCA (order) than it is in there most likely, no need to check against all other species 
                        out_species_extra = out_species_extra+[[element[1],element[1]]]
                    print('Dissolving {}'.format(element))
                elif lca_temp == groups:#no need to check at the lca level in a 1VSmany, if it is not a subset at this point than it will never be
                    lca_done = lca_done+element
                added = added + element 
            #we do not want to split a species that is already assigned before, otherwise even if it is a subset it will be retained
            out_species = [el for el in out_species if el not in remove_combo]
            out_species = out_species +out_species_extra
            out_species = [el for el in out_species if len(set(el))==1 and any(el[0] in ele for ele in out_species if len(set(ele))>1)==False]+[el for el in out_species if len(set(el))!=1]
            #if nothing to compare anymore, than we stop the run
            if (len(out_species)==0 and len([el for el in out_species if el not in already_compared])==0) or len([element for element in out_species if element[0]!=element[1]])==0:
                print('Stopping before next node allocation')
                if len(set(removal)^set(considered_groups))>0:#IF 1 outcome, meaning 1 species/mixture, we need to output this
                    for x in set(removal)^set(considered_groups):
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                else:
                    for x in set(keep+mixture_one_side_keep+mixed_keep):
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                if len([el[0] for el in out_species_extra])>0:
                    for x in set([el[0] for el in out_species_extra]):
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                for element in out_species:
                    for x in element:
                        if 'theoretical' in x:
                            x = tuple(theoretical_trace_back[int(x.split('_')[-1])])+(int(x.split('_')[-1]),)
                        else:
                            x = (x,x,0)
                        final_output_genetic_mixtures[x]=1
                stop_making_mixtures = True
                continue
            
            #remove all candidates that are no longer possible, so their 'unique' peptides do not have an influence on the scoring of the true hits.
            taxa_iteration = set()
            remove_from_df  = set()
            for i in total_consensus_df.keys():
                taxa_iteration = taxa_iteration|set([ele for ele in total_consensus_df[i].index if 'theoretical' not in ele or any(ele in elements for elements in out_species)])
                remove_from_df = remove_from_df|set([ele for ele in total_consensus_df[i].index if 'theoretical' in ele and any(ele in elements for elements in out_species)==False])
            if len(remove_from_df)>0:
                print('The following theoretical sequences were removed from the possibility list: {}'.format(remove_from_df))
            for i in total_consensus_df.keys():
                total_consensus_df[i]=total_consensus_df[i].loc[list(taxa_iteration)]
            
            mixtures_all,pep_locs_all = further_iteration(cpu_count, total_consensus_df,total_insilico,out_species,species_in,taxonomy,df)
        
        final_output_genetic_mixtures = {k:v for k,v in final_output_genetic_mixtures.items() if v==1}
        if len(final_output_genetic_mixtures)>0:
            combined_file_output[groups]=[final_output_genetic_mixtures,total_consensus_df,total_insilico]
    if len(combined_file_output)>0:
        seqs,sps = plot_taxonomy(file_extinct,theoretical_trace_back,species_in,sequences,single_outcome,taxonomy,combined_file_output,all_starting_peptides,df,path,sample_path,taxonomy_missing)#plot the theoreticals into the taxonomy, detected peptide based
        all_sp_in = list(set(all_sp_in+[el for el in sps if el in species_in]))
        output_sequences[file_extinct] = seqs
    else:
        print('Not enough evidence to look for mixtures or species not in the database!!')
    
    
    return output_sequences, missed_outcome, all_starting_peptides

def return_missed(m,path,sample_path):
    with open(path/'Output_Classicol'/sample_path/'mixture_plots' /'Species_deleted_from_mixture_analysis.txt', "a") as f:
        f.write("The following species were deleted from the mixture analysis of ClassiCOL2.0 due to a lack of peptides:\n")
        for sp in m:
            f.write(sp+'\n')
    return

def return_measured_seqs(output_sequences,path,sample_path,peptide_list):
    for title,sequences in output_sequences.items():
        if 'Taxonomic_results_after_rescoring_' in title:
            title = title.replace('Taxonomic_results_after_rescoring_','')
        else:
            title = title.replace('Taxonomic_results_after_rescoring','')
        all_seqs = {}
        peptide_list = sorted(peptide_list, key=lambda x:len(x))[::-1]#first larger peptides than short ones
        if type(sequences)==type(dict()):
            for th,res in sequences.items():
                name = th.split('_')[:-1]
                name = '_'.join(name)
                close_r = res[2]
                name2 = 'Group: '+name+' measured data and non-mutational space '+close_r
                name = 'Group: '+name+' measured data only matching: '+close_r
                seq = ''.join([el if res[3][loc]==1 else 'X' for loc,el in enumerate(res[0])])
                #Can be that some of the AA residues creep in fasta randomly, so adjust before returning
                seq_adjust = ''.join([el for el in seq])
                for p in peptide_list:
                    if p in seq_adjust:
                        seq_adjust = seq_adjust.replace(p,'!'*len(p))
                    elif p in seq:
                        first_occurance_list = [el+seq.index(p) for el in range(0,len(p))]
                        seq_adjust = ''.join([el if loc not in first_occurance_list else '!' for loc,el in enumerate(seq_adjust)])
                seq = ''.join([el if seq_adjust[loc]=='!' else 'X' for loc,el in enumerate(seq)])
                all_seqs[name]=seq
                all_seqs[name2]=res[0]
                all_seqs[close_r]=res[1]
        else:
            for i in range(0,len(sequences)):
                all_seqs[title+'_'+str(i)]=sequences[i]
        sequence_list = []
        for k,s in all_seqs.items():
            s=Seq(s)
            sequence_list.append(SeqRecord(s,k,'',''))
        
        print('writing fasta {}'.format(title))
        
        name_file = title+'.fasta'
        try:
            with open(path/'Output_Classicol'/sample_path/'mixture_plots' /'fastas'/name_file,'w') as output_handle:
                SeqIO.write(sequence_list,output_handle,'fasta')
        except:
            with open(path/'Output_Classicol'/sample_path/'mixture_plots' /'fastas'/'output.fasta','w') as output_handle:
                SeqIO.write(sequence_list,output_handle,'fasta')
    return

if __name__ == "__main__":#rescore output in csv, reduced info output file, summary output
    ##############PROMPT INPUT ARGUMENTS###############################
    parser = argparse.ArgumentParser(# -e ERROR \n   -p PEPTIDE_TABLE [PEPTIDE_TABLE] |    [-l LIMIT]\n   [-t TAXONOMY]\n   [-n NEIGHBOURING] [-a]
                                         usage="py ClassiCOL_version_2_0_0.py \n -f Path to additional Database in FASTA format (not required)\n -d DIRECTORY to classicol \n -l folderNAME with search engine output files\n -s Search engine (MASCOT, MaxQuant, manual, winnow, peaks)\n -t limitation of taxonomy (e.g. Pecora or for species: Bos_taurus or Homo_sapiens/Canis)\n -c CPUs maximum to be used \n -b Single species (S) or Mixture (M)",
                                         description="ClassiCOL species classification via peptide ambiguation")


    inputs = parser.add_argument_group('\nVariable inputs')
    inputs.add_argument("-f", dest="Manual_fasta", help="path to fasta document other than standard classicol database (ends in .fasta or .txt)", type=str)
    inputs.add_argument("-d", dest="Directory",  help="Directory where Classicol is located", type=str, required=True)
    inputs.add_argument("-s", dest="Search_engine", help="Search engine used, allowed types are: MASCOT .cvs and MaxQuant .txt", type=str, required=True)
    inputs.add_argument("-t",dest="limited_taxonomy", help="Subset of taxonomy e.g. Pecora or Pecora/Primates", type=str)
    inputs.add_argument("-l",dest="File_location", help="path to folder containing search engine files", type=str)
    inputs.add_argument("-m", dest="Fixed_modification", help="Fixed modification e.g. C,45.98/M,...", type=str)
    inputs.add_argument("-c", dest="CPUs_to_be_used", help="Amount of CPUs the algorithm can maximally use", type=str, required=True)
    inputs.add_argument("-b", dest="Single_Mixed", help="Single species sample or mixture? When not speciefied mixture is considered.", type=str)
    args = parser.parse_args()
    path = pathlib.Path(args.Directory)
    
    add_fasta = None if args.Manual_fasta is None else pathlib.Path(args.Manual_fasta)
    
    if args.limited_taxonomy == "":
        args.limited_taxonomy = None
    if args.Single_Mixed == "" or args.Single_Mixed == None: 
        Mixture_analysis = 'M'
    else:
        Mixture_analysis = args.Single_Mixed
    demo_tf: bool = args.File_location.lower() == "demo"
    location_searchfiles = path / "Demo" if demo_tf else pathlib.Path(args.File_location)
    lim_tax = "Pecora" if demo_tf and args.limited_taxonomy is None else args.limited_taxonomy
    Search_engine = args.Search_engine.lower()
    
    fixed_mod: list[tuple[str, float]] = []
    if args.Fixed_modification is not None:
        for modification_str in args.Fixed_modification.split("|"):
            amino_acid, mass = modification_str.split(',')
            fixed_mod.append((amino_acid, float(mass)))

    # Set CPU count
    cpu_count = max(1, min(int(args.CPUs_to_be_used), multiprocessing.cpu_count()))
    if cpu_count != int(args.CPUs_to_be_used):
        print(f"WARNING: The application only supports CPU counts between {1} and {multiprocessing.cpu_count()}.")
        print(f"WARNING: CPU count has been changed to {cpu_count}.")
    Search_engine = Search_engine.lower()
    if Search_engine not in ['mascot', 'maxquant', 'manual','winnow','peaks']:
        print(f"ERROR: '{Search_engine}' is not a valid search engine.")
        
    outpath = path / 'Output_Classicol'
    outpath.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    all_files_to_analyse: list[tuple[pathlib.Path, typing.Any]] = []
    sequence_db = crap_f(path, add_fasta)
    ############################################################################
    data_to_remember = pd.DataFrame(columns=['animal','protein','mascot_peptide',
                                             'found_match','switch','location','type','PTM','PTM_loc','PTM_og'])
    ############################################################################
    # Find files to analyze
    for file in location_searchfiles.iterdir():
        match Search_engine:
            case 'mascot':
                if file.suffix == '.csv':
                    all_files_to_analyse.append((file, None))
            case 'maxquant':
                if file.suffix == '.txt':
                    for b in maxquant_bulk(file):
                        all_files_to_analyse.append((file, b))
            case 'manual':
                if file.suffix == '.csv':
                    all_files_to_analyse.append((file, None))
            case 'winnow':
                if file.suffix == '.xlsx':
                    all_files_to_analyse.append((file,None))
    match Search_engine:
        case 'peaks':
            for root, _, files in os.walk(location_searchfiles):
                for file in files:
                    if file.endswith('.csv') and 'peptide' in file and 'peptides' not in file:
                        all_files_to_analyse.append((pathlib.Path(root+'/'+file),None))
    # Output data
    general_summary_output_file: dict[str, list[typing.Any]] = {}
    
    for test_file in all_files_to_analyse:
        AA_codes = retreive_AA_codes()
        if fixed_mod != None:
            for AA,f in fixed_mod:
                AA_codes[AA]=AA_codes[AA]+f
         # Set filenames
        file_name = test_file[0].name
        sample_path = test_file[0].stem
        match Search_engine:
            case 'mascot':
                df2, unimod_db, ids,AA_codes = load_files_mascot(path, test_file[0], AA_codes)
            case 'maxquant':
                exp=test_file[1]
                file_name = file_name.replace('.','_')
                file_name = f"{exp}"
                sample_path = f"{exp}"
                df2, unimod_db, ids,AA_codes = load_files_maxquant(path, test_file[0], exp, AA_codes)
            case 'manual':
                df2, unimod_db, ids,AA_codes = load_manual_files(path, test_file[0], AA_codes)
            case 'winnow':
                overall_ptm = {'ox':'Oxidation'}
                df2, unimod_db, ids,AA_codes,overall_ptm = load_files_winnow(path,test_file[0],AA_codes, overall_ptm)
            case 'peaks':
                file_name = str(test_file).replace('\\','/').split('/')[-2]
                file_name = file_name.replace('.','_')
                sample_path = file_name
                df2, unimod_db, ids,AA_codes = load_files_peaks(path, test_file[0], AA_codes)
        outpath = path / 'Output_Classicol' / sample_path
        outpath.mkdir(parents=True, exist_ok=True)
        outpath_mix = path / 'Output_Classicol' / sample_path / 'mixture_plots'
        outpath_mix.mkdir(parents=True, exist_ok=True)
        outpath_mix_aligned = path / 'Output_Classicol' / sample_path / 'mixture_plots' / 'aligned'
        outpath_mix_aligned.mkdir(parents=True, exist_ok=True)
        outpath_mix_fastas = path / 'Output_Classicol' / sample_path / 'mixture_plots' / 'fastas'
        outpath_mix_fastas.mkdir(parents=True, exist_ok=True)
        outpath_mix_taxonomic_output = path / 'Output_Classicol' / sample_path / 'mixture_plots' / 'taxonomic_output'
        outpath_mix_taxonomic_output.mkdir(parents=True, exist_ok=True)
        
        
        # Analyze file and retrieve summary
        file_name = file_name.split('.')[0]
        classicol_done,data_to_remember,taxonomy = ClassiCOL_analysis(
            path,
            sample_path,
            file_name,
            sequence_db,
            df2,
            unimod_db, 
            ids,
            data_to_remember,
            AA_codes,
            lim_tax, 
            demo_tf,
            cpu_count,
        )
        general_summary_output_file[file_name]=classicol_done
        
        files_out = []
        for files in os.walk(path / 'Output_Classicol' / sample_path):
            for i in files[-1]:
                files_out.append(i)
        if any('Taxonomic_results_after_rescoring' in i for i in files_out):
            
            #mixture analysis starts here
            output_sequences,missed,peptide_list = ClassiCOL_mixture_analysis(
                path,
                sample_path,
                file_name,
                sequence_db,
                cpu_count,
                Mixture_analysis,
                taxonomy
                )
            
            if len(output_sequences)>0:
                print('Making fastas from measured data')
                return_measured_seqs(output_sequences,path,sample_path,peptide_list)
            if len(missed)>0:
                return_missed(missed,path,sample_path)
        
        
        
    #create summary output file
    today = str(date.today())
    summ_out = 'Summary_taxonomic_classification_'+today+'.csv'
    with open(path / 'Output_Classicol'/summ_out, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['Experiment','Taxonomic_restriction','Total_peptide_count_in_input_file',
                         'Species(1)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count',
                         'Species(2)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count',
                         'Species(3)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count'])
        for k,v in general_summary_output_file.items():
            writer.writerow([k]+v)
    
    
# try:
#     name_file = taxa_name+'_consensus_search_space_visualization_'+protein_name+'.html'
#     fig.write_html(path/'Output_Classicol'/sample_path/ 'mixture_plots' / 'aligned' /name_file)
# except:
#     fig.write_html(path/'Output_Classicol'/sample_path/'mixture_plots' / 'aligned' / 'multiple_alignment.html')    
    
    
    
    
    
    
    
    
    