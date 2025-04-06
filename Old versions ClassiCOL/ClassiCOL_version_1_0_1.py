#!/usr/bin/env python3
"""
Created on Thu Aug 22 09:57:08 2024

@author: iaengels
"""
#Classicol_version_1_0_0.py
import argparse
import time
import os
import numpy as np
import pandas as pd
import csv
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
import taxoniq
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
import maxquant
import random
from scipy.spatial.distance import braycurtis
from Bio.Align import substitution_matrices
import plotly.express as px
from datetime import date
import warnings
warnings.filterwarnings("ignore")     
warnings.filterwarnings("ignore", category=DeprecationWarning)

def crap_f(path,add_fasta):#done
    print('making database')
    crap = {}
    files_own = []
    for fastafile in os.walk(path+'/BoneDB'):
        for i in fastafile[-1]:
            if i.endswith('.txt') or i.endswith('.fasta') or i.endswith('.fa'):
                files_own.append(path+'/BoneDB/'+i)
    if add_fasta != None:
        for fastafile in os.walk(add_fasta):
            for i in fastafile[-1]:
                if i.endswith('.txt') or i.endswith('.fasta') or i.endswith('.fa'):
                    files_own.append(add_fasta+'/'+i)
    done = ''
    for i in files_own:
        print(i.split('/')[-1])
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

def do_unimod(path,ptms):#done
    raw_unimod = pd.DataFrame()
    ptm_list = []
    for ptm in ptms:
        if len(ptm)>0:
            temp=ptm.split('(')
            if 'N-term' in temp[1].split(')')[0]:
                ptm_list.append((temp[0].replace(' ',''),'N-term'))
            if 'C-term' in temp[1].split(')')[0]:
                ptm_list.append((temp[0].replace(' ',''),'C-term'))
            for a in temp[1].split(')')[0]:
                ptm_list.append((temp[0].replace(' ',''),a))
    with open(path+'/MISC/unimod.txt') as r:
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
        for ptm,aas in ptm_list: #We only check for ptms that mascot searched for
            if ptm==p and aas==aa:
                drop.append(True)
                added=True
                break
        if added==False:
            drop.append(False)
    unimod_db['digit']=np.array(drop)
    unimod_db=unimod_db[unimod_db['digit'] == True]
    return unimod_db, unimod

def find_mass(peptide,AA_codes): #sums the masses of the amino acid sequence inputed #done
    mass = 0
    for AA in peptide:
        mass += AA_codes[AA]
    return mass    

def make_matrix(codes,uni):#done
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

def load_files_mascot(path, name_file):#done
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
    protein = {}
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
    adj_pep=[]
    for i, peptide in enumerate(df2['pep_seq'].values):
        adjusted_pept = ''
        for AA in peptide:
            adjusted_pept += AA+'|'
        adj_pep.append(adjusted_pept[:-1])
    df2['adj']=np.array(adj_pep)
    df2['charge']=charges
    return df, df2, unimod_db, protein, ids, adj_pep

def maxquant_bulk(name_file):#done
    print('Checking for bulk maxquant')
    #name_file = name of the file benchmark
    MQ_file = maxquant.io.read_maxquant(name_file)
    bulk = set(list(MQ_file['Raw file'].values))
    return list(bulk)

def load_files_maxquant(path, name_file,experiment):#done
    print('open maxquant file')
    #name_file = name of the file benchmark
    MQ_file = maxquant.io.read_maxquant(name_file)
    MQ_file = MQ_file[['Sequence','Modified sequence','Raw file','Charge','Modifications']][MQ_file['Raw file']==experiment]
    MQ_file.columns = ['pep_seq','pep_var_mod_pos','pep_scan_title','charge','mods']
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
    
    adj_pep=[]
    for i, peptide in enumerate(df2['pep_seq'].values):
        adjusted_pept = ''
        for AA in peptide:
            adjusted_pept += AA+'|'
        adj_pep.append(adjusted_pept[:-1])
    df2['adj']=np.array(adj_pep)
    return df, df2, unimod_db, protein, ids, adj_pep

def animals_from_db_input(sequence_db, lim_t,demo):#done
    if lim_t == None:
        print('searching against all species in the database')
    else:
        print('Making selection of {} and random other species'.format(lim_t))
        lim_t = lim_t.split('/')
        lim_sp = []
        for x in lim_t:
            lim_sp.append(' '.join(x.split('_')))
        lim_t=lim_sp
    input_animals = ['Pseudomonas aeruginosa','Sus scrofa']
    skip_animals = []
    Class = {}
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
        
        if anim not in input_animals and anim not in skip_animals:
            if lim_t == None:
                input_animals.append(anim)
                continue
            ncbi_animal = anim
            found = False
            while found == False:
                try:
                    taxon = taxoniq.Taxon(scientific_name=ncbi_animal)
                    found = True
                    taxon = [(t.rank.name, t.scientific_name) for t in taxon.lineage]
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
                except:
                    ncbi_animal = ' '.join(ncbi_animal.split(' ')[:-1])
                if len(ncbi_animal)==0:
                    found= True
                    print('Species {} has no taxonomy'.format(anim))
                    skip_animals.append(anim)
    if lim_t != None and demo==False:
        for k,v in Class.items():
            print(k)
            random.shuffle(v)
            random_animals =v[:15]
            input_animals = list(set(input_animals)|set(random_animals))
            skip_animals = list(set(skip_animals)^set(random_animals))
    return list(set(input_animals)), list(set(skip_animals))

def find_mass_matches(sequence, p_mass,mods,pep,unimod_db,AA_codes,uncertain):
    start = 0
    end = len(sequence)
    possible = []
    temp = 1
    unimod_masses = [0]
    # for mass,mod,aa in unimod_db[['mass','PTM','AA']].values:#seach for peptide seq that explain masses of ptms+seq
    #     if mod+'_'+aa in mods.keys() and aa in pep: #only if PTM also in sequence
    #         for i in range(1,mods[mod+'_'+aa]+1):
    #             unimod_masses.append(float(mass)*i)#PTM-> no PTM
    
    unimod_masses = unimod_masses+[num - t for num in unimod_masses for t in unimod_db['mass'].values] #ptms added is lower backbone mass
    unimod_masses = unimod_masses+[num - t for num in unimod_masses for t in unimod_db['mass'].values] #2 ptms added is lower backbone mass 
    
    unimod_masses=set(unimod_masses)
    while start+temp <= end:
        test_seq = sequence[start:start+temp]#stepwise window slide
        if len(test_seq)>len(pep)+1:
            start += 1
            temp -= 1
            continue
        unkown = False
        if str(re.search('['+''.join(uncertain.keys())+']',str(test_seq))) != 'None' or len(set(test_seq)&set(uncertain.keys()))>0:
            keep_seq = test_seq
            other_seq = ''.join([el for el in test_seq if el in uncertain.keys()])
            test_seq = ''.join([el for el in test_seq if el not in uncertain.keys()])
            unkown = True
        test = find_mass(test_seq,AA_codes)
        if unkown == True :
            missing_mass = p_mass-test
            missing_mass = [missing_mass-num for num in unimod_masses]
            checks = other_seq
            other_seq = [uncertain[el] for el in other_seq]
            poss_seqs = []
            if len(other_seq)<2 or checks.count('X')<=1:#only allow 1 X else everything will start to fit and to many possibilities
                for l in range(0,(len(other_seq))):
                    if len(poss_seqs)==0:
                        poss_seqs = [num for num in other_seq[l]]
                    else:
                        poss_seqs = [num+el for num in poss_seqs for el in other_seq[l]]
            added = False
            mass_too_much = False
            for x in poss_seqs:
                if True in [True if num-0.015<=(test-p_mass+find_mass(x,AA_codes))<=num+0.015 else False for num in unimod_masses]:
                    new_seq = ''.join([num if num in AA_codes.keys() else '!' for num in keep_seq])
                    add_seq = ''
                    count_x = 0
                    for m in new_seq:
                        if m=='!':
                            add_seq+=x[count_x]
                            count_x +=1
                        else:
                            add_seq+=m
                    possible.append((add_seq,(start,start+temp)))
                    if added == False:
                        start += 1
                        temp = 1
                    added = True
                elif p_mass>(test+find_mass(x,AA_codes)):
                    mass_too_much = True
            if added == False and mass_too_much == True:
                temp += 1
            elif added == False:
                start += 1
                temp -= 1
                        
        elif True in [True if num-0.015<=test-p_mass<=num+0.015 else False for num in unimod_masses]: #go if same mass or with a new ptm or without an existing one
            possible.append((test_seq,(start,start+temp)))
            start += 1
            temp = 1
        elif test < p_mass:#slide like a caterpilar
            temp += 1 
        else:
            start += 1
            temp -= 1
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

def program(seq_db,seq_real,peptides,mass_matrix): #compare in-silico peptide with the found peptide
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
    if len(adaptation_real)==0 or len(adaptation_db)==0:
        return False, [],[], []
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
    return (s1_align, s2_align)

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

def ptm_mass(ptm,unimod_db): #add this mass to the mass of the normal sequence
    mass = 0
    PTMs = {}
    if len(ptm)>0:
        for i in ptm.split(';'):
            for p, a,m in unimod_db[['PTM','AA','mass']].values:
                if p in i and a in i.split('(')[1]:
                    temp = [x for x in i if x.isdigit()==True]
                    if len(temp)>0:
                        mass += int(temp[0])*float(m)
                        PTMs[p+'_'+a]=int(temp[0])
                        break
                    else:
                        mass += float(m)
                        PTMs[p+'_'+a]=1
                        break       
    return mass, PTMs

def massblast(db,ids,PTMs, peptides,p_adjs,unimod_db,mascot_pos,title,charge,remember_good,remember_bad,mass_matrix,cap,remember_peptides, animal,AA_codes,uncertain):#done
    done= []
    final_result = []
    done_title = []
    found_peptide = []
    for ilocs,p in enumerate(peptides):
        if title[ilocs] in done_title:
            continue
        if p+PTMs[ilocs] in remember_peptides.keys():
            find = p+PTMs[ilocs]
            if animal not in remember_peptides[find]:
                continue #peptide already tested and not found in test animal
        done_title.append(title[ilocs])
        if len(final_result) > 0 and cap==True:
            start_temp = [num for num in final_result if p in num and PTMs[ilocs] in num]
            if len(start_temp)>1:
                start_temp=start_temp[0]
        #we need to do this step to include PSMs of peptides we already allocated
            if len(start_temp)>0 and cap==True: #need all PSMs for downstream filtering with MS2PIP and DeepLC
                for adding in start_temp:
                    final_result.append((p, adding[1],adding[2], adding[3],'extra_PSM',adding[6],title[ilocs],charge[ilocs]))
                continue
        if (p,PTMs[ilocs]) in done or len(p)>=len(db)-2 or len(p)<6: #do not want to do double work
            continue
        
        if cap==False and p in found_peptide: #If you are only interested in the peptide, than you need only not all possibles
            continue
        done.append((p,PTMs[ilocs]))
        ptm_masses, PTMs_insearch = ptm_mass(PTMs[ilocs],unimod_db)
        p_mass= find_mass(p,AA_codes)+ptm_masses
        p_adj=p_adjs[ilocs]
        
        #remake the peptide as it can have other uncertainty
        exact_match_with_uncertainty = ''
        for l in p:
            u = ''
            for k,v in uncertain.items():
                if l in v:
                    u += k
            if len(u)>0:
                u = '['+l+u+']' #need for this because the B,Z,X can be in the database sequence
            else:
                u = l
            exact_match_with_uncertainty += u
        #Add all the exact matches, als oif X,B or Z in sequence. Only 1 uncertain allowed in output
        added_seq=False
        if str(re.search(exact_match_with_uncertainty,str(db))) != 'None' and len([num for num in re.search(exact_match_with_uncertainty,str(db)).group() if num not in AA_codes.keys()])<=1:#match all exact matches en remind the location
            match = re.search(exact_match_with_uncertainty,str(db)).group()
            location = re.finditer(exact_match_with_uncertainty,str(db))
            locs = []
            for i in location:
                locs.append(i.span())
                ptm_loc,extra_mass = find_ptm_location(PTMs[ilocs],p,unimod_db,mascot_pos[ilocs],ids)
            if len([num for num in match if num not in AA_codes.keys()])>0:
                final_result.append((p,p,PTMs[ilocs],locs,'With uncertainty',ptm_loc,title[ilocs],charge[ilocs]))
                found_peptide.append(p)
            else:
                found_peptide.append(p)
                final_result.append((p,match,PTMs[ilocs],locs,'Original',ptm_loc,title[ilocs],charge[ilocs]))
            added_seq=True
        #if isobaric switch was found, no need to check this again for another animal
        if p in remember_good:
            for testing in remember_good[p]:
                if testing[1] in str(db) and testing[-1]==PTMs[ilocs]:
                    match = re.search(testing[1],str(db)).group()
                    location = re.finditer(testing[1],str(db))
                    locs = []
                    for i in location:
                        locs.append(i.span())
                        ptm_loc,extra_mass = find_ptm_location(PTMs[ilocs],p,unimod_db,mascot_pos[ilocs],ids)
                    final_result.append((p,match,testing[2],locs,'isobaric_r',ptm_loc,title[ilocs],charge[ilocs]))
                    added_seq=True
                    found_peptide.append(p)
        #start isoblast
        if added_seq==False: #if no exact match, than start isobaric switches
            final_addition = False
            possible = find_mass_matches(db, p_mass,PTMs_insearch,p,unimod_db,AA_codes,uncertain)
            for test_seq,location in possible:
                if final_addition == True: #1 match/ sequence is enough
                    continue
                seq1 = Seq(p) #original sequence
                seq2 = Seq(test_seq)#database sequence
                alignment = do_alignment(seq1, seq2)
                alteration1 = alignment[1] 
                alteration2 = alignment[0]
                if alteration1.count('-')<=3:# three isobaric mistakes allowed, else too much possibilities, and overall almost never correct #math.ceil(len(p)/5) and alteration2.count('-')<=math.ceil(len(p)/5):#allow switches based in peptide length, longer peptides are allowed to have multiple mistakes
                    test, adapt_real, adapt_db,combo = program(alteration1, alteration2, p_adj,mass_matrix)
                    if test!=False:
                        adapt, addition_ok,ptms_loc = locate_switches(adapt_real,adapt_db,alteration2,alteration1,combo)
                        if addition_ok == True:
                            addition_ok = check_ptms_mascot(adapt,PTMs[ilocs],ids)
                            if addition_ok == True:
                                adapt.append(PTMs[ilocs])
                                ptm_loc,extra_mass = find_ptm_location(PTMs[ilocs],test_seq,unimod_db,mascot_pos[ilocs],ids,ptms_loc)
                                tryptic = False
                                if (p[-1] in ['R','K'] and test_seq[-1] in ['R','K']) or (p[-1] not in ['R','K'] and test_seq[-1] not in ['R','K']) or (p[-1] not in ['R','K'] and test_seq[-1] in ['R','K']):
                                    tryptic = True
                                if ptm_loc != False and abs(p_mass-(find_mass(str(test_seq),AA_codes)+extra_mass))<0.05 and tryptic == True:
                                    final_result.append((p,str(test_seq),adapt,[location],'isobaric',ptm_loc,title[ilocs],charge[ilocs]))
                                    final_addition = True #this could give a problem
                                    found_peptide.append(p)
                                    if p in remember_good:
                                        remember_good[p]=remember_good[p]+[(p,str(test_seq),adapt,[location],'isobaric',ptm_loc,title[ilocs],charge[ilocs],PTMs[ilocs])]
                                    else:
                                        remember_good[p]=[(p,str(test_seq),adapt,[location],'isobaric',ptm_loc,title[ilocs],charge[ilocs],PTMs[ilocs])]
    return final_result, remember_good, remember_bad 

def thread_worker(process_executor,x,df2,rg,rb,unimod_db,mass_matrix,ids,test_animal,cap, remember_peptides,AA_codes,uncertain):#done
    name = x[1]
    sequence=x[0]
    skip_animals = x[2]
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
    if anim in skip_animals or anim != test_animal:
        return 'skip'
    print(name)
    future = process_executor.submit(massblast, sequence,ids,df2['pep_var_mod'].values,df2['pep_seq'].values,adj_pep,unimod_db,df2['pep_var_mod_pos'].values,df2['pep_scan_title'].values, df2['charge'].values,rg,rb,mass_matrix,cap,remember_peptides, anim,AA_codes,uncertain)
    thread_out = future.result() 
    for k,v in thread_out[2].items():
        if k in rb:
            rb[k]=list(set(rb[k]+v))
        else:
            rb[k]=v
    for k,v in thread_out[1].items():
        if k in rg:
            rg[k]=rg[k]+[num for num in v if num not in rg[k]]
        else:
            rg[k]=v
    thread_out = thread_out[0]
    if len(thread_out)>0:
        array = [(anim,name,num[0],num[1],num[2],num[3],num[4],num[5],num[6],num[7]) for num in thread_out]
    else:
        array = [False]
    return array, rg,rb

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

def thread_align(i,seq,m,calculator,sdb):
    seq_keep = sdb[seq]
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

def find_taxons(columns):
    taxons={}
    for i in columns:
        ncbi_animal = i
        
        found = False
        while found == False:
            try:
                taxon = taxoniq.Taxon(scientific_name=ncbi_animal)
                found = True
            except:
                ncbi_animal = ' '.join(ncbi_animal.split(' ')[:-1])
            if len(ncbi_animal)<=1:
                found= True
                ncbi_animal = 'unkown'            
        if ncbi_animal != 'unkown':
            taxon = [(t.rank.name, t.scientific_name) for t in taxon.lineage]
        else:
            taxon=[('species','unkown')]
        taxons[i]=taxon
    return taxons

def make_plots_coverage_per_animal(df, sequence_db, all_animals, data_array,labels,index, locs_dict, columns, z,names, file_name,sample_path):    
    df_plot=pd.DataFrame(np.array(z),index=index,columns=columns)
    print('calculating taxonomy tree')
    taxons = find_taxons(columns)

    taxon_distance=[]
    taxon_col = []
    for i in columns:
        j = taxons[i]
        taxon_col.append(i)
        temp = []
        for t in columns:
            u = taxons[t]
            count = 0
            if i==t or j==u:
                temp.append(count)
                continue
            added = False
            z = set(u)&set(j)
            if len(j)<len(u) and len(j) != 1: #need for minimum, because subspecies and unkown species to ncbi are allocated as 'species' only 
                test_t = j
            else:
                test_t = u
            for y in test_t:
                if y in z:
                    temp.append(count)
                    added = True
                    break
                count += 1
            if added==False:
                temp.append(1)
        taxon_distance.append(temp)
    print('end taxonomy tree')
    taxon_distance = np.array(taxon_distance)
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
    fig.write_html(path+'/Output_Classicol/'+sample_path+'/Heatmap_'+file_name+'.html')
    plotly.offline.plot(fig)
    
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
            # out = knapzak(div)
            # if len(out) <= 1:
                #cannot take out any clusters 
            end = True
            
            results.append((list(div[3]),list(div[1]),False))
            # else:
            #     del_next_it = []
            #     del_next_it_proteins = []
            #     done = []
            #     for i in out:
            #         if i in done:
            #             continue
            #         done.append(out)
            #         t_i = div[4].loc[i]
            #         t_i = t_i.loc[:, (t_i != 0).any(axis=0)]
            #         t_i = list(t_i.columns)
            #         i = [i,t_i]
            #         results.append(i)#add cluster to outcome
            #         del_next_it = del_next_it+i[0]
            #         df_distance = df_distance.drop(i[1], axis = 0)
            #         df_distance = df_distance.drop(i[1], axis = 1)
            #         del_next_it_proteins = del_next_it_proteins + i[1]
            #         df_distance_taxon = df_distance_taxon.drop(i[0],axis=1)
            #         df_distance_taxon = df_distance_taxon.drop(i[0],axis=0)
            #     #prep for next iteration
            #     new_h = div[4]
            #     new_h = new_h[~new_h.index.isin(del_next_it)]
            #     new_h = new_h.drop(del_next_it_proteins, axis=1)
            #     t = []
            #     for x in new_h.columns:
            #         if np.sum(new_h[x]) != 0:
            #             t.append(x)
            #     new_X = df_distance[df_distance.index.isin(t)]
            #     new_X = new_X[new_h.columns]
            #     new_names = new_X.columns
            #     new_X = new_X[t].values
                
            #     new_Y = df_distance_taxon[~df_distance_taxon.index.isin(del_next_it)]

            #     new_animal = new_Y.index
            #     new_Y = new_Y[new_animal].values
            #     divide.append([new_X,new_names,new_Y,new_animal,new_h])
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
            print('reiteration ...')
            for i in score:
                contain_i = set(dfs['found_match'][dfs['animal'].isin(i)].values)
                re = list(set(dfs['found_match'][dfs['animal'].isin(i)].values))
                keep[tuple(i)]=True
                for t in score:
                    if t==i:
                        continue
                    contain_t = set(dfs['found_match'][dfs['animal'].isin(t)].values)
                    if contain_i.issubset(contain_t)==True and len(contain_i^contain_t)!=0:
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
                            if contain_x.issubset(contain_y)==True and len(contain_x^contain_y)!=0:
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
                    print('deleting {} because lower than 10% of total peptides'.format(i))
                    delete.append(i)
            score = [num for num in score if num not in delete]
            
            combos = []
            for el in score:
                for num in score:
                    if el != num and tuple(sorted(list(el+num))) not in combos:
                        combos.append(tuple(sorted(list(el+num))))
            print('Catching sneaky ones ...')
            output_score = []
            keep = {}
            deleted = []
            for i in score:
                keep[tuple(i)]=True
                contain_i = set(dfs['found_match'][dfs['animal'].isin(i)].values)
                for x in combos:
                    if len([num for num in deleted if len(set(list(num))&set(list(x)))>0]):
                        continue
                    if len([num for num in i if num not in list(x)])==0:
                        continue
                    contain_x = set(dfs['found_match'][dfs['animal'].isin(x)].values)
                    if (contain_i.issubset(contain_x) and (len(contain_i)<len(contain_x)*0.85 or sorted(list(contain_x))==sorted(list(contain_i)))):
                        delete = True
                        deleted.append(i)
                        for animal_check in x:
                            check_temp = set(dfs['found_match'][dfs['animal'].isin([animal_check])].values)
                            if check_temp.issubset(contain_i) or len(contain_i^check_temp)<=1:
                                delete = False
                        if delete==True:
                            keep[tuple(i)]=False
                            print('deleting',i)
                            break
                output_score=[el for el in keep.keys() if keep[el]!=False]
            if len(output_score)>0:
                output2.append(output_score)
            else:
                output2.append(score)
    return output2, recall


def make_sunburst(dfs,all_animals,df_output,file_name,ip_animals,df_plot,path,sample_path):
    #sunburst plot to retrace the track 
    all_taxonomy = find_taxons(list(set(dfs['animal'].values)))
    
    labels = []
    values = []
    parents = []
    done = []
    values_iso = []
    values_uni = []
    braycurtis_dist = []
    all_peptides = list(set(list(dfs['mascot_peptide'].values)))
    array1 = [1 if pep in all_peptides else 0 for pep in all_peptides]
    seq_concat = '_'.join(all_peptides)
    weights_bc = [1/seq_concat.count(val) for val in all_peptides]
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
        combo_taxon = []
        for key,val in all_taxonomy.items():
            if key != animal:
                combo_taxon = combo_taxon+val
        combo_taxon=set(combo_taxon)
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
            if t in combo_taxon or t[1] in parents:
                
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
                temp_peptides = set(dfs['mascot_peptide'][dfs['animal'].isin(animals_involved)].values)
                #we want only the shared peptides/level. At species level the uniqueness does count
                unique_peptides = []
                for peptide_test in temp_peptides:
                    test_pep = set(dfs['animal'][(dfs['found_match']==peptide_test)&(dfs['animal'].isin(animals_involved))].values)
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
        if (x[1],x[2]) not in done and x[2] not in labels: #if species and sub-species in database than it can be that it is in there double
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
    hovertemplate='<b>%{label} </b> <br> Taxon score: %{color:.3f}<br> Total peptide count: %{value:.0f} <br> isoBLAST peptides: %{meta:.0f}'))
    fig.update_layout(title = 'Score of output species to sample '+file_name)
    fig.update_layout(autosize=True,margin = dict(t=30, l=0, r=0, b=0))
    fig.update_layout(coloraxis_colorbar_title='Score')
    try:
        fig.write_html(path+'/Output_Classicol/'+sample_path+'/sunburst_'+file_name+'.html')
    except:
        fig.write_html(path+'/Output_Classicol/'+sample_path+'/sunburst.html')
    plotly.offline.plot(fig)
    time.sleep(2)
    f_out = make_output_file(path, df_plot, df_output,file_name, braycurtis_dist, labels,sample_path)
    make_sunburst_with_missing2(dfs,ip_animals,list(set(df_output['animal'].values)),path,sample_path)
    return braycurtis_dist, labels, f_out

def make_output_file(path,df_plot, df_output, file_name, bc, l,sample_path):
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
    with open(path+'/Output_Classicol/'+sample_path+'/ZooMS_results_'+file_name+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['ZooMS analysis output (~: isoBLAST match, *: Unique PSM)'])
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
    return path+'/Output_Classicol/'+sample_path+'/ZooMS_results_'+file_name+'.csv'

def make_connection_graph(df_output, file_name,final_output,found_animals):
    taxon = []
    for i in df_output['animal'].values:
        taxon.append(final_output[i])
    df_output['taxon']=taxon

    df_plot = pd.DataFrame(columns = ['taxon','found_match','mascot_peptide','protein'])
    temp = [tuple(el) for el in df_output[['taxon','found_match','mascot_peptide','protein']].values]
    for i in set(temp):
        df_temp = pd.DataFrame(np.array(i).reshape(1,-1),columns = ['taxon','found_match','mascot_peptide','protein'])
        df_plot = pd.concat([df_plot,df_temp],ignore_index=True)

    pro = []
    added = []
    for i in df_plot['protein'].values:
        temp = i.split('[')[0]
        if '=' in temp:
            temp = i.split('=')[1]
            pro.append('_'.join(temp.split(' ')[0]))
            added.append('_'.join(temp.split(' ')[0]))
        else:
            pro.append('_'.join(temp.split(' ')[1:]))
            added.append('_'.join(temp.split(' ')[1:]))
    df_plot['protein']=pro
    
    unique = []
    for i in df_plot['mascot_peptide'].values:
        temp = sorted(list(set(df_plot['taxon'][df_plot['mascot_peptide']==i])))
        if len(temp)==1:
            unique.append('U_'+' + '.join(temp))
        else:
            all_taxonomy = find_taxons(found_animals)
            if temp[0] in all_taxonomy:
                start = all_taxonomy[temp[0]]
                for x in all_taxonomy:
                    start = set(all_taxonomy[x])&set(start)
                for t in all_taxonomy[temp[0]]:
                    if t in start:
                        unique.append('Shared_'+t[1])
                        break
            else:
                unique.append('other')
    df_plot['uniqueness']=unique
    return df_plot, df_output
 
def rescore(path, file_loc,f_name,sample_path):
    file_name = file_loc
    file_name_print = file_name.split('/')[-1]
    file_name_print = file_name_print.split('.csv')[0]
    df = pd.read_csv(file_name,names=['a','b','c'])
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
            i = i.split('OS=')[-1]
            i = i.split(' OX')[0]
            score.append(s)
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
    
    if max(list(df['score'].values))>0.25:
        df = df[df['score']>0.25]
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
            fig = px.bar(temp, x="peptides", y='count', color="species", title=file_name_print+": Uncommon peptide to species cluster "+str(k))
            fig.update_layout(height=600, width=1500, barmode = 'stack', xaxis={'categoryorder':'total descending'})
            fig.update_xaxes(showticklabels=False)
            fig.write_html(path+'/Output_Classicol/'+sample_path+'/Barplot_uniquePeptides_'+file_name_print+'.html')
            plotly.offline.plot(fig)
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
                fig.write_html(path+'/Output_Classicol/'+sample_path+'/Rescore_'+file_name_print+'.html')
            except:
                fig.write_html(path+'/Output_Classicol/'+sample_path+'/Rescore.html')
            plotly.offline.plot(fig)
            time.sleep(3)
        else:
            df_scores = pd.DataFrame()
            df_scores['species']=['']
            df_scores['score original']=['']
            df_scores['Rescore']=['']
            fig = px.bar(df, x=['all peptides'], y=[1],color=species, title=file_name_print+": Only 1 species in cluster "+str(k))
            try:
                fig.write_html(path+'/Output_Classicol/'+sample_path+'/Barplot_uniquePeptides_'+file_name_print+'.html')
            except:
                fig.write_html(path+'/Output_Classicol/'+sample_path+'/Barplot_uniquePeptides.html')
            plotly.offline.plot(fig)
            time.sleep(3)
    to_summ = make_output_file_after_rescoring(path,df,df_scores, f_name,sample_path)
    return to_summ

def go_to_species(c,species_list):
    for x in c:
        if x=='environmental samples' or 'sp.' in x:
            continue
        t = taxoniq.Taxon(scientific_name=x)
        t_lin = [(tax.rank.name, tax.scientific_name) for tax in t.lineage]
        if ('kingdom', 'Metazoa') not in t_lin:
            continue
        if 'species' in t_lin[0][0]:
            species_list.append(x)
            try:
                gtsp = [el.scientific_name for el in t.child_nodes if 'sp.' not in el.scientific_name and 'unclassified' not in el.scientific_name and el.scientific_name.isdigit()==False]
                species_list = species_list +gtsp
            except:
                continue
        else:
            try:
                gtsp = [el.scientific_name for el in t.child_nodes if 'sp.' not in el.scientific_name and 'unclassified' not in el.scientific_name and el.scientific_name.isdigit()==False]
                gts = []
                for g in gtsp:
                    g_no_numbers = ''.join([num for num in g if num.isdigit()==False])
                    if g_no_numbers==g:
                        gts.append(g)
                species_list = species_list + go_to_species(gts,[])
            except:
                continue
    return species_list

def make_output_file_after_rescoring(path,df_og,df_rescore, file_name,sample_path):
    print('generating final output file')
    with open(path+'/Output_Classicol/'+sample_path+'/Taxonomic_results_after_rescoring_'+file_name+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['Species','Species_Rank','Original_score','Rescored','Unique_peptide_amongst_candidates',
                         'Ambiguous_Peptide','Protein','Peptides','PTMs','Representative_PSM_title'])
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
            if ptm=='Nan' or ptm=='nan':
                ptm = ''
            writer.writerow([sp,gr,sc,rescore,u,amb,pro,p,ptm,t])
    summary=[]
    for gr in range(1,4):
        temp = df_og[df_og['group']==gr]
        if len(temp)==0:
            summary = summary+['']*6
            continue
        temp_sp=sorted(list(set(temp['species'].values)))
        taxons=find_taxons(temp_sp)
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

def background_check(c,ip):
    out = {}
    for q in c:
        m = go_to_species([q],[])
        if len(set(m)&set(ip))==0:
            out[q] = m
    return out

def make_sunburst_with_missing2(dfs,ip_animals,found_animals,path,sample_path):
    print('start making sunburst with missingness')
    df_missing = dfs[['mascot_peptide','found_match','animal','type']]
    #ip_animals are all animals that were searched against
    additional = {}
    additional_animals = []
    for i in found_animals:#for all the species we found from database
        print('Looking for related species of {}'.format(i))
        try:
            taxon = taxoniq.Taxon(scientific_name=i)
        except:
            continue
        t_lin = [(tax.rank.name, tax.scientific_name) for tax in taxon.lineage]#get taxonomic lineage
        #find children nodes for each considered lineage uptil the order
        for find_child in t_lin:
            if any(el in find_child[0] for el in ['species','genus','family','order']):#subspecies also need to be plotted if any
                f_child = find_child[1]
                try:
                    taxon = taxoniq.Taxon(scientific_name=f_child)
                    children = taxon.child_nodes
                    children = [el.scientific_name for el in children if 'sp.' not in el.scientific_name and el.scientific_name.isdigit()==False]
                    missing = background_check(children,ip_animals)
                    for k,aa in missing.items():
                        # additional_animals=additional_animals + aa
                        if 'unclassified' not in k:#len(aa)==0 and 
                            additional_animals.append(k)
                    additional[f_child]=list(missing.keys())#only keep 1 level lower
                except:
                    continue
    
    all_missing_animals = found_animals+additional_animals
    for m in set(additional_animals):
        df_missing_add = np.array(['','',m,'Original'])
        df_missing_add = pd.DataFrame(df_missing_add.reshape(1,-1),columns=df_missing.columns)
        df_missing = pd.concat([df_missing,df_missing_add],ignore_index=True)
    all_taxonomy = find_taxons(all_missing_animals)
    print('done looking')
    labels = []
    values = []
    parents = []
    done = []
    values_iso = []
    braycurtis_dist = []
    
    all_peptides = list(set(list(df_missing['mascot_peptide'].values)))
    array1 = [1 if pep in all_peptides else 0 for pep in all_peptides]
    seq_concat = '_'.join(all_peptides)
    weights_bc = [1/seq_concat.count(val) for val in all_peptides]          
    print('start assigning branch scores missing')
    for animal in all_missing_animals:#first the found animals are mapped, than the missing species are mapped to the existing tree
        
        a_taxon = all_taxonomy[animal]
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
            temp_peptides = set(df_missing['mascot_peptide'][df_missing['animal']==animal].values)
            array2 = [1 if pep in temp_peptides else 0 for pep in all_peptides]
            score = braycurtis(array1, array2,w=weights_bc)
            braycurtis_dist.append(1-score)#columns need to be the same
        combo_taxon = []
        for key,val in all_taxonomy.items():
            if key != animal:
                combo_taxon = combo_taxon+val
        combo_taxon=set(combo_taxon)
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
            if t in combo_taxon or t[1] in parents:
                
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
                    else:
                        continue
                if [t[1]] in done and 'species' not in t:
                    parents.append(t[1])
                    
                    break
                done.append([t[1]])
                parents.append(t[1])
                labels.append(t[1])
                values.append(add) 
                values_iso.append(add_iso)
                temp_peptides = set(df_missing['mascot_peptide'][df_missing['animal'].isin(animals_involved)].values)
                #we want only the shared peptides/level. At species level the uniqueness does count
                unique_peptides = []
                animals_involved_og = [key for key,val in all_taxonomy.items() if t in val and key in dfs['animal'].values]
                for peptide_test in temp_peptides:
                    test_pep = set(dfs['animal'][(dfs['found_match']==peptide_test)&(dfs['animal'].isin(animals_involved_og))].values)
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
        if (x[1],x[2]) not in done:
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
    
    hovertemplate='<b>%{label} </b> <br> Taxon score: %{color:.3f}<br> Total peptide count: %{value:.0f} <br> isoBLAST peptides: %{meta:.0f}'))
    fig.update_layout(title = 'Score of output species to sample with species not in DB: '+file_name)
    fig.update_layout(autosize=True,margin = dict(t=30, l=0, r=0, b=0))
    fig.update_layout(coloraxis_colorbar_title='Score')
    try:
        fig.write_html(path+'/Output_Classicol/'+sample_path+'/sunburst_including_missing_species_'+file_name+'.html')
    except:
        fig.write_html(path+'/Output_Classicol/'+sample_path+'/sunburst_including_missing_species.html')
    plotly.offline.plot(fig)
    time.sleep(2)
    return temporary_dataframe
    
def load_manual_files(path,test_file):
    #peprec file columns = ['peptide','ptm'] ptm= 0|N-term|1|Oxidation|-1|C-term
    print('open file {}'.format(test_file))
    df = pd.read_csv(test_file, header=0)
    df = df.fillna('')
    df.columns = ['pep_seq','pep_var_mod']
    charges = [2]*len(df)#add charges for completeness, could be deleted later
    file_name = test_file.split('/')[-1]
    file_name = file_name.split('.')[0]
    df['pep_scan_title']=[file_name+'_'+str(num) for num in range(0,len(df))]#add_unique_title
    df['prot_desc']=['unkown_protein']*len(df)
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
    protein = {}
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
    adj_pep=[]
    for i, peptide in enumerate(df2['pep_seq'].values):
        adjusted_pept = ''
        for AA in peptide:
            adjusted_pept += AA+'|'
        adj_pep.append(adjusted_pept[:-1])
    df2['adj']=np.array(adj_pep)
    df2['charge']=charges
    return df, df2, unimod_db, protein, ids, adj_pep

if __name__ == "__main__":#rescore output in csv, reduced info output file, summary output
    AA_codes = {"A" : 71.03711, "C" : 103.00919, "D" : 115.02694, "E" : 129.04259, "F" : 147.06841,
                "G" : 57.02146, "H" : 137.05891, "L" : 113.08406, "M" : 131.04049, "K" : 128.09496,
                "N" : 114.04293,"P" : 97.05276, "Q" : 128.05858, "S" : 87.03203,"R" : 156.10111, 
                "T" : 101.04768, "V" : 99.06841, "W" : 186.07931, "Y" : 163.06333,"I" : 113.08406,'U':168.05,'&':1000}

    uncertain = {'B':['D','N'],'Z':['Q','E'],'X':['A','C','D','E','F','G','H','I','L','K','M','N','P','Q','R','S','T','V','W']}

    ##############PROMPT INPUT ARGUMENTS###############################
    parser = argparse.ArgumentParser(# -e ERROR \n   -p PEPTIDE_TABLE [PEPTIDE_TABLE] |    [-l LIMIT]\n   [-t TAXONOMY]\n   [-n NEIGHBOURING] [-a]
                                         usage="py ClassiCOL_version_1_0_0.py \n -f Path to additional Database in FASTA format (not required)\n -d DIRECTORY to classicol \n -l folderNAME with search engine output files\n -s Search engine (MASCOT or MaxQuant)\n -t limitation of taxonomy (e.g. Pecora or for species: Bos_taurus or Homo_sapiens/Canis)\n -c CPUs maximum to be used",
                                         description="ClassiCOL species classification via peptide ambiguation")


    inputs = parser.add_argument_group('\nVariable inputs')
    inputs.add_argument("-f", dest="Manual_fasta", help="path to fasta document other than standard classicol database (ends in .fasta or .txt)", type=str)
    inputs.add_argument("-d", dest="Directory",  help="Directory where Classicol is located", type=str, required=True)
    inputs.add_argument("-s", dest="Search_engine", help="Search engine used, allowed types are: MASCOT .cvs and MaxQuant .txt", type=str, required=True)
    inputs.add_argument("-t",dest="limited_taxonomy", help="Subset of taxonomy e.g. Pecora or Pecora/Primates", type=str)
    inputs.add_argument("-l",dest="File_location", help="path to folder containing search engine files", type=str)
    inputs.add_argument("-m", dest="Fixed_modification", help="Fixed modification e.g. C,45.98/M,...", type=str)
    # inputs.add_argument("-v", dest="Variable_modification", help="Variable modification e.g. C,45.98/M,...", type=str)
    inputs.add_argument("-c", dest="CPUs_to_be_used", help="Amount of CPUs the algorithm can maximally use", type=str, required=True)
    args = parser.parse_args()
    add_fasta, path, Search_engine, lim_tax, location_searchfiles, fixed_mod,CPUs_to_be_used = args.Manual_fasta, args.Directory, args.Search_engine, args.limited_taxonomy, args.File_location, args.Fixed_modification, args.CPUs_to_be_used
    if location_searchfiles.lower() == 'demo':
        location_searchfiles = path+'/Demo'
        demo_tf = True
        if lim_tax==None:
            lim_tax = 'Pecora'
    else:
        demo_tf = False
        
    if fixed_mod != None:
        fixed_mod = fixed_mod.split('|')
        fixed_mod = [num.split(',') for num in fixed_mod]
        fixed_mod = [(el[0],float(el[1])) for el in fixed_mod]
    ##
    CPUs_to_be_used = int(CPUs_to_be_used)
    CPUs = min(CPUs_to_be_used, multiprocessing.cpu_count()-3)
    Search_engine = Search_engine.lower()
    outpath = path+'/Output_Classicol'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    os.chdir(path)
    all_files_to_analyse = []
    sequence_db = crap_f(path, add_fasta)
    AA_codes_start = AA_codes.copy()
    remember_peptides = {} #remember for batch search
    rg= {} #remember for within search
    ############################################################################
    data_to_remember = pd.DataFrame(columns=['animal','protein','mascot_peptide',
                                             'found_match','switch','location','type','PTM','PTM_loc'])
    ############################################################################
    for mascotfile in os.walk(location_searchfiles):
        for i in mascotfile[-1]:
            if i.endswith('.csv') and Search_engine == 'mascot':
                all_files_to_analyse.append(location_searchfiles+'/'+i)
            ###################################################
            elif i.endswith('.txt') and Search_engine == 'maxquant':
                bulk = maxquant_bulk(location_searchfiles+'/'+i)
                for b in bulk:
                    all_files_to_analyse.append([location_searchfiles+'/'+i,b])
            ###################################################
            elif i.endswith('.csv') and Search_engine == 'manual':
                all_files_to_analyse.append(location_searchfiles+'/'+i)
            ###################################################
    general_summary_output_file = {}
    for test_file in all_files_to_analyse:
        
        AA_codes = AA_codes_start.copy()
        if fixed_mod != None:
            for AA,f in fixed_mod:
                AA_codes[AA]=AA_codes[AA]+f

        if Search_engine == 'mascot':
            df, df2, unimod_db, protein, ids, adj_pep = load_files_mascot(path,test_file)
        #######################################################
        elif Search_engine == 'maxquant':
            exp=test_file[1]
            test_file = test_file[0]
            df, df2, unimod_db, protein, ids, adj_pep = load_files_maxquant(path,test_file,exp)
        ######################################################
        elif Search_engine == 'manual':
            df, df2, unimod_db, protein, ids, adj_pep = load_manual_files(path,test_file)
        else:
            print('Invalid search engine')
            break
        file_name = test_file.split('/')[-1]
        file_name = file_name.split('.')[0]
        if Search_engine == 'maxquant':
            file_name = file_name+'_'+exp
        sample_path=file_name
        outpath = path+'/Output_Classicol/'+sample_path
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        ######################################################
        print('getting rid of keratins, Trypsin, Lys-C')
        contamination = []
        for i in df2['prot_desc'].values:
            if 'keratin' in i.lower():
                contamination.append(True)
            else:
                contamination.append(False)
        df2['contamination']=contamination
        df2 = df2[df2['contamination']==False]
        
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
        
        ###########################################################################
        print('Adjusting data for batch search')
        
        if len(data_to_remember) == 0:
            already = [False]*len(df2)
        else:
            data_to_remember['switch']=[', '.join(num) for num in data_to_remember['switch'].values]
            data_to_remember['location']=[str(num) for num in data_to_remember['location'].values]
            data_to_remember = data_to_remember.drop_duplicates()
            data_to_remember['location']=[eval(num) for num in data_to_remember['location'].values]
            already = []
            for k,v in df2[['pep_seq','pep_var_mod_pos']].values:
                comb = [k,''.join([el if (el == '0' or el=='.') else '1' for el in v])]
                found = True
                for dtm in data_to_remember[['mascot_peptide','PTM_loc']].values:
                    dtm = [dtm[0],''.join([el if (el == '0' or el=='.') else '1' for el in dtm[1]])]
                    if dtm == comb:
                        already.append(True)
                        found = True
                        break
                    else:
                        
                        found = False
                if found == False:
                    already.append(False)
        
        df2['already']=already#THIS WAS ADDED DURING DEBUGGING
        df_already_found = df2[df2['already']==True]
        print('Reduced the dataframe with {} peptides because they were already matched'.format(len(df_already_found)))
        df2 = df2[df2['already']==False]
        ###########################################################################
        
        print('making matrix')
        mass_matrix = make_matrix(AA_codes,unimod_db)
        print('Binary matrix creation successful')
        
        print('starting isoBLAST search')
        end_result = {}
        columns=['animal','protein','mascot_peptide','found_match','switch','location','type','PTM','title','charge']
        df_out = pd.DataFrame(columns=columns)
        count = 0
        input_animals,skip_animals = animals_from_db_input(sequence_db, lim_tax,demo_tf)
        input_animals = sorted(input_animals)
        rb = {}
        
        print('using {} CPUs'.format(CPUs))
        total_peptide_begin = len(df2)
        check_all_psms = False
        with ProcessPoolExecutor() as process_executor:
            for a in input_animals:
                count += 1
                print('{}/{} starting on {} isoBLASTing {} different peptides'.format(count, len(input_animals),a,len(df2)))
                start_time = time.time()
                with ThreadPoolExecutor(CPUs) as thread_executor:
                    futures = [thread_executor.submit(thread_worker, process_executor, (sequence,name,skip_animals),df2,rg,rb,unimod_db,mass_matrix,ids,a,check_all_psms, remember_peptides,AA_codes,uncertain) for sequence,name in sequence_db.items()]
                    results = [future.result() for future in as_completed(futures)]
                    print('Took {} minutes'.format((time.time()-start_time)/60))
                    for output in results:
                        if 'skip' in output:
                            continue
                        if False in output and len(output)==1:
                            continue
                        for k,v in output[2].items():
                            if k in rb:
                                rb[k]=list(set(rb[k]+v))
                            else:
                                rb[k]=v
                        for k,v in output[1].items():
                            if k in rg:
                                rg[k]=rg[k]+[num for num in v if num not in rg[k]]
                            else:
                                rg[k]=v
                        if False in output[0]:
                            continue
                        df_add = pd.DataFrame(data=output[0],columns=columns)
                        df_out=pd.concat([df_out,df_add], ignore_index = True)
        #######################################################
        print('adding position information')
        temp = {(a1,a2):a3 for a1,a2,a3 in df2[['pep_seq','pep_scan_title','pep_var_mod_pos']].values}
        add = [temp[(a1,a2)] for a1,a2 in df_out[['mascot_peptide','title']].values]
        df_out['PTM_loc']=add
        #######################################################
        print('Remove spectra isoBLAST link to trypsin and Lyc-C')
        linked_to_contaminants = []
        for p,t in df_out[['protein','title']].values:
            if 'Pseudomonas' in p or 'trypsin' in p.lower():
                linked_to_contaminants.append(t)
        dfs = df_out[df_out['title'].isin(linked_to_contaminants) == False]
        ######################################################################
        title_already_found = df_already_found['pep_scan_title'].values
        retreived = pd.DataFrame(columns=dfs.columns)
        for t in title_already_found:
            temp = df_already_found[df_already_found['pep_scan_title']==t]
            for k,v in temp[['pep_seq','pep_var_mod_pos']].values:
                temp_found = data_to_remember[(data_to_remember['mascot_peptide']==k) &
                                              (data_to_remember['PTM_loc']==v)]
                temp_found['title']=[t]*len(temp_found)
                retreived = pd.concat([retreived,temp_found],ignore_index=True)
                break
        print('Adding {} peptides'.format(len(retreived)))
        dfs = pd.concat([dfs,retreived],ignore_index=True)
        ######################################################################
        print('reducing the data to 1 peptide sequence/peptide')
        drop = []
        added = []
        for p,fm in dfs[['protein','found_match']].values:
            if (p,fm) in added:
                drop.append(True)
            else:
                drop.append(False)
                added.append((p,fm))
        dfs['double']=drop
        dfs = dfs[dfs['double']==False]
        
        print('Remove single hit wonders')#only keep proteins with >1 peptide evidence in other words we do not trust single hit wonders!
        keep = []
        for pro in set(dfs['protein'].values):
            temp = dfs['found_match'][dfs['protein']==pro].values
            a = set(dfs['animal'][dfs['protein']==pro].values)
            if len(set(temp))>1: #>2 peptides per proteins
                keep.append(pro)
        dfs = dfs[dfs['protein'].isin(keep)]
        index, locs_dict, columns, z = calc_distance_collagen(sequence_db, dfs)
        
        #separate the uncertain ones from the analysis
        print('removing uncertain peptides that are single hit wonders')
        nots = list(dfs['found_match'][dfs['type']!='With uncertainty'].values)
        real_uncertain = [el for el in dfs['found_match'][dfs['type']=='With uncertainty'].values if el not in nots]
        keep_uncertain = dfs[(dfs['type']=='With uncertainty') & (dfs['found_match'].isin(real_uncertain))]
        dfs = dfs[dfs['found_match'].isin(real_uncertain)==False]
        if len(set(dfs['protein'].values))>1:#need at least 1 protein for continued search
            print('starting on the heatmap and collagen sequence multiple alignment')
            if not os.path.exists(path+'/MISC/collagen_distance.csv'):
                print('No distance file found ... \n creating csv file')
                if not os.path.exists(path+'/MISC'):
                    os.makedirs(path+'/MISC')
                with open(path+'/MISC/collagen_distance.csv', 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                    writer.writerow(['index','collagen_seq1','collagen_seq2','distance'])
                    writer.writerow(['test',0,0,0])
                    writer.writerow(['test2',0,0,0])
                    writer.writerow(['test3',0,0,0])
                
            ##########################################################
            Z_distance_csv = pd.read_csv(path+'/MISC/collagen_distance.csv', header=0)      
            Z_distance_csv.index = Z_distance_csv['index'].values
            Z_distance_csv = Z_distance_csv[list(Z_distance_csv.columns)[1:]]
            
            #save each output so this step can go faster in the future
            matrix = substitution_matrices.load('BLOSUM90')
            print('calculating distances with {} CPUs'.format(CPUs))
            
            names = index
            temp_db = {val:key for key,val in sequence_db.items()}
            all_seqs= [temp_db[key] for key in names]
            seq_to_name = {v:k for k,v in temp_db.items()}
            calculator = DistanceCalculator('blosum90')
            done = []
            X = []
            # with ProcessPoolExecutor() as process_executor:
            for i in all_seqs:
                print('Starting on alignment of {}'.format(seq_to_name[i]))
                temp = []
                print(all_seqs.index(i)+1,'/',len(all_seqs))
                # with ThreadPoolExecutor(CPUs) as thread_executor:
                align_seqs = []
                for t in all_seqs:
                    if t not in done and (seq_to_name[i] not in Z_distance_csv.columns
                                          or seq_to_name[t] not in Z_distance_csv.columns):#if not considered this run and not in previous runs
                        align_seqs.append(t)
                    elif seq_to_name[i] in Z_distance_csv.columns and seq_to_name[t] in Z_distance_csv.columns:#if considered previously
                        if list(Z_distance_csv[seq_to_name[i]].values)[list(Z_distance_csv.columns).index(seq_to_name[t])]=='nothing':#if it was not matched previously
                            if t not in done:#if it has not been calculated this run already
                                align_seqs.append(t)
                print(len(align_seqs), 'to do')
                if len(align_seqs)>0:
                    results = [thread_align(i,t,matrix,calculator,sequence_db) for t in set(align_seqs)]
                    print('collected')
                    results = {key:val for val,key in results}
                    
                print('Aligned!')
                for t in all_seqs:
                    if t in done:
                        temp.append(float(X[all_seqs.index(t)][all_seqs.index(i)]))
                    elif seq_to_name[i] == seq_to_name[t]:
                        temp.append(0)
                    elif t not in align_seqs and i not in align_seqs:
                        if list(Z_distance_csv[seq_to_name[i]].values)[list(Z_distance_csv.columns).index(seq_to_name[t])]!='nothing':
                            temp.append(float(list(Z_distance_csv[seq_to_name[i]].values)[list(Z_distance_csv.columns).index(seq_to_name[t])]))
                        else:
                            temp.append(float(results[sequence_db[t]]))
                    else:
                        temp.append(float(results[sequence_db[t]]))
                done.append(i)
                X.append(temp)
                print('done')
        
            distance=np.array(X)
            labels = names
            print('end calculating distances')
            
            ########################################################
            all_animals = {el:el for el in input_animals}
            taxon_distance, df_heatmap_values, heat_data,names,animals = make_plots_coverage_per_animal(dfs, sequence_db, all_animals, distance, labels, index, locs_dict, columns, z, names, file_name,sample_path)
            print('find outliers')
            output = find_animals(distance,names,taxon_distance,animals,df_heatmap_values)#eiwit niveau, adapt to filter out bullshit
            #find locations in tree based on protein
            #####################################################
            print('Saving distances ...')
            distance =pd.DataFrame(np.array(distance),columns=names,index=names)
            print('Start retrieving calculations ...')
            in_dist_csv = list(Z_distance_csv.columns)
            print('Adding new columns')
            new_to_add=[]
            for col in distance.columns:
                if col not in Z_distance_csv.columns:
                    temp = []
                    new_to_add.append(col)
                    for col2 in in_dist_csv:
                        if col2 in names:
                            temp.append(float(distance[col][col2]))#freshly calculated
                        else:
                            temp.append('nothing')#not calculated this time
                    #adding new calculations to dataframe
                    Z_distance_csv[col]=temp
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
                if (nums/len(save_cols))%10==0:
                    print('{} of {} done'.format(nums,len(save_cols)))
                if col not in distance.columns or col in new_to_add:#already added or not calculated this round
                    continue
                temp = Z_distance_csv[col].values
                if 'nothing' in temp:
                    change = []
                    for num,el in enumerate(temp):
                        if el == 'nothing':
                            if col in names and save_cols[num] in names:
                                change.append(float(distance[col][save_cols[num]]))#freshly calculated
                                continue
                        change.append(el)
                    Z_distance_csv[col]=change
        
            Z_distance_csv = Z_distance_csv.drop(['collagen_seq1','collagen_seq2','distance'],axis=0,errors='ignore')
            Z_distance_csv = Z_distance_csv.drop(['collagen_seq1','collagen_seq2','distance'],axis=1,errors='ignore')
            print('Save as csv file')
            with open(path+'/MISC/collagen_distance.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                writer.writerow(['index']+list(Z_distance_csv.columns))
                for q in Z_distance_csv.columns:
                    writer.writerow([q]+list(Z_distance_csv[q].values))
            #######################################################################################
            print('Done saving')
            taxon_distance = pd.DataFrame(np.array(taxon_distance),columns=animals,index=animals)
            print('Starting on the walk down the taxonomic tree')
            
            output_final, recall = find_animals2_1(distance,names,taxon_distance,animals,df_heatmap_values,output,dfs)#find animal at protein level
            
            print('End of the walk')
            #find uniqueness for each location in tree
            df_output = pd.DataFrame(columns=dfs.columns)
            
            #rescore the isoblasts to the df_output
            df_isobaric = df_output[df_output['type'].str.contains('isobaric')]
            df_other = df_output[df_output['type'].str.contains('isobaric')==False]
            delete = []
            for i in df_isobaric['mascot_peptide'].values:
                if i in df_other['mascot_peptide'].values:
                    delete.append(True)
                else:
                    delete.append(False)
            df_isobaric['del']=delete
            df_isobaric = df_isobaric[df_isobaric['del']==False]
            df_isobaric = df_isobaric.drop(['del'],axis=1)
            df_output = pd.concat([df_other,df_isobaric],ignore_index=True)
        
            for t in output_final:
                for i in t:
                    i = tuple(i)
                    if i not in recall:
                        continue
                    animal_pep = dfs[dfs['animal'].isin(i)]
                    df_output = pd.concat([df_output,animal_pep],ignore_index=True)
                    
            delete_animals = []
            for i in set(df_output['animal'].values):
                if len(set(df_output['found_match'][df_output['animal']==i]))==1:
                    delete_animals.append(i)
            df_output=df_output[df_output['animal'].isin(delete_animals)==False]
            #set the taxonomy
            found_animals = []
            for t in output_final:
                for i in t:
                    found_animals = found_animals+list(i)
        
            all_taxonomy = find_taxons(found_animals)
            if len(all_taxonomy)>0:
                final_output = {} #check tis block if it is correct!
                for ts in output_final:
                    for i in ts:
                        taxon = all_taxonomy[i[0]]
                        for t in i:
                            taxon = set(taxon)&set(all_taxonomy[t])
                        for t in all_taxonomy[i[0]]:
                            if t in taxon:
                                for x in i:
                                    final_output[x]=t[1]
                                break
                #If an animal is found that had an uncertain peptide match, than add this peptide to the output. All other uncertain peptides can be discarded
                keep_uncertain = keep_uncertain[keep_uncertain['animal'].isin(set(df_output['animal'].values))]
                keep_uncertain = keep_uncertain[df_output.columns]
                df_output = pd.concat([df_output,keep_uncertain],ignore_index=True)
                df_plot,df_output = make_connection_graph(df_output, file_name,final_output,found_animals)
                bc, labels,file_out_final = make_sunburst(dfs,all_animals,df_output,file_name,input_animals,df_plot,path,sample_path)
                ################################################################################
                df_adding = dfs[data_to_remember.columns]#keep outcome
                data_to_remember = pd.concat([data_to_remember,df_adding],ignore_index=True)
                ################################################################################
                to_summary=rescore(path,file_out_final,file_name,sample_path)
                general_summary_output_file[file_name]=[lim_tax,total_peptide_begin]+to_summary
            else:
                with open(path+'/Output_Classicol/'+sample_path+'/ZooMS_results_'+file_name+'.csv', 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                    writer.writerow(['ZooMS analysis found nothing'])
                    general_summary_output_file[file_name]=[lim_tax,total_peptide_begin]+['nothing_found']+['']*17
        else:
            with open(path+'/Output_Classicol/'+sample_path+'/ZooMS_results_'+file_name+'.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
                writer.writerow(['ZooMS analysis found nothing'])
                general_summary_output_file[file_name]=[lim_tax,total_peptide_begin]+['nothing_found']+['']*17
    
    #create summary output file
    today = str(date.today())
    with open(path+'/Output_Classicol/Summary_taxonomic_classification_'+today+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
        writer.writerow(['Experiment','Taxonomic_restriction','Total_peptide_count_in_input_file',
                         'Species(1)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count',
                         'Species(2)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count',
                         'Species(3)','Score','Rescore','Pep_count','isoBLAST_count','unique_peptide_count'])
        for k,v in general_summary_output_file.items():
            writer.writerow([k]+v)
    
    
    
    
    
    
    
    
    
    
    
    
