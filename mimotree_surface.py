import os
import pandas as pd
import math
from tqdm import tqdm
from biopandas.pdb import PandasPdb

#surface_dict = {}
GXG = {'ALA':113,'ARG':241,'ASN':158,'ASP':151,'CYS':140,
       'GLN':189,'GLU':183,'GLY':85,'HIS':194,'ILE':182,
       'LEU':180,'LYS':211,'MET':204,'PHE':218,'PRO':143,
       'SER':122,'THR':146,'TRP':259,'TYR':229,'VAL':160}

def not_empty(s):
    return s and s.strip()

#proteinSASA_txt is achieved from Chimera, the function returns the surface residues as a dictionary
'''def surface_residue_id(proteinSASA_txt,chain_identifier):
    residueSASA = {}
    surface_residue = {}
    SASA = open(proteinSASA_txt,'r')
    for x in SASA:
        y = x.split()
        r = y[0].split('.')
        if y[0] != 'attribute:' and y[0] != 'recipient:' and r[1] == chain_identifier:
            name = y[0]
            sasa = float(y[1])
            residueSASA[name] = sasa
    max_sasa = max(residueSASA.values())
    th = float(max_sasa)*0.05 #surface exposed residue surface area threshold
    for key,value in residueSASA.items():
            if value > th:
                surface_residue[key] = value
    return(surface_residue)'''

def surface_residue_id(pdb, chain_identifier, proteinSASA_txt):    
    ppdb = PandasPdb().read_pdb(pdb)
    pdb_df = ppdb.df['ATOM']
    pdb_df = pdb_df[pdb_df['chain_id'] == chain_identifier]
    residueSASA = {}
    surface_residue = {}
    SASA = open(proteinSASA_txt,'r')
    for x in SASA:
        y = x.split()
        r = y[0].split('.')
        if y[0] != 'attribute:' and y[0] != 'recipient:' and r[1] == chain_identifier:
            name = y[0]
            sasa = float(y[1])
            residueSASA[name] = sasa
    max_sasa = max(residueSASA.values())
    for key,value in residueSASA.items():
        key_n = int(''.join(list(filter(str.isdigit, key))))
        AA = pdb_df[pdb_df['residue_number']==key_n]['residue_name'].tolist()[0]
        if value > max_sasa*0.05:
            surface_residue[key] = value
    return surface_residue, pdb_df

def generate_df(pdb_df, surface_residue, chain_identifier):
    surface_residue_list = list(surface_residue.keys())
    surface_residue_num = [x.strip(':.'+chain_identifier) for x in surface_residue_list]
    surface_residue_num = [int(y) for y in surface_residue_num]
    pdb_df = pdb_df.loc[pdb_df['residue_number'].isin(surface_residue_num)]
    return pdb_df

def initial_graph(surface_residue, pdb_df, chain_identifier):
    surface_dict = {}
    surface_residue_list = list(surface_residue.keys())
    surface_residue_num = [x.strip(':.'+chain_identifier) for x in surface_residue_list]
    for y in surface_residue_num:
        residue_n = list(pdb_df.loc[pdb_df['residue_number'] == float(y)]['residue_name'])[0]
        name = residue_n + y
        surface_dict[name] = []
    return surface_dict

'''
#surface_residue = surface_residue_id('1e6jSASA.txt','P')
# pdb_txt: original pdb text file
# surface_residue: surface residue dictionary achieved by func surface_residue_id()
# chain_identifier: e.q. 'P'
def Initialize_graph(pdb_txt,surface_residue,chain_identifier):
    surface_dict = {}
    surface_residue_pdb_txt = open('surface_residue_pdb.txt','w')
    with open(pdb_txt,'r') as r:
        lines = r.readlines()
    for key in surface_residue.keys():
        y = key.strip(':.'+chain_identifier)
        for l in lines:
            row = l.split(" ")
            if row[0] == 'ATOM':
                f = filter(not_empty,row)
                row = list(f)
                if row[4] == chain_identifier and row[5] == y:
                    surface_residue_pdb_txt.write(l)
    surface_residue_pdb_txt.close()
    residue = open('surface_residue_pdb.txt', "r")
    for x in residue:
        y = x.split(" ")
        if y[0] == 'ATOM':
            f = filter(not_empty,y)
            y = list(f)
            name = y[3]+y[5]
            surface_dict[name] = []
    residue.close()
    return(surface_dict)
'''

'''
def generate_df(surface_residue_pdb):
    txt_file = surface_residue_pdb
    base = os.path.splitext(txt_file)[0]
    os.rename(txt_file, base + '.pdb')
    ppdb = PandasPdb().read_pdb('surface_residue_pdb.pdb')
    pdb_df = ppdb.df['ATOM']
    return pdb_df
'''

#the minimum distance between the atoms of two surface residues
def residue_distance(n1,n2,pdb_df):
    distance = []
    residue1 = pdb_df[pdb_df['residue_number']==n1]
    residue2 = pdb_df[pdb_df['residue_number']==n2]
    r1 = len(residue1)
    r2 = len(residue2)
    for i in range(r1):
        content1 = residue1.iloc[i].tolist()
        x1, y1, z1 = float(content1[11]), float(content1[12]), float(content1[13])
        for j in range(r2):
            content2 = residue2.iloc[j].tolist()
            x2, y2, z2 = float(content2[11]), float(content2[12]), float(content2[13])
            result = math.sqrt(math.pow(x1-x2, 2) + math.pow(y1-y2, 2) + math.pow(z1-z2, 2))
            distance.append(result)
    return content1[5]+str(content1[8])+' '+content2[5]+str(content2[8])+' '+ '{0:.2f}'.format(min(distance))

#identify neighbors with distance threshold
#distance_param = 4
def surface_map(pairs_distance_txt, distance_param, surface_dict):
    pairs_distance = open(pairs_distance_txt, "r")
    for x in pairs_distance:
        y = x.split()
        if float(y[2]) < distance_param:
            surface_dict[y[0]].append(y[1])
            surface_dict[y[1]].append(y[0])
    return surface_dict


