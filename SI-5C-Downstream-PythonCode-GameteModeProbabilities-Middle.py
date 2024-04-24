import sympy as sp 
import itertools

i_m = sp.symbols('i_m')
j_m = sp.symbols('j_m')
k_m = sp.symbols('k_m')
l_m = sp.symbols('l_m')
G = sp.symbols('G')
g = sp.symbols('g')
a1,a2,b1,b2,a3,a4,b3,b2,b4,c1,c2,c3,c4 = sp.symbols('a1,a2,b1,b2,a3,a4,b3,b2,b4,c1,c2,c3,c4')

modes = [[[i_m, G, i_m], [i_m, G, i_m]],
 [[i_m, G, j_m], [i_m, G, j_m]],
 [[i_m, G, i_m], [i_m, G, j_m]],
 [[i_m, G, j_m], [i_m, G, k_m]],
 [[i_m, G, i_m], [j_m, G, i_m]],
 [[i_m, G, j_m], [k_m, G, j_m]],
 [[i_m, G, i_m], [j_m, G, j_m]],
 [[i_m, G, i_m], [j_m, G, k_m]],
 [[i_m, G, j_m], [j_m, G, i_m]],
 [[i_m, G, j_m], [j_m, G, k_m]],
 [[i_m, G, j_m], [k_m, G, l_m]],
 [[i_m, G, i_m], [i_m, g, j_m]],
 [[i_m, G, i_m], [j_m, g, i_m]],
 [[i_m, G, i_m], [j_m, g, k_m]],
 [[i_m, G, j_m], [j_m, g, k_m]],
 [[i_m, g, i_m], [i_m, G, i_m]],
 [[i_m, g, j_m], [i_m, G, j_m]],
 [[i_m, g, i_m], [i_m, G, j_m]],
 [[i_m, g, j_m], [i_m, G, k_m]],
 [[i_m, g, i_m], [j_m, G, i_m]],
 [[i_m, g, j_m], [k_m, G, j_m]],
 [[i_m, g, i_m], [j_m, G, j_m]],
 [[i_m, g, i_m], [j_m, G, k_m]],
 [[i_m, g, j_m], [j_m, G, i_m]],
 [[i_m, g, j_m], [j_m, G, k_m]],
 [[i_m, g, j_m], [k_m, G, l_m]],
 [[i_m, g, i_m], [i_m, g, i_m]],
 [[i_m, g, j_m], [i_m, g, j_m]],
 [[i_m, g, i_m], [i_m, g, j_m]],
 [[i_m, g, j_m], [i_m, g, k_m]],
 [[i_m, g, i_m], [j_m, g, i_m]],
 [[i_m, g, j_m], [k_m, g, j_m]],
 [[i_m, g, i_m], [j_m, g, j_m]],
 [[i_m, g, i_m], [j_m, g, k_m]],
 [[i_m, g, j_m], [j_m, g, i_m]],
 [[i_m, g, j_m], [j_m, g, k_m]],
 [[i_m, g, j_m], [k_m, g, l_m]]]
 
def num_letters_not_Gg(m):
    tmp = [m[0][0],m[0][1],m[0][2],m[1][0],m[1][1],m[1][2]]
    tmp_set = set(tmp) - set([G,g])
    tmp_lst = list(tmp_set)
    
    return len(tmp_lst)
 
modes_w_numbers = []
modes_w_numbers_id = []
mode_count = -1
for mode in modes:
    chr_count = num_letters_not_Gg(mode)
    combs = list(itertools.permutations(range(4),chr_count))
    mode_count = mode_count + 1
    
    for comb in combs:
        if chr_count == 1:
            tmp_mode = [[mode[0][0],mode[0][1],mode[0][2]],[mode[1][0],mode[1][1],mode[1][2]]]
            for sub_mode in tmp_mode:
                indexes = [i for i in range(len(sub_mode)) if sub_mode[i] in [i_m]]
                for index in indexes:
                    sub_mode[index] = comb[0]
            modes_w_numbers.append(tmp_mode)
            modes_w_numbers_id.append(mode_count)
        elif chr_count == 2:
            tmp_mode = [[mode[0][0],mode[0][1],mode[0][2]],[mode[1][0],mode[1][1],mode[1][2]]]
            for sub_mode in tmp_mode:
                indexes_i_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [i_m]]
                indexes_j_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [j_m]]
                for index in indexes_i_m:
                    sub_mode[index] = comb[0]
                for index in indexes_j_m:
                    sub_mode[index] = comb[1]
            modes_w_numbers.append(tmp_mode)
            modes_w_numbers_id.append(mode_count)
        elif chr_count == 3:
            tmp_mode = [[mode[0][0],mode[0][1],mode[0][2]],[mode[1][0],mode[1][1],mode[1][2]]]
            for sub_mode in tmp_mode:
                indexes_i_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [i_m]]
                indexes_j_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [j_m]]
                indexes_k_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [k_m]]
                for index in indexes_i_m:
                    sub_mode[index] = comb[0]
                for index in indexes_j_m:
                    sub_mode[index] = comb[1]
                for index in indexes_k_m:
                    sub_mode[index] = comb[2]
            modes_w_numbers.append(tmp_mode)
            modes_w_numbers_id.append(mode_count)
        elif chr_count == 4:
            tmp_mode = [[mode[0][0],mode[0][1],mode[0][2]],[mode[1][0],mode[1][1],mode[1][2]]]
            for sub_mode in tmp_mode:
                indexes_i_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [i_m]]
                indexes_j_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [j_m]]
                indexes_k_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [k_m]]
                indexes_l_m = [i for i in range(len(sub_mode)) if sub_mode[i] in [l_m]]
                for index in indexes_i_m:
                    sub_mode[index] = comb[0]
                for index in indexes_j_m:
                    sub_mode[index] = comb[1]
                for index in indexes_k_m:
                    sub_mode[index] = comb[2]
                for index in indexes_l_m:
                    sub_mode[index] = comb[3]
            modes_w_numbers.append(tmp_mode)
            modes_w_numbers_id.append(mode_count)
        else:
            print('Error')

modes_diallelic = []

for gam in modes_w_numbers:
    tmp_gam = []
    for chr in gam:
        tmp_chr = []
        locus = -1
        for allele in chr:
            tmp_allele = allele
            locus = locus + 1
            if locus == 0:
                if allele == 0:
                    tmp_allele = a1
                elif allele == 1:
                    tmp_allele = a2
                elif allele == 2:
                    tmp_allele = a3
                elif allele == 3:
                    tmp_allele = a4
                else:
                    print('error')
            elif locus == 1:
                if allele == G:
                    tmp_allele = G
                elif allele == g:
                    tmp_allele = g
                else:
                    print('error')
            elif locus == 2:
                if allele == 0:
                    tmp_allele = c1
                elif allele == 1:
                    tmp_allele = c2
                elif allele == 2:
                    tmp_allele = c3
                elif allele == 3:
                    tmp_allele = c4
                else:
                    print('error')
            else:
                print('error')
            tmp_chr.append(tmp_allele)
        tmp_gam.append(tmp_chr)
    modes_diallelic.append(tmp_gam)
 
r= sp.symbols('r')
q=sp.symbols('q')
v=sp.symbols('v')
vp=sp.symbols('vp')
rp=sp.symbols('rp')
qp=sp.symbols('qp')
t=sp.symbols('t')
(ppc)=sp.symbols('(ppc)')
(pdm)=sp.symbols('(pdm)')
(pmp)=sp.symbols('(pmp)')
p_pair=sp.symbols('p_pair')

G, g = sp.symbols('G, g')

import copy
from sympy import *
from itertools import permutations 
from itertools import chain
from itertools import product

#SWAPS OUTTER SET -C
def swap_r1(z, i, j, x, y=2):
    tmp=copy.deepcopy(z)
    tmp[i][x][y], tmp[j][x][y] =  tmp[j][x][y], tmp[i][x][y]
    return tmp

#SWAPS MIDDLE SET-B
def swap_r2(z, i,j,x,y=1):
    tmp=copy.deepcopy(z)
    tmp[i][x][y:], tmp[j][x][y:] =  tmp[j][x][y:], tmp[i][x][y:]
    return tmp

#SWAPS INNER SET-A
def swap_r3(z, i,j,x):
    tmp=copy.deepcopy(z)
    tmp[i][x], tmp[j][x] =  tmp[j][x], tmp[i][x]
    return tmp

def biv_D(arrng,p_pair,Gg_arrng):
    events_D_list=[[]]
    dg_D=[arrng]
    
    r_events = [[],[r,0,2,0],[r,0,2,1],[q,0,2,0],[q,0,2,1],[v,0,2,0],[v,0,2,1], \
         [r,1,3,0],[r,1,3,1],[q,1,3,0],[q,1,3,1],[v,1,3,0],[v,1,3,1]]
    
    for event in [1,2,3,4,5,6,7,8,9,10,11,12]:
        tmp_events_D_list0=copy.deepcopy(events_D_list)
        tmp_events_D_list1=copy.deepcopy(events_D_list)
        
        for recomb in [0,1]:
            if (recomb == 0):
                if (event <= 2 or event == 7 or event == 8):
                    p = (1 - r_events[event][0])
                elif (3 <= event <= 6 or event >= 9) and Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = (1 - r_events[event][0])
                else:
                    p = (1 - r_events[event][0] * p_pair) 
                
                dg_0=dg_D 
                for i in tmp_events_D_list0:
                    i.append(p)
                    
            else:
                if (event <= 2 or event == 7 or event == 8):
                    p = r_events[event][0]
                elif (3 <= event <= 6 or event >= 9) and Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = r_events[event][0]
                else:
                    p = r_events[event][0] * p_pair
                
                for i in range(len(dg_D)):
                    tmp = dg_D[i]
                    
                    if r_events[event][0] == r:
                        dg = swap_r1(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_D.append(dg)
                    elif r_events[event][0] == q:
                        dg = swap_r2(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_D.append(dg)
                    elif r_events[event][0] == v:
                        dg = swap_r3(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_D.append(dg)    
                for i in tmp_events_D_list1:
                    i.append(p)
                    
        events_D_list=[]
        for i in tmp_events_D_list0:
            events_D_list.append(i)
        for i in tmp_events_D_list1:
            events_D_list.append(i)

    product_D_list=[]
    for i in events_D_list:
        SUM = 1
        for j in i:
            SUM *= j
        product_D_list.append([SUM])
    product_D_list=[[item * ((1-t) + t*(1-pmp-ppc-pdm))/24 for item in subl][0] for subl in product_D_list]
    
    return dg_D, product_D_list

def quad_A_D(arrng,p_pair,Gg_arrng):
    
    events_A_list = [[]]
    dg_A =[arrng]
    
    r_events = [[],[r,0,2,0],[r,0,2,1],[q,0,2,0],[q,0,2,1],[vp,0,2,0],\
        [vp,0,2,1],[r,1,3,0],[r,1,3,1],[q,1,3,0],[q,1,3,1],\
        [vp,1,3,0],[vp,1,3,1],[vp,0,1,0],[vp,0,1,1],[vp,2,3,0],\
        [vp,2,3,1]]
    
    for event in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        tmp_events_A_list0=copy.deepcopy(events_A_list)
        tmp_events_A_list1=copy.deepcopy(events_A_list)
        
        for recomb in [0,1]:
            if (recomb == 0):
                if (event <= 2 or event == 7 or event == 8):
                    p = (1 - r_events[event][0])
                elif (3 <= event <= 6 or 9 <= event <= 12) and \
                            Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = (1 - r_events[event][0])
                elif event > 12:
                    p = (1 - r_events[event][0])
                else:
                    p = (1 - r_events[event][0] * p_pair)
                    
                dg_0 = dg_A
                for i in tmp_events_A_list0:
                    i.append(p)
            else:
                if (event <= 2 or event == 7 or event == 8):
                    p = r_events[event][0]
                elif (3 <= event <= 6 or 9 <= event <= 12) and \
                            Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = r_events[event][0]
                elif event > 12:
                    p = r_events[event][0]
                else:
                    p = r_events[event][0] * p_pair
                    
                for i in range(len(dg_A)):
                    tmp=dg_A[i]
                    
                    if r_events[event][0] == r:
                        dg = swap_r1(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_A.append(dg)
                    elif r_events[event][0] == q:
                        dg = swap_r2(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_A.append(dg)
                    elif r_events[event][0] == v or r_events[event][0] == vp:
                        dg = swap_r3(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_A.append(dg)
                        
                for i in tmp_events_A_list1:
                    i.append(p)
                    
        events_A_list=[]
        for i in tmp_events_A_list0:
            events_A_list.append(i)
        for i in tmp_events_A_list1:
            events_A_list.append(i)

    product_A_list=[]
    for i in events_A_list:
        SUM = 1
        for j in i:
            SUM *= j
        product_A_list.append([SUM])
    product_A_list1=[[item *(t)*(ppc)/24 for item in subl][0] for subl in product_A_list]

    return dg_A,  product_A_list1

def quad_B_D(arrng,p_pair,Gg_arrng):
    events_b_list = [[]]
    dg_B = [arrng]
    
    r_events = [[],[rp,0,2,0],[rp,0,2,1],[rp,1,3,0],[rp,1,3,1],\
        [rp,0,1,0],[rp,0,1,1],[q,0,1,0],[q,0,1,1],\
        [v,0,1,0],[v,0,1,1],[rp,2,3,0],[rp,2,3,1],\
        [q,2,3,0],[q,2,3,1],[v,2,3,0],[v,2,3,1]]
    
    for event in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        tmp_events_b_list_0=copy.deepcopy(events_b_list)
        tmp_events_b_list_1=copy.deepcopy(events_b_list)
        
        for recomb in [0,1]:
            if (recomb == 0):
                if (event <= 6 or event == 11 or event == 12):
                    p = (1 - r_events[event][0])
                elif (7 <= event <= 10 or event > 12) and \
                                Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = (1 - r_events[event][0])
                else:
                    p = (1 - r_events[event][0] * p_pair)
                    
                dg_0 = dg_B
                for i in tmp_events_b_list_0:
                    i.append(p)
            else:
                if (event <= 6 or event == 11 or event == 12):
                    p = r_events[event][0]
                elif (7 <= event <= 10 or event > 12) and \
                                Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = r_events[event][0]
                else:
                    p = r_events[event][0] * p_pair
                    
                for i in range(len(dg_B)):
                    tmp=dg_B[i]
                    
                    if r_events[event][0] == r or r_events[event][0] == rp:
                        dg = swap_r1(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_B.append(dg)
                    elif r_events[event][0] == q:
                        dg = swap_r2(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_B.append(dg)
                    elif r_events[event][0] == v:
                        dg = swap_r3(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_B.append(dg)
                    
                for i in tmp_events_b_list_1:
                    i.append(p)
                    
        events_b_list=[]
        for i in tmp_events_b_list_0:
            events_b_list.append(i)
        for i in tmp_events_b_list_1:
            events_b_list.append(i)
    
    product_b_list=[]
    for i in events_b_list:
        SUM = 1
        for j in i:
            SUM *= j
        product_b_list.append([SUM])
    product_b_list1=[[item *(t)*(pdm)/24 for item in subl][0] for subl in product_b_list]
    
    return dg_B, product_b_list1

def quad_C_D(arrng,p_pair,Gg_arrng):
    events_c_list = [[]]
    dg_C =[arrng]
    
    r_events = [[],[r,0,2,0],[r,0,2,1],[qp,0,2,0],[qp,0,2,1],\
        [r,1,3,0],[r,1,3,1],[qp,1,3,0],[qp,1,3,1],\
        [qp,0,1,0],[qp,0,1,1],[v,0,1,0],[v,0,1,1],\
        [qp,2,3,0],[qp,2,3,1],[v,2,3,0],[v,2,3,1]]
    
    for event in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]:
        tmp_events_c_list_0=copy.deepcopy(events_c_list)
        tmp_events_c_list_1=copy.deepcopy(events_c_list)
        for recomb in [0,1]:
            if (recomb == 0):
                if (event <= 2 or event == 5 or event == 6):
                    p = (1 - r_events[event][0])
                elif (3 <= event <= 4 or 7 <= event <= 8) and \
                            Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = (1 - r_events[event][0])
                elif (9 <= event <= 12 or 13 <= event <= 16):
                    p = (1 - r_events[event][0])
                else:
                    p = (1 - r_events[event][0] * p_pair)
                
                dg_0 = dg_C
                for i in tmp_events_c_list_0:
                    i.append(p)
            else:
                if (event <= 2 or event == 5 or event == 6):
                    p = r_events[event][0]
                elif (3 <= event <= 4 or 7 <= event <= 8) and \
                            Gg_arrng[r_events[event][1]] == Gg_arrng[r_events[event][2]]:
                    p = r_events[event][0]
                elif (9 <= event <= 12 or 13 <= event <= 16):
                    p = r_events[event][0]
                else:
                    p = r_events[event][0] * p_pair
                    
                for i in range(len(dg_C)):
                    tmp=dg_C[i]
                    
                    if r_events[event][0] == r:
                        dg = swap_r1(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_C.append(dg)
                    elif r_events[event][0] == q or r_events[event][0] == qp:
                        dg = swap_r2(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_C.append(dg)
                    elif r_events[event][0] == v:
                        dg = swap_r3(tmp, r_events[event][1], r_events[event][2],r_events[event][3])
                        dg_C.append(dg)  
                    
                for i in tmp_events_c_list_1:
                    i.append(p)
                    
        events_c_list=[]
        for i in tmp_events_c_list_0:
            events_c_list.append(i)
        for i in tmp_events_c_list_1:
            events_c_list.append(i)
    
    product_c_list=[]
    for i in events_c_list:
        SUM = 1
        for j in i:
            SUM *= j
        product_c_list.append([SUM])
    product_c_list1=[[item *(t)*(pmp)/24 for item in subl][0] for subl in product_c_list]
    
    return dg_C, product_c_list1

def ana(r_list,is_biv):
    
    if is_biv == 0:
        ana_events = [[0,0,2,0],[0,1,2,1],[0,0,2,1],[0,1,2,0],[1,0,3,0],\
            [1,1,3,1],[1,0,3,1],[1,1,3,0],[0,0,3,0],[0,1,3,1],\
            [0,0,3,1],[0,1,3,0],[1,0,2,0],[1,1,2,1],[1,0,2,1],\
            [1,1,2,0]]
    else:
        ana_events = [[0,0,1,0],[0,1,1,1],[0,0,1,1],[0,1,1,0],[2,0,3,0],\
            [2,1,3,1],[2,0,3,1],[2,1,3,0],[0,0,3,0],[0,1,3,1],\
            [0,0,3,1],[0,1,3,0],[2,0,1,0],[2,1,1,1],[2,0,1,1],\
            [2,1,1,0]]
    
    pair = [[] for event in range(len(ana_events))]
    
    for event in range(len(ana_events)):
        for i in range(len(r_list)):
            pair[event].append([copy.deepcopy(r_list[i][ana_events[event][0]][ana_events[event][1]]),\
                                            copy.deepcopy(r_list[i][ana_events[event][2]][ana_events[event][3]])])
    
    total_gamete_list = [copy.deepcopy(pair[event][gamete]) for event in range(len(ana_events)) for gamete in range(len(pair[event]))]
    
    return total_gamete_list

def total_list(lists):
    prob=[[item / 16 for item in [subl]] for subl in lists]
    prod=list(map(sum, prob))
    total_prob=prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod+prod
    return total_prob

a, b, c = sp.symbols('a, b, c', cls=Function)

def gam_mode(gametes, gamete_probabilities, mode_list):
    gametes_in_mode = [[] for i in range(len(mode_list))]
    prob_modes = [[] for i in range(len(mode_list))]

    for j in range(len(mode_list)):
        for i in range(len(gametes)):
            #print(j,i)
            if gametes[i] in [mode_list[j]]:
                gametes_in_mode[j].append(gametes[i])
                prob_modes[j].append(gamete_probabilities[i])
    return [sum(prob_modes[i]) for i in range(len(prob_modes))], \
                        [len(gametes_in_mode[i]) for i in range(len(gametes_in_mode))]

def mode_prob(x, y, z, perm, p_pair, Gg_perm):
    
    if x==0:
        
        D=biv_D(perm[z],p_pair, Gg_perm[z])
        A=ana(D[0],1)
        T=total_list(D[1])
        G=gam_mode(A, T, modes_diallelic)
        return G
                           
    else:
        if y==0:
            
            A=quad_A_D(perm[z],p_pair, Gg_perm[z])
            AA=ana(A[0],0)
            T=total_list(A[1])
            G=gam_mode(AA, T, modes_diallelic)
            return G
        
        if y==1:
            
            B=quad_B_D(perm[z],p_pair, Gg_perm[z])
            A=ana(B[0],0)
            T=total_list(B[1])
            G=gam_mode(A, T, modes_diallelic)
            return G
        
        if y==2:
            
            C=quad_C_D(perm[z],p_pair, Gg_perm[z])
            A=ana(C[0],0)
            T=total_list(C[1])
            G=gam_mode(A, T, modes_diallelic)
            return G

B_base=[[[a1,c1], [a1,c1]],[[a2,c2], [a2,c2]], [[a3,c3], [a3,c3]],[[a4,c4], [a4,c4]]]

#Gg_lst = list(product([[G,G],[g,g]], repeat=4))
Gg_lst = [([G, G], [G, G], [g, g], [g, g])]

B_D = []
B_indx = 0

for geno in Gg_lst:
    B_D.append([])
    item_indx_1 = 0
    for i in range(len(B_base)):
        B_D[B_indx].append([])
        item_indx_2 = 0
        for j in range(len(B_base[i])):
            B_D[B_indx][i].append([B_base[i][j][0],geno[item_indx_1][item_indx_2],B_base[i][j][1]])
            item_indx_2 = item_indx_2 + 1
        item_indx_1 = item_indx_1 + 1
    B_indx = B_indx + 1
    
perm_D=[[list(p) for p in itertools.permutations(item)] for item in B_D ]
Gg_perm = [[list(p) for p in itertools.permutations(item)] for item in Gg_lst ]

comb_biv = [[i,j] for i in range(len(Gg_lst)) for j in range(24) ]
comb_quad = [[i,j] for i in range(len(Gg_lst)) for j in range(24) ]

import pickle

biv_res = [mode_prob(0, 0, comb_biv[i][1], perm_D[comb_biv[i][0]], p_pair, \
                     Gg_perm[comb_biv[i][0]]) \
                                                for i in range(len(comb_biv))]

with open("biv_res.out", 'wb') as outf:
    pickle.dump(biv_res, outf)

#quad_A_res = [mode_prob(1, 0, comb_quad[i][1], perm_D[comb_quad[i][0]], p_pair, \
#                     Gg_perm[comb_quad[i][0]]) \
#                                                for i in range(len(comb_quad))]

#with open("quad_A_res.out", 'wb') as outf:
#    pickle.dump(quad_A_res, outf)
                                                
#quad_B_res = [mode_prob(1, 1, comb_quad[i][1], perm_D[comb_quad[i][0]], p_pair, \
#                     Gg_perm[comb_quad[i][0]]) \
#                                                for i in range(len(comb_quad))]

#with open("quad_B_res.out", 'wb') as outf:
#    pickle.dump(quad_B_res, outf)

#quad_C_res = [mode_prob(1, 2, comb_quad[i][1], perm_D[comb_quad[i][0]], p_pair, \
#                     Gg_perm[comb_quad[i][0]]) \
#                                                for i in range(len(comb_quad))]

#with open("quad_C_res.out", 'wb') as outf:
#    pickle.dump(quad_C_res, outf)
