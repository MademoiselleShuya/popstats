# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 16:49
@Author : Shuya Zhang, Ngenie Liang
@File : computing_functions.py
@Project : popstats
"""
import math
import sys
import random

import numpy as np

validmarkers = ['A', 'T', 'G', 'C', '0']


def FST_W_pairwise(col):
    NpopA = float(col[0])
    NpopB = float(col[2])

    popAcount = int(col[1])
    popBcount = int(col[3])

    npops = 2.0
    # nsamples = float(NpopA + NpopB) useless
    n_bar = (NpopA / npops) + (NpopB / npops)
    samplefreq = ((popAcount + popBcount) / (NpopA + NpopB))
    pop1freq = popAcount / float(NpopA)
    pop2freq = popBcount / float(NpopB)
    Npop1 = NpopA
    Npop2 = NpopB
    S2A = (1 / ((npops - 1.0) * n_bar)) * (
            ((Npop1) * ((pop1freq - samplefreq) ** 2)) + ((Npop2) * ((pop2freq - samplefreq) ** 2)))
    nc = 1.0 / (npops - 1.0) * ((Npop1 + Npop2) - (((Npop1 ** 2) + (Npop2 ** 2)) / (Npop1 + Npop2)))
    T_1 = S2A - ((1 / (n_bar - 1)) * ((samplefreq * (1 - samplefreq)) - ((npops - 1) / npops) * S2A))
    T_2 = (((nc - 1) / (n_bar - 1)) * samplefreq * (1 - samplefreq)) + (
            1.0 + (((npops - 1) * (n_bar - nc)) / (n_bar - 1))) * (S2A / npops)

    return T_1, T_2


def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)


def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)


def FST_H_pairwise(col):
    NpopA = float(col[0])
    NpopB = float(col[2])

    popAcount = int(col[1])
    popBcount = int(col[3])

    pop1freq = popAcount / float(NpopA)
    pop2freq = popBcount / float(NpopB)
    Npop1 = NpopA
    Npop2 = NpopB
    T_1 = (pop1freq - pop2freq) ** 2 - ((pop1freq * (1.0 - pop1freq)) / (Npop1 - 1)) - (
            (pop2freq * (1.0 - pop2freq)) / (Npop2 - 1))
    T_2 = (pop1freq * (1.0 - pop2freq)) + (pop1freq * (1.0 - pop2freq))

    return T_1, T_2


def selectiongenotypes(allgenos, targetlist, chromosome, position):
    returnstring = ''
    for i in targetlist:
        try:
            geno = allgenos[i]
        except IndexError:
            print('IndexError', i, len(allgenos) - 1)
            print(chromosome, position, allgenos)
            exit(0)
        if geno != '0' and geno != 'N':
            returnstring += geno
    return returnstring


def selectmultigenotypes(allgenos, targetlist, chromosome, position):
    returnstring = []
    for i in targetlist:
        try:
            geno = allgenos[i]
        except IndexError:
            print('IndexError', i, len(allgenos) - 1)
            print(chromosome, position, allgenos)
            exit(0)
        if '0' not in geno and 'N' not in geno:
            returnstring.append(geno)
    return returnstring


def selectionhaplotypes(allgenos, targetlist, chromosome, position):
    returnstring = []
    for i in targetlist:
        try:
            geno = allgenos[i]
        except IndexError:
            print('IndexError', i, len(allgenos) - 1, file=sys.stderr)
            print(chromosome, position, allgenos, file=sys.stderr)
            exit(0)
        if geno in validmarkers:
            returnstring.append(geno)
        else:
            returnstring.append('0')
    return returnstring


def scramblefun(toscramble):
    returngenos = []
    for i in range(0, len(toscramble), 2):
        sg = toscramble[i:i + 2]
        sg = random.sample(sg, 2)
        returngenos.append(sg[0])
        returngenos.append(sg[1])
    return returngenos


def haploidizefun(tohaploidize):
    hreturngenos = []
    for i in range(0, len(tohaploidize), 2):
        hg = tohaploidize[i:i + 2]
        hg = random.sample(hg, 2)
        hreturngenos.append(hg[0])
    return hreturngenos


def Bcorrfun(inputlist):
    """dictionaries for t and n where keys are B-stat"""
    tdict = {}
    ndict = {}
    for line in inputlist:
        # print line
        bstat = line[4]
        tstat = line[2]
        nstat = line[3]
        if bstat in list(tdict.keys()):
            addition = tdict[bstat]
            addition += tstat
            tdict[bstat] = addition
            addition = ndict[bstat]
            addition += nstat
            ndict[bstat] = addition
        else:
            tdict[bstat] = tstat
            ndict[bstat] = nstat

    blist = []
    Dlist = []
    bnsum = 0
    for i in list(tdict.keys()):
        blist.append(i)
        bDstat = 1.0 * tdict[i] / ndict[i]
        bnsum = ndict[i]
        Dlist.append(bDstat)
    bmean = sum(blist) / len(blist)
    x = [float(b - bmean) for b in blist]
    y = Dlist
    # corrstat = pearson_def(x, y)

    corrstat = np.polyfit(x, y, 1)
    corrstat = corrstat[0]
    return [corrstat, bnsum]


def mutate(mutateinp, mutationrate, mutalleles):
    mutateoutp = ''
    for mbase in mutateinp:
        if random.random() <= mutationrate:
            newbase = [mbas for mbas in mutalleles if mbas != mbase]
            mutateoutp += newbase[0]
        else:
            mutateoutp += mbase
    return mutateoutp


def Bstratfun(inputlist):
    """dictionaries for t and n where keys are B-stats"""
    hight = 0
    highn = 0
    lowt = 0
    lown = 0

    for line in inputlist:
        # print line
        bstat = line[4]
        tstat = line[2]
        nstat = line[3]
        if bstat in [0, 1]:
            lowt += tstat
            lown += nstat
        elif bstat in [8, 9]:
            hight += tstat
            highn += nstat

    highstat = hight / highn
    lowstat = lowt / lown
    Bdiff = lowstat - highstat
    return Bdiff
