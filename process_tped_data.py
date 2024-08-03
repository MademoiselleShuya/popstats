# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 16:58
@Author : Shuya Zhang, Ngenie Liang
@File : process_tped_data.py
@Project : popstats
"""
import random
import math

from utils import progress
from computing_functions import scramblefun, selectmultigenotypes, selectionhaplotypes, \
    selectiongenotypes, haploidizefun, mutate, FST_H_pairwise, FST_W_pairwise

ILSdevoid = """10241177 12619185
16946047 18747389
19303480 22198160
38344992 41272675
45930478 77954462
99459295 111145964
128232540 136796526
151519514 155156362"""
ILSdevoid = [i.split() for i in ILSdevoid.split('\n')]

Neantroughs = """10000000 20000000
40000000 50000000
70000000 80000000
100000000 110000000
120000000 130000000
130000000 140000000
20000000 30000000"""
Neantroughs = [i.split() for i in Neantroughs.split('\n')]


def process_tped_data(data, options, regionchoice, popcount, targetpop1, targetpop2, targetpop3, targetpop4,
                      targetpop1count, targetpop2count, targetpop3count, targetpop4count, targetpop21,
                      targetpop22, targetpop23, targetpop24, ancpop, testpop, ascpop, ascpopcount, ascpop2,
                      ascpopcount2, ascpop3, ascpopcount3, mutationclass):
    line_counter = 0
    previouschromosome = '0'
    threshold = 0.1
    minsamplesize = 6

    previouslines = []
    t_list = []
    n_list = []
    choice_list = []
    abbalist = []
    babalist = []
    maf_list = []

    maf = False
    notransitions = options.notransitions

    for line in open(data):
        line_counter += 1
        if options.clock:
            if options.outfile == False and line_counter > 999 and str(line_counter)[-3:] == '000':
                progress(line_counter)
        col = line.split()
        chromosome = int(col[0])
        rsid = col[1]
        position = int(float(col[3].split('-')[0]))
        if options.chromosome:
            if chromosome != options.chromosome:
                # print line,
                continue
            elif chromosome == 23 and (options.not23 == False):
                if position < 2699520:
                    continue  # Filters out PAR1 in GRCh37
                elif (154931044 < position) and (position < 155260560):
                    continue  # Filters out PAR2 in GRCh37

                if options.Dutheilfilter:
                    ILSdevoidstatus = False
                    for st, end in ILSdevoid:
                        if position > int(st) and position < int(end):
                            ILSdevoidstatus = True
                            break
                    if ILSdevoidstatus: continue

                if options.Sriramfilter:
                    Sriramstatus = False
                    for st, end in Neantroughs:
                        if position > int(st) and position < int(end):
                            Sriramstatus = True
                            break
                    if Sriramstatus: continue
        else:
            if (chromosome < 1) or (chromosome > 22) and (options.not23 == False): continue

        if options.region:
            if position < regionchoice[0]: continue
            if position > regionchoice[1]: continue

        if options.Bcorr or options.Bstrat:
            Bstat = int(rsid.split(',')[1])
            if options.verbose:
                print(Bstat, end=' ')
        if options.Bchoice:
            Bstat = int(rsid.split(',')[1])
            if Bstat != int(options.Bchoice):
                continue
        if options.SNPfreq:
            if position != int(options.SNPfreq):
                continue
        if options.morgan:
            position = float(col[2])

        genotypes = col[4:]

        if options.scrambleall:
            genotypes = scramblefun(genotypes)

        genocount = len(genotypes)
        if genocount != popcount:
            print('genocount != popcount,', genocount, popcount)
            exit(0)

        if options.multi:
            genotypes1 = selectmultigenotypes(genotypes, targetpop1, chromosome, position)
            genotypes2 = selectmultigenotypes(genotypes, targetpop2, chromosome, position)
            genotypes3 = selectmultigenotypes(genotypes, targetpop3, chromosome, position)
            genotypes4 = selectmultigenotypes(genotypes, targetpop4, chromosome, position)
            if len(genotypes1) < 1: continue
            if len(genotypes2) < 1: continue
            if len(genotypes3) < 1: continue
            if len(genotypes4) < 1: continue

            alleles = list(set(genotypes1 + genotypes2 + genotypes3 + genotypes4))
            alleles = [a for a in alleles if '0' not in a]
            if len(alleles) < 2: continue

            t_locus = []
            n_locus = []
            p1, p2, p3, p4 = 0, 0, 0, 0
            for a in alleles[1:]:
                ref_allele = a
                p1 = (1.0000 * genotypes1.count(ref_allele)) / float(len(genotypes1))
                p2 = (1.0000 * genotypes2.count(ref_allele)) / float(len(genotypes2))
                p3 = (1.0000 * genotypes3.count(ref_allele)) / float(len(genotypes3))
                p4 = (1.0000 * genotypes4.count(ref_allele)) / float(len(genotypes4))
                if options.fixed1:
                    if p1 not in [0.0]: continue
                    if p2 in [0.0]: continue
                    if p3 == p4: continue
                    if options.verbose:
                        print(chromosome, position, ref_allele, p1, p2, p3, p4, genotypes1,
                              ''.join(genotypes2),
                              ''.join(genotypes3), ''.join(genotypes4), (p1 - p2) * (p3 - p4))

                if options.fixed3:
                    if p3 not in [1.0, 0.0]: continue
                    if p4 in [1.0, 0.0]: continue
                    if options.verbose:
                        print(p1, p2, p3, p4, ''.join(genotypes1), ''.join(genotypes2), ''.join(genotypes3),
                              ''.join(genotypes4), (p1 - p2) * (p3 - p4))

                if options.fixed4:
                    if p4 not in [0.0]: continue
                    if options.verbose:
                        print(p1, p2, p3, p4, ''.join(genotypes1), ''.join(genotypes2), ''.join(genotypes3),
                              ''.join(genotypes4), (p1 - p2) * (p3 - p4))

                # D-test Patterson et al. 2012
                t = (p1 - p2) * (p3 - p4)
                n = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))

                if options.f3:
                    n = 1.0
                    nW = float(len(genotypes3))
                    p3n = float(len(genotypes3))
                    p3count = genotypes3.count(ref_allele)

                    t = (p3 - p1) * (p3 - p2)  #
                    n = (p3count * (p3n - p3count)) / (p3n * (p3n - 1))
                    n = n / p3n
                    t = t - n
                    B = 2.0 * p3 * (1.0 - p3)
                    n = B

                elif options.f3vanilla:
                    n = 1.0
                    t = (p3 - p1) * (p3 - p2)  #

                elif options.f4:
                    n = 1.0

                t_locus.append(t)
                n_locus.append(n)
                if options.fixed1:
                    t_list.append(t)
                    n_list.append(n)
                    choice_list.append((chromosome, position, t, n))
                    continue

            t = sum(t_locus)  # /len(t_locus)
            n = sum(n_locus)  # /len(n_locus)

            t_list.append(t)
            n_list.append(n)
            # print genotypes1,genotypes2,genotypes3,genotypes4,t,n,p1,p2,p3,p4
            choice_list.append((chromosome, position, t, n))
            if options.verbose and options.fixed1 == False and options.fixed3 == False and options.fixed4 == False:
                print(''.join(genotypes1), ''.join(genotypes2), ''.join(genotypes3), ''.join(genotypes4), t,
                      n, p1, p2, p3, p4)
            continue

        ####
        "LD"
        ####

        if options.LD:
            upperdist = options.LD + options.LDwindow
            lowerdist = options.LD - options.LDwindow

            genotypes1 = selectionhaplotypes(genotypes, targetpop1, chromosome, position)
            genotypes2 = selectionhaplotypes(genotypes, targetpop2, chromosome, position)
            genotypes3 = selectionhaplotypes(genotypes, targetpop3, chromosome, position)
            genotypes4 = selectionhaplotypes(genotypes, targetpop4, chromosome, position)
            if len(genotypes1) < 1: continue
            if len(genotypes2) < 1: continue
            if len(genotypes3) < 1: continue
            if len(genotypes4) < 1: continue

            if options.polymorphic:
                if len(list(set(genotypes1))) != 2: continue
                if len(list(set(genotypes2))) != 2: continue
                if len(list(set(genotypes3))) != 2: continue
                if len(list(set(genotypes4))) != 2: continue

            desiredcount = min([targetpop3count, targetpop4count])
            if targetpop3count != targetpop4count and options.equaln:
                genotypes3 = genotypes3[0:desiredcount]
                genotypes4 = genotypes4[0:desiredcount]

            if options.nomissing:
                if '0' in genotypes1 + genotypes2 + genotypes3 + genotypes4:
                    continue

            if options.scramble:
                genotypes1 = scramblefun(genotypes1)
                genotypes2 = scramblefun(genotypes2)
                genotypes3 = scramblefun(genotypes3)
                genotypes4 = scramblefun(genotypes4)

            if options.popscramble:
                random.shuffle(genotypes1)
                random.shuffle(genotypes2)
                random.shuffle(genotypes3)
                random.shuffle(genotypes4)
            if options.haploidize:
                genotypes1 = haploidizefun(genotypes1)
                genotypes2 = haploidizefun(genotypes2)
                genotypes3 = haploidizefun(genotypes3)
                genotypes4 = haploidizefun(genotypes4)

            siteABBA = False
            siteBABA = False
            if options.autocorr or options.onlyBABA or options.onlyABBA or options.ABBAorBABA:
                tgenotypes1 = genotypes1[0]
                tgenotypes2 = genotypes2[0]
                tgenotypes3 = genotypes3[0]
                tgenotypes4 = genotypes4[0]
                a, b, c, d = tgenotypes1[0], tgenotypes2[0], tgenotypes3[0], tgenotypes4[0]
                if '0' in [a, b, c, d]: continue

                if a == d and b == c and a != b:
                    siteABBA = True
                elif a == c and b == d and a != b:
                    siteBABA = True

                if options.onlyBABA:
                    if siteBABA == False:
                        continue
                if options.onlyABBA:
                    if siteABBA == False:
                        continue
                if options.ABBAorBABA:
                    if siteABBA == False and siteBABA == False:
                        continue

            genotypes1 = ''.join(genotypes1)
            genotypes2 = ''.join(genotypes2)
            genotypes3 = ''.join(genotypes3)
            genotypes4 = ''.join(genotypes4)

            if chromosome != previouschromosome:
                previouschromosome = chromosome
                previouslines = [[position, genotypes1, genotypes2, genotypes3, genotypes4]]
                continue

            if len(previouslines) < 1:
                previouslines = [[position, genotypes1, genotypes2, genotypes3, genotypes4]]
                continue

            newpreviouslinesstart = 0
            for lx in range(0, len(previouslines)):
                l = previouslines[lx]
                thepos = l[0]
                distance = position - thepos
                if distance > upperdist:
                    newpreviouslinesstart = lx
                    continue
                elif distance <= upperdist and distance >= lowerdist:
                    alleles = list(set(genotypes1 + genotypes2 + genotypes3 + genotypes4))
                    alleles = [a for a in alleles if '0' not in a]
                    if len(alleles) > 2: continue

                    lgenotypes1 = l[1]
                    lgenotypes2 = l[2]
                    lgenotypes3 = l[3]
                    lgenotypes4 = l[4]

                    if targetpop3count != targetpop4count and options.equaln:
                        lgenotypes3 = lgenotypes3[0:desiredcount]
                        lgenotypes4 = lgenotypes4[0:desiredcount]

                    lalleles = list(set(lgenotypes1 + lgenotypes2 + lgenotypes3 + lgenotypes4))
                    # alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
                    lalleles = [a for a in lalleles if '0' not in a]
                    if len(lalleles) > 2: continue
                    if len(lgenotypes1) < 1: continue
                    if len(lgenotypes2) < 1: continue
                    if len(lgenotypes3) < 1: continue
                    if len(lgenotypes4) < 1: continue

                    if options.polymorphic:
                        if len(list(set(lgenotypes1))) != 2: continue
                        if len(list(set(lgenotypes2))) != 2: continue
                        if len(list(set(lgenotypes3))) != 2: continue
                        if len(list(set(lgenotypes4))) != 2: continue

                    if options.onlyBABA or options.onlyABBA or options.ABBAorBABA:
                        lgenotypes1 = lgenotypes1[0]
                        lgenotypes2 = lgenotypes2[0]
                        lgenotypes3 = lgenotypes3[0]
                        lgenotypes4 = lgenotypes4[0]
                        if '0' in lgenotypes1 + lgenotypes2 + lgenotypes3 + lgenotypes4: continue
                        ref_allele = lalleles[0]
                        p1 = (1.0000 * lgenotypes1.count(ref_allele)) / float(len(lgenotypes1))
                        p2 = (1.0000 * lgenotypes2.count(ref_allele)) / float(len(lgenotypes2))
                        p3 = (1.0000 * lgenotypes3.count(ref_allele)) / float(len(lgenotypes3))
                        p4 = (1.0000 * lgenotypes4.count(ref_allele)) / float(len(lgenotypes4))
                        # D-test Patterson et al. 2012
                        t = (p1 - p2) * (p3 - p4)
                        n = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))

                        if t == 0.0:
                            continue
                        if siteABBA:
                            t = t - 1.0
                        elif siteBABA:
                            t = t + 1.0
                        choice_list.append((chromosome, position, t, n))
                        t_list.append(t)
                        n_list.append(n)
                        if options.verbose:
                            print(t, '\t\t\t\t', position, distance, tgenotypes1, tgenotypes2, tgenotypes3,
                                  tgenotypes4, t_list)
                            print(n, '\t\t\t\t', position, distance, lgenotypes1, lgenotypes2, lgenotypes3,
                                  lgenotypes4)
                            print('---')
                        continue
                    if options.autocorr2:
                        for ai in range(0, len(genotypes1)):
                            for bi in range(0, len(genotypes2)):
                                for ci in range(0, len(genotypes3)):
                                    for di in range(0, len(genotypes4)):

                                        tgenotypes1 = genotypes1[ai]
                                        tgenotypes2 = genotypes2[bi]
                                        tgenotypes3 = genotypes3[ci]
                                        tgenotypes4 = genotypes4[di]

                                        tlgenotypes1 = lgenotypes1[ai]
                                        tlgenotypes2 = lgenotypes2[bi]
                                        tlgenotypes3 = lgenotypes3[ci]
                                        tlgenotypes4 = lgenotypes4[di]
                                        a, b, c, d = tgenotypes1[0], tgenotypes2[0], tgenotypes3[0], \
                                            tgenotypes4[
                                                0]
                                        if '0' in [a, b, c, d]: continue
                                        # print a,b,c,d
                                        siteABBA = False
                                        siteBABA = False
                                        if a == d and b == c and a != b:
                                            siteABBA = True
                                        elif a == c and b == d and a != b:
                                            siteBABA = True
                                        if '0' in lgenotypes1 + lgenotypes2 + lgenotypes3 + lgenotypes4: continue
                                        ref_allele = tlgenotypes1[0]
                                        p1 = (1.0000 * tlgenotypes1.count(ref_allele)) / float(
                                            len(tlgenotypes1))
                                        p2 = (1.0000 * tlgenotypes2.count(ref_allele)) / float(
                                            len(tlgenotypes2))
                                        p3 = (1.0000 * tlgenotypes3.count(ref_allele)) / float(
                                            len(tlgenotypes3))
                                        p4 = (1.0000 * tlgenotypes4.count(ref_allele)) / float(
                                            len(tlgenotypes4))
                                        # D-test Patterson et al. 2012
                                        t = (p1 - p2) * (p3 - p4)
                                        n = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))

                                        if t == 0.0:
                                            continue
                                        # print t
                                        if siteABBA == True:
                                            t = t - 1.0
                                        elif siteBABA == True:
                                            t = t + 1.0
                                        choice_list.append((chromosome, position, t, n))
                                        t_list.append(t)
                                        n_list.append(n)
                                        if options.verbose:
                                            print(t, '\t\t\t\t', position, distance, tgenotypes1,
                                                  tgenotypes2,
                                                  tgenotypes3, tgenotypes4, t_list)
                                            print(n, '\t\t\t\t', position, distance, lgenotypes1,
                                                  lgenotypes2,
                                                  lgenotypes3, lgenotypes4)
                                            print('---')
                                        continue
                    elif options.autocorr:
                        """
                        the autocorrelation statistic is
                        """
                        lgenotypes1 = random.choice(lgenotypes1)
                        lgenotypes2 = random.choice(lgenotypes2)
                        lgenotypes3 = random.choice(lgenotypes3)
                        lgenotypes4 = random.choice(lgenotypes4)
                        a, b, c, d = genotypes1[0], genotypes2[0], genotypes3[0], genotypes4[0]
                        x, y, z, v = lgenotypes1[0], lgenotypes2[0], lgenotypes3[0], lgenotypes4[0]
                        if '0' in [a, b, c, d, x, y, z, v]: continue
                        state = '0'
                        # abba
                        if a == d and b == c and a != b:
                            if x == v and y == z and x != y:
                                state = 'dABBA'
                                abbalist.append(state)
                                if options.verbose:
                                    print(position)
                                    print(state)
                                    print(a, b, c, d)
                                    print(x, y, z, v)
                            elif x == z and y == v and x != y:
                                state = 'ndABBA'
                                abbalist.append(state)
                                if options.verbose:
                                    print(position)
                                    print(state)
                                    print(a, b, c, d)
                                    print(x, y, z, v)
                        elif a == c and b == d and a != b:
                            if x == z and y == v and x != y:
                                state = 'dBABA'
                                babalist.append(state)
                                if options.verbose:
                                    print(position)
                                    print(state)
                                    print(a, b, c, d)
                                    print(x, y, z, v)
                            elif x == v and y == z and x != y:
                                state = 'ndBABA'
                                babalist.append(state)
                                if options.verbose:
                                    print(position)
                                    print(genotypes1, genotypes2, genotypes3, genotypes4)
                                    print(lgenotypes1, lgenotypes2, lgenotypes3, lgenotypes4)
                                    print(state)
                                    print(a, b, c, d)
                                    print(x, y, z, v)

                        continue

                    pfreqs = []
                    qfreqs = []
                    pqfreqs = []
                    pqfreqs2 = []
                    pqfreqs3 = []
                    pcounts = []
                    qcounts = []
                    for pop in [[genotypes1, lgenotypes1], [genotypes2, lgenotypes2],
                                [genotypes3, lgenotypes3],
                                [genotypes4, lgenotypes4]]:
                        haps = []
                        nomissing = ''
                        lnomissing = ''
                        for a, b in zip(pop[0], pop[1]):
                            if '0' in a or '0' in b: continue
                            nomissing += a
                            lnomissing += b
                            haps.append(a + b)
                        # print a+b

                        if len(nomissing) < 2 or len(lnomissing) < 2: continue

                        pqfreq = 1.0 * haps.count(alleles[0] + lalleles[0]) / len(haps)
                        pqfreqs.append(pqfreq)

                        if options.hapD:
                            allhaps = [alleles[0] + lalleles[0]]
                            if len(lalleles) > 1:
                                allhaps.append(alleles[0] + lalleles[1])
                                pqfreq2 = 1.0 * haps.count(alleles[0] + lalleles[1]) / len(haps)
                                pqfreqs2.append(pqfreq2)
                            else:
                                pqfreqs2.append(0.0)
                            if len(alleles) > 1:
                                allhaps.append(alleles[1] + lalleles[0])
                                pqfreq3 = 1.0 * haps.count(alleles[1] + lalleles[0]) / len(haps)
                                pqfreqs3.append(pqfreq3)
                            else:
                                pqfreqs3.append(0.0)

                        pfreqs.append(1.0 * nomissing.count(alleles[0]) / len(nomissing) * 1.0)
                        qfreqs.append(1.0 * lnomissing.count(lalleles[0]) / len(lnomissing) * 1.0)
                        pcounts.append(len(nomissing) * 1.0)
                        qcounts.append(len(lnomissing) * 1.0)
                    if len(pcounts) != 4: continue
                    # print pfreqs,qfreqs
                    # print pcounts,qcounts
                    # print '--'
                    #### for when allele frequency for P1,P2,P3,P4 is computed in P1,P2 and P3,P4 instead of each separately
                    #### haplotype frequency is still computed separately
                    if options.f3 == False and options.withinfreq == False:
                        newpfreqs = []
                        newqfreqs = []
                        stat = (pfreqs[0] * pcounts[0] + pfreqs[1] * pcounts[1]) / (pcounts[0] + pcounts[1])
                        newpfreqs.append(stat)
                        newpfreqs.append(stat)

                        stat = (pfreqs[2] * pcounts[2] + pfreqs[3] * pcounts[3]) / (pcounts[2] + pcounts[3])
                        newpfreqs.append(stat)
                        newpfreqs.append(stat)

                        stat = (qfreqs[0] * qcounts[0] + qfreqs[1] * qcounts[1]) / (qcounts[0] + qcounts[1])
                        newqfreqs.append(stat)
                        newqfreqs.append(stat)

                        stat = (qfreqs[2] * qcounts[2] + qfreqs[3] * qcounts[3]) / (qcounts[2] + qcounts[3])
                        newqfreqs.append(stat)
                        newqfreqs.append(stat)
                        """
                        print pfreqs,qfreqs
                        print pcounts,qcounts
                        print newpfreqs,newqfreqs
                        print '--'
                        """
                        qfreqs = newqfreqs
                        if options.withinoutgroupsfreq == False:
                            pfreqs = newpfreqs

                    # print alleles[0],lalleles[0]
                    # print pqfreqs,pfreqs,qfreqs
                    if len(pqfreqs) < 4:
                        continue
                    Ds = []
                    for p, q, pq in zip(pfreqs, qfreqs, pqfreqs):
                        Dstat = 1.0 * pq - (p * q)
                        if options.Dprim:
                            if Dstat < 0.0:
                                Dmax = min([1.0 * p * q, (1.0 - p) * (1.0 - q)])
                            elif Dstat > 0.0:
                                Dmax = min([1.0 * p * (1.0 - q), (1.0 - p) * q])
                            elif Dstat != 0.0:
                                Dmax = 1.0
                            if Dstat != float(0.0):
                                # print Dstat,Dmax
                                Dstat = Dstat / Dmax

                        elif options.Dr:
                            # print Dstat,p,q
                            Dstat = Dstat / math.sqrt(p * (1.0 - p) * q * (1.0 - q))
                        Ds.append(Dstat)
                    # print alleles[0]+lalleles[0],pq,p,q,Dstat
                    # print Ds
                    if len(Ds) < 1:
                        continue
                    # print Ds,Ds.count(0.0)
                    # if 0.0 in Ds[0:2] and 0.0 in Ds[2:]:continue
                    # if 0.0 in Ds:continue
                    # if 0.0 in Ds[2:]:continue
                    try:
                        D4 = (Ds[0] - Ds[1]) * (Ds[2] - Ds[3])
                    except IndexError:
                        print(Ds)
                        print(genotypes1, genotypes2, genotypes3, genotypes4)
                        print(lgenotypes1, lgenotypes2, lgenotypes3, lgenotypes4)
                        print('---')
                    t = D4  # /distance

                    if options.countDzero:
                        if Ds[1] == 0.0 and Ds[0] != 0.0:
                            t = 1.0
                        else:
                            t = 0.0
                        D4 = t

                    if options.doubletest:
                        f4_1 = (pfreqs[0] - pfreqs[1]) * (pfreqs[2] - pfreqs[3])
                        f4_2 = (qfreqs[0] - qfreqs[1]) * (qfreqs[2] - qfreqs[3])
                        # f4=sum([f4_1,f4_2])/2.0
                        t = D4 + f4_1  # +f4_2

                    n = 1.0
                    if options.f3:
                        D3 = (Ds[0] - Ds[3]) * (Ds[1] - Ds[3])
                        t = D3
                        n = 1.0
                    elif options.Dcorr:
                        D4 = Ds[0] * Ds[1]
                        t = D4
                        n = math.sqrt((Ds[0] * Ds[0]) * (Ds[1] * Ds[1]))


                    elif options.hapD:
                        t_locus = []
                        n_locus = []
                        if options.twohaps:
                            if len(allhaps) != 2: continue
                        for hapallele in [pqfreqs, pqfreqs2, pqfreqs3]:
                            p1 = hapallele[0]
                            p2 = hapallele[1]
                            p3 = hapallele[2]
                            p4 = hapallele[3]
                            t = (p1 - p2) * (p3 - p4)
                            n = 1.0
                            n = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))
                            Ds = [p1, p2, p3, p4]
                            D4 = t
                            t_locus.append(t)
                            n_locus.append(n)

                        t = sum(t_locus)  # /len(t_locus)
                        n = sum(n_locus)  # /len(n_locus)

                    # if t==0.0 or t==-0.0 and options.hapD ==False:continue
                    if options.verbose:
                        print(Ds[0], Ds[1], Ds[2], Ds[3], D4, '\t\t\t\t', position, alleles[0], genotypes1,
                              genotypes2, genotypes3, genotypes4)
                        print(Ds[0], Ds[1], Ds[2], Ds[3], D4, '\t\t\t\t', position, lalleles[0], lgenotypes1,
                              lgenotypes2, lgenotypes3, lgenotypes4)
                        print('---')
                    choice_list.append((chromosome, position, t, n))
                    t_list.append(t)
                    n_list.append(n)
                # break
                else:
                    continue

            previouslines.append([position, genotypes1, genotypes2, genotypes3, genotypes4])
            previouslines = previouslines[newpreviouslinesstart:]
            # choice_list.append((chromosome,position,t,n))
            continue

        genotypes1 = selectiongenotypes(genotypes, targetpop1, chromosome, position)
        genotypes2 = selectiongenotypes(genotypes, targetpop2, chromosome, position)
        genotypes3 = selectiongenotypes(genotypes, targetpop3, chromosome, position)
        genotypes4 = selectiongenotypes(genotypes, targetpop4, chromosome, position)

        if options.haploidize:
            genotypes1 = haploidizefun(genotypes1)
            genotypes2 = haploidizefun(genotypes2)
            genotypes3 = haploidizefun(genotypes3)
            genotypes4 = haploidizefun(genotypes4)

        if len(genotypes1) < options.mincount: continue
        if len(genotypes2) < options.mincount: continue
        if len(genotypes3) < options.mincount: continue
        if len(genotypes4) < options.mincount: continue

        if options.pops2:
            genotypes21 = selectiongenotypes(genotypes, targetpop21, chromosome, position)
            genotypes22 = selectiongenotypes(genotypes, targetpop22, chromosome, position)
            genotypes23 = selectiongenotypes(genotypes, targetpop23, chromosome, position)
            genotypes24 = selectiongenotypes(genotypes, targetpop24, chromosome, position)

            if len(genotypes21) < 1: continue
            if len(genotypes22) < 1: continue
            if len(genotypes23) < 1: continue
            if len(genotypes24) < 1: continue

        if options.fixed1:
            if len(list(set(genotypes1))) != 1:
                continue
            allele1 = genotypes1[0]
        # allele2=[a for a in genotypes2 if a != allele1]
        # if len(allele2) <1:continue
        # print genotypes2
        # genotypes2=allele2[0]
        # if len(list(set(genotypes3+genotypes4))) ==1:continue

        if options.allhaploid:
            genotypes1 = random.choice(genotypes1)
            genotypes2 = random.choice(genotypes2)
            genotypes3 = random.choice(genotypes3)
            genotypes4 = random.choice(genotypes4)
            if options.informative:
                if genotypes3 == genotypes4:
                    continue
                if genotypes1 == genotypes2:
                    continue

        if options.haploidinclade:
            genotypes3 = random.choice(genotypes3)
            genotypes4 = random.choice(genotypes4)
            if genotypes3 == genotypes4:
                continue

        if options.nomissing:
            if len(genotypes1) != targetpop1count: continue
            if len(genotypes2) != targetpop2count: continue
            if len(genotypes3) != targetpop3count: continue
            if len(genotypes4) != targetpop4count: continue

        if options.ancestor:
            ancgenotypes = selectiongenotypes(genotypes, ancpop, chromosome, position)
            if len(ancgenotypes) < 1 or len(list(set(ancgenotypes))) != 1: continue

        if options.testpop:
            testgenotypes = selectiongenotypes(genotypes, testpop, chromosome, position)
            if len(testgenotypes) < 1: continue
        # print genotypes1,genotypes2,genotypes3,genotypes4,testgenotypes

        if options.SNPfreq:
            dercount = genotypes1.count(ancgenotypes[0])
            totn = len(genotypes1)
            derfreq = 1.0 * dercount / totn
            print(dercount, '\t', derfreq, '\t', totn)
            exit(0)

        if options.ascertainfreq != -1:
            ascgenotypes = selectiongenotypes(genotypes, ascpop, chromosome, position)
            if options.downsampleasc:
                if len(ascgenotypes) < options.downsampleasc: continue
                ascgenotypes = ''.join(random.sample(ascgenotypes, options.downsampleasc))

            if len(ascgenotypes) < 1: continue
            derivedfreq = ascpopcount - ascgenotypes.count(ancgenotypes[0])
            if len(ascgenotypes) != ascpopcount or derivedfreq != options.ascertainfreq: continue
        elif options.ascertain:
            ascgenotypes = selectiongenotypes(genotypes, ascpop, chromosome, position)
            if len(ascgenotypes) < 1 or len(list(set(ascgenotypes))) != 2: continue

        if options.ascertainfreq2 != -1:
            ascgenotypes2 = selectiongenotypes(genotypes, ascpop2, chromosome, position)
            if options.downsampleasc2:
                if len(ascgenotypes2) < options.downsampleasc2: continue
                ascgenotypes2 = ''.join(random.sample(ascgenotypes2, options.downsampleasc2))

            if len(ascgenotypes2) < 1: continue
            derivedfreq2 = ascpopcount2 - ascgenotypes2.count(ancgenotypes[0])
            if len(ascgenotypes2) != ascpopcount2 or derivedfreq2 != options.ascertainfreq2: continue
        elif options.ascertain2:
            ascgenotypes2 = selectiongenotypes(genotypes, ascpop2, chromosome, position)
            if len(ascgenotypes2) < 1 or len(list(set(ascgenotypes2))) != 2: continue

        if options.ascertainfreq3 != -1:
            ascgenotypes3 = selectiongenotypes(genotypes, ascpop3, chromosome, position)
            if options.downsampleasc3:
                if len(ascgenotypes3) < options.downsampleasc3: continue
                ascgenotypes3 = ''.join(random.sample(ascgenotypes3, options.downsampleasc3))

            if len(ascgenotypes3) < 1: continue
            derivedfreq3 = ascpopcount3 - ascgenotypes3.count(ancgenotypes[0])
            if len(ascgenotypes3) != ascpopcount3 or derivedfreq3 != options.ascertainfreq3: continue
        elif options.ascertain3:
            ascgenotypes3 = selectiongenotypes(genotypes, ascpop3, chromosome, position)
            if len(ascgenotypes3) < 1 or len(list(set(ascgenotypes3))) != 2: continue

        if options.pop2freq:
            if len(genotypes2) != targetpop2count: continue
            obspop2freq = targetpop2count - genotypes2.count(ancgenotypes[0])
            if obspop2freq != options.pop2freq: continue

        alleles = list(set(genotypes))
        # alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
        alleles = [a for a in alleles if a != '0']

        if options.mutatepop2:
            if len(alleles) != 2: continue
            genotypes2 = mutate(genotypes2, options.mutatepop2, alleles)

        if options.mutatepop4:
            if len(alleles) != 2: continue
            genotypes4 = mutate(genotypes4, options.mutatepop4, alleles)

        alleles = list(set(genotypes1 + genotypes2 + genotypes3 + genotypes4))
        if options.testpop:
            alleles = list(set(genotypes1 + genotypes2 + genotypes3 + genotypes4 + testgenotypes))
        if len(alleles) > 2: continue
        if options.polymorphic:
            if len(alleles) != 2:
                continue
        if options.mutationclass:
            badclass = True
            # print mutationclass[0],mutationclass[1]
            if mutationclass[0] in alleles and mutationclass[1] in alleles:
                badclass = False
            if badclass == True:
                continue
        if options.informative:
            alleles1 = list(set(genotypes1 + genotypes2))
            alleles2 = list(set(genotypes3 + genotypes4))
            if len(alleles1) != 2: continue
            if len(alleles2) != 2: continue

            if options.testpop:
                alleles1 = list(set(genotypes1 + genotypes2))
                alleles2 = list(set(testgenotypes + genotypes4))
                if len(alleles1) != 2: continue
                if len(alleles2) != 2: continue

        if options.ancestor:
            ref_allele = random.choice(ancgenotypes)
        else:
            try:
                ref_allele = random.choice(alleles)
            except IndexError:
                print(alleles)

        if notransitions:
            if 'C' in genotypes and 'T' in genotypes:
                continue
            elif 'G' in genotypes and 'A' in genotypes:
                continue

        if maf:
            # print genotypes1,genotypes2,genotypes3,genotypes4

            af1 = 1.000 * genotypes1.count(genotypes1[0]) / len(genotypes1)
            af2 = 1.000 * genotypes2.count(genotypes2[0]) / len(genotypes2)
            af3 = 1.000 * genotypes3.count(genotypes3[0]) / len(genotypes3)
            af4 = 1.000 * genotypes4.count(genotypes4[0]) / len(genotypes4)

            maf1 = min([af1, (1.00 - af1)])
            maf2 = min([af2, (1.00 - af2)])
            maf3 = min([af3, (1.00 - af3)])
            maf4 = min([af4, (1.00 - af4)])
            maf_list.append((maf1, maf2, maf3, maf4))
            # print maf1,maf2,maf3,maf4
            if (maf1 < threshold) and (len(genotypes1) > minsamplesize): continue
            if (maf2 < threshold) and (len(genotypes2) > minsamplesize): continue
            if (maf3 < threshold) and (len(genotypes3) > minsamplesize): continue
            if (maf4 < threshold) and (len(genotypes4) > minsamplesize): continue

        p1 = (1.0000 * genotypes1.count(ref_allele)) / float(len(genotypes1))
        p2 = (1.0000 * genotypes2.count(ref_allele)) / float(len(genotypes2))
        p3 = (1.0000 * genotypes3.count(ref_allele)) / float(len(genotypes3))
        p4 = (1.0000 * genotypes4.count(ref_allele)) / float(len(genotypes4))


        if options.outdiff:
            outdiff = abs(p1 - p2)
            if outdiff < options.outdiff: continue
        if options.outdiffexact:
            outdiff = abs(p1 - p2)
            if str(outdiff)[0:3] != str(options.outdiffexact)[0:3]: continue

        if options.indiff:
            indiff = abs(p3 - p4)
            if indiff < options.indiff: continue

        # D-test Patterson et al. 2012
        t = (p1 - p2) * (p3 - p4)
        n = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))

        if options.pops2:
            p1 = (1.0000 * genotypes21.count(ref_allele)) / float(len(genotypes21))
            p2 = (1.0000 * genotypes22.count(ref_allele)) / float(len(genotypes22))
            p3 = (1.0000 * genotypes23.count(ref_allele)) / float(len(genotypes23))
            p4 = (1.0000 * genotypes24.count(ref_allele)) / float(len(genotypes24))
            t2 = (p1 - p2) * (p3 - p4)
            n2 = (p1 + p2 - (2.0 * p1 * p2)) * (p3 + p4 - (2.0 * p3 * p4))
            t = t
            n = t2  # 1.0#+n2

        if options.pop2weight:
            # t=t*((1.0-p2)*0.06-0.02)
            t = (math.exp(p1 - p2)) * (p3 - p4)

        if options.f4:
            n = 1.0
        elif options.f5:
            if options.ancestor:
                outgrfreq = 1.0
            if options.testpop:
                p5 = (1.0000 * testgenotypes.count(ref_allele)) / float(len(testgenotypes))
                outgrfreq = p5

            f3weight = (p2 - p5) * (p2 - p1)
            t = t * f3weight
            n = 1.0

        elif options.simpleD:
            if genotypes3[0] == ancgenotypes[0]:
                t = 1.0
            elif genotypes4[0] == ancgenotypes[0]:
                t = -1.0
            n = 1.0

        elif options.symmetry:
            choicepop1 = random.choice(genotypes1)
            choicepop2 = random.choice(genotypes2)
            ancchoice = random.choice(ancgenotypes)
            if choicepop1 == choicepop2:
                continue
            elif choicepop1 != ancchoice:
                t = -1.0
                n = 1.0
            elif choicepop2 != ancchoice:
                t = 1.0
                n = 1.0
            else:
                continue

        elif options.simpleDfreq:
            freq3 = (1.0 * len(genotypes3) - genotypes3.count(ancgenotypes[0])) / (1.0 * len(genotypes3))
            freq4 = (1.0 * len(genotypes4) - genotypes4.count(ancgenotypes[0])) / (1.0 * len(genotypes4))

            expBABA = float(freq3) * (1.0 - float(freq4))
            expABBA = (1.0 - float(freq3)) * float(freq4)
            t = expABBA - expBABA
            n = 1.0

        elif options.simpleDtailtest:
            freq3 = (1.0 * len(genotypes3) - genotypes3.count(ancgenotypes[0])) / (1.0 * len(genotypes3))
            freq4 = (1.0 * len(genotypes4) - genotypes4.count(ancgenotypes[0])) / (1.0 * len(genotypes4))

            expBABA = float(freq3) * (1.0 - float(freq4))
            expABBA = (1.0 - float(freq3)) * float(freq4)
            t = expABBA - expBABA
            if p2 == 0.0:
                t = t * 1.0
            elif p2 == 1.0:
                t == t * -1.0
            else:
                continue
            n = 1.0

        elif options.testpop and options.f5 == False:
            pt = (1.0000 * testgenotypes.count(ref_allele)) / float(len(testgenotypes))
            t = (p1 - p2) * (pt - p4)
            n = (p1 - p2) * (p3 - p4)

        elif options.LiReich:
            pop1ind = random.sample(genotypes1, 2)
            alleles = list(set(pop1ind))
            if len(alleles) != 2:
                continue
            else:
                pop2copy = random.sample(genotypes2, 1)
                pop2copy = pop2copy[0]
                if pop2copy != ancgenotypes[0]:
                    t = 1.0
                elif pop2copy == ancgenotypes[0]:
                    t = 0.0
                # print ancgenotypes[0],pop1ind,pop2copy,t
                n = 1.0

        elif options.FSTWC or options.Tdiv:
            FST_res = FST_W_pairwise(
                [len(genotypes1), genotypes1.count(ref_allele), len(genotypes2),
                 genotypes2.count(ref_allele)])
            t = FST_res[0]
            n = FST_res[1]

        elif options.FST:
            FST_res = FST_H_pairwise(
                [len(genotypes1), genotypes1.count(ref_allele), len(genotypes2),
                 genotypes2.count(ref_allele)])
            t = FST_res[0]
            n = FST_res[1]

        elif options.f3:
            n = 1.0
            nW = float(len(genotypes3))
            # print p1,p2,p3

            p3n = float(len(genotypes3))
            p3count = genotypes3.count(ref_allele)

            t = (p3 - p1) * (p3 - p2)  #
            # n= 2.0*( (p3*(1.0-p3) ))
            n = (p3count * (p3n - p3count)) / (p3n * (p3n - 1))
            n = n / p3n
            t = t - n
            if options.nohzcorrection == False:
                B = 2.0 * p3 * (1.0 - p3)
                n = B
            else:
                n = 1.0

        elif options.f3vanilla:
            n = 1.0
            nW = float(len(genotypes3))
            # print p1,p2,p3

            p3n = float(len(genotypes3))
            p3count = genotypes3.count(ref_allele)

            t = (p3 - p1) * (p3 - p2)
            n = 1.0

        elif options.f2:
            t = (p1 - p2) ** 2
            n = 1.0

        elif options.mSFS:
            if len(genotypes1) != targetpop1count: continue
            if len(genotypes2) != targetpop2count: continue
            ancestor = list(set(ancgenotypes))[0]
            targetcount1 = len(genotypes1) - genotypes1.count(ancestor)
            targetcount2 = len(genotypes2) - genotypes2.count(ancestor)

            r1 = random.choice(genotypes3)
            r2 = random.choice(genotypes4)
            if r2 == ancestor and r1 != ancestor:
                t_list.append((targetcount1, targetcount2))
            elif r2 != ancestor and r1 == ancestor:
                n_list.append((targetcount1, targetcount2))
            continue

        elif options.SFS:
            if len(genotypes2) != targetpop2count: continue
            ancestor = list(set(genotypes1))[0]
            targetcount = len(genotypes2) - genotypes2.count(ancestor)

            r1 = random.choice(genotypes3)
            r2 = random.choice(genotypes4)
            if r2 == ancestor and r1 != ancestor:
                t_list.append(targetcount)
            elif r2 != ancestor and r1 == ancestor:
                n_list.append(targetcount)
            continue

        elif options.SFS2:
            if len(genotypes2) != targetpop2count: continue
            ancestor = list(set(genotypes1))[0]
            targetcount = len(genotypes2) - genotypes2.count(ancestor)
            samplechoice = min(len(genotypes3), len(genotypes4))
            r1 = 10 - random.sample(genotypes3, samplechoice).count(ancestor)
            r2 = 10 - random.sample(genotypes4, samplechoice).count(ancestor)

            if r1 == 1 and r2 == 0:
                t_list.append(targetcount)
            elif r1 == 0 and r2 == 1:
                n_list.append(targetcount)
            continue

        elif options.anc_test:
            if len(genotypes2) != targetpop2count: continue
            ancestor = list(set(genotypes1))[0]
            if genotypes2.count(ancestor) != len(genotypes2): continue

            t_temp = []
            then = 100
            abba = 0
            baba = 0
            for num in range(0, then):
                r1 = random.choice(genotypes3)
                r2 = random.choice(genotypes4)
                if r2 == ancestor and r1 != ancestor:
                    abba += 1.0
                elif r2 != ancestor and r1 == ancestor:
                    baba += 1.0

            t = (abba - baba)
            n = 1.0

        hz1 = p1 * (1.0 - p1)
        hz2 = p2 * (1.0 - p2)
        hz3 = p3 * (1.0 - p3)
        hz4 = p4 * (1.0 - p4)
        avhz = (hz1 + hz2 + hz3 + hz4) / 4.0
        t_list.append(t)
        n_list.append(n)
        if options.verbose:
            print(chromosome, position, genotypes1, genotypes2, genotypes3, genotypes4, p1, p2, p3, p4, t,
                  hz1,
                  hz2, hz3, hz4)

        if options.nojackknife: continue
        if options.Bcorr or options.Bstrat:
            choice_list.append((chromosome, position, t, n, Bstat))
            # b_list.append(Bstat)
            continue
        choice_list.append((chromosome, position, t, n))  # ,avhz

    choice_list.append((999, 1, 0, 0, 0))
    return choice_list, t_list, n_list, abbalist, babalist
