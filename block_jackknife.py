# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 18:17
@Author : Shuya Zhang, Ngenie Liang
@File : block_jackknife.py
@Project : popstats
"""

from computing_functions import Bstratfun, Bcorrfun


def make_dp_main(options, poplabel1, poplabel2, poplabel3, poplabel4, t_list, n_list, choice_list):
    """
    :param options:
    :param poplabel1:
    :param poplabel2:
    :param poplabel3:
    :param poplabel4:
    :param t_list:
    :param n_list:
    :param choice_list:
    :return:
    """
    Dp_main = sum(t_list) / sum(n_list)
    Dp_main2 = 0
    blist = []
    # print Dp_main
    num_SNPs = len(t_list)
    if options.positivestat and (Dp_main < 0.0):
        Dp_main = Dp_main * -1.0
        poplabel3temp = poplabel4
        poplabel4 = poplabel3
        poplabel3 = poplabel3temp

    if options.Bstrat:
        Dp_main = Bstratfun(choice_list)
        blist = []
        Dp_main2 = sum(t_list) / sum(n_list)
    if options.Bcorr:
        Bcorrstat, nsum = Bcorrfun(choice_list)
        # print Bcorrstat
        Dp_main = Bcorrstat

        for i in range(0, 10):
            tdat = []
            ndat = []
            for dat in choice_list:
                if int(dat[4]) == i:
                    tdat.append(dat[2])
                    ndat.append(dat[3])
            tdat_b = sum(tdat)
            ndat_b = sum(ndat)
            try:
                datstat = tdat_b / ndat_b
            except ZeroDivisionError:
                datstat = 0.0
            blist.append(datstat)
        # print ''
        if options.nickverbose:
            blist = []
            Dp_main2 = sum(t_list) / sum(n_list)

    if options.nojackknife:
        print('+'.join(poplabel1), '\t', '+'.join(poplabel2), '\t', '+'.join(poplabel3), '\t',
              '+'.join(poplabel4), '\t', float(Dp_main), '\t', 'NA', '\t', 'NA', '\t', 'NA', '\t', num_SNPs)
        exit(0)

    return Dp_main, blist, Dp_main2, num_SNPs


def block_jackknife(options, t_list, n_list, choice_list):
    """
    :param options:
    :param poplabel1:
    :param poplabel2:
    :param poplabel3:
    :param poplabel4:
    :param t_list:
    :param n_list:
    :param choice_list:
    :return:
    """
    t_list_b = []
    n_list_b = []

    current_chromosome = 0

    prev_position = 0
    if options.chromosome == 23 and options.not23 == False:
        prev_position = 2700000

    if options.numSNPblocks > 0:
        block_size = int(1.0 * len(t_list) / options.numSNPblocks)

    sumt_main = sum(t_list)
    sumn_main = sum(n_list)

    if options.excludeoutliers:
        blockDs = {}
        for line in choice_list:
            col = line
            chromosome = int(col[0])
            position = int(col[1])
            location = str(chromosome) + '_' + str(position)
            t = col[4]
            n = col[3]
            t_list_b.append(t)
            n_list_b.append(n)
            if True:  # basepair blocks
                if chromosome != current_chromosome:
                    t_list_b = []
                    n_list_b = []
                    current_chromosome = chromosome
                    prev_position = 0
                    continue

                elif position > (prev_position + options.block_size):
                    sumt = sum(t_list_b)
                    sumn = len(t_list_b)
                    D_pop = sumt / sumn
                    blockDs[D_pop] = location
                    t_list_b = []
                    n_list_b = []
                    prev_position = position

        blockDvals = list(blockDs.keys())
        totblocks = len(blockDvals)
        percentage = 0.05
        toexclude = int(percentage * totblocks)
        sortedDs = sorted(blockDvals)
        print(sortedDs, toexclude, totblocks - toexclude - 1)
        minDs = sortedDs[toexclude]
        maxDs = sortedDs[int(totblocks - toexclude - 1)]
        print(minDs, maxDs)
        bannedlocations = []
        for k in blockDvals:
            if k <= minDs or k >= maxDs:
                bannedlocations.append(blockDs[k])

        current_chromosome = 0

        prev_position = 0
        if options.chromosome == 23 and options.not23 == False:
            prev_position = 2700000
        for line in choice_list:
            col = line
            chromosome = int(col[0])
            position = int(col[1])
            location = str(chromosome) + '_' + str(position)
            t = col[2]
            n = col[3]
            t_list_b.append(t)
            n_list_b.append(n)
            if True:  # basepair blocks
                if chromosome != current_chromosome:
                    t_list_b = []
                    n_list_b = []
                    current_chromosome = chromosome
                    prev_position = 0
                    continue

                elif position > (prev_position + options.block_size):
                    sumt = sum(t_list_b)
                    sumn = sum(n_list_b)
                    D_pop = sumt / sumn
                    if location in bannedlocations:
                        # print sumt_main/sumn_main
                        sumt_main = sumt_main - sumt
                        sumn_main = sumn_main - sumt
                    t_list_b = []
                    n_list_b = []
                    prev_position = position

        Dp_main = sumt_main / sumn_main

    t_list_b = []
    n_list_b = []

    jackknife_Dp = []
    current_chromosome = 1  ### ASSUMES SORTED CHROMOSOMES IN TPED
    mjlist = []
    t_list_bootstrap = []
    n_list_bootstrap = []

    prev_position = 0
    if options.chromosome == 23 and options.not23 == False:
        prev_position = 2700000
    counter = 0
    if options.morgan:
        prev_position = 0.0

    nickcounter = 0
    bcounter = 0
    for line in choice_list:
        # col= line.split()
        col = line
        chromosome = int(col[0])
        position = int(col[1])
        if options.morgan:
            position = float(col[1])

        t = col[2]
        n = col[3]
        bcounter += 1
        if bcounter == 1:
            prev_position = position
            blockstartindex = bcounter
        t_list_b.append(t)
        n_list_b.append(n)

        if options.numSNPblocks:
            if (len(t_list_b)) >= block_size:
                sumt = sumt_main - sum(t_list_b)
                sumn = sumn_main - sum(n_list_b)
                if options.bootstrap != 0:
                    sumt = sum(t_list_b)
                    sumn = sum(n_list_b)
                    t_list_bootstrap.append(sumt)
                    n_list_bootstrap.append(sumn)
                D_pop = sumt / sumn
                jackknife_Dp.append(D_pop)
                t_list_b = []
                n_list_b = []
            continue

        elif options.chromblocks:
            if chromosome != current_chromosome:
                sumt = sumt_main - sum(t_list_b)
                sumn = sumn_main - sum(n_list_b)
                if options.bootstrap != 0:
                    sumt = sum(t_list_b)
                    sumn = sum(n_list_b)
                    t_list_bootstrap.append(sumt)
                    n_list_bootstrap.append(sumn)
                D_pop = sumt / sumn
                jackknife_Dp.append(D_pop)
                mjlist.append(len(t_list_b))
                if options.verboseblocks: print(chromosome, '\t', sum(t_list_b) / sum(n_list_b), '\t',
                                                len(t_list_b), '\t', sum(t_list_b), '\t', sum(n_list_b))

                t_list_b = []
                n_list_b = []
                current_chromosome = chromosome
            continue

        elif options.Bcorr or options.Bstrat:
            if chromosome != current_chromosome:
                choice_list_b = choice_list[0:blockstartindex] + choice_list[bcounter:]

                mjcount = bcounter - blockstartindex
                if options.Bstrat:
                    D_pop = Bstratfun(choice_list_b)
                else:
                    D_pop, nsum = Bcorrfun(choice_list_b)
                    mjcount = nsum
                jackknife_Dp.append(D_pop)
                mjlist.append(mjcount)

                if options.nickverbose:
                    nickcounter += 1
                    tdat = []
                    ndat = []
                    for dat in choice_list_b:
                        tdat.append(dat[2])
                        ndat.append(dat[3])
                    tdat_b = sum(tdat)  # /len(tdat)
                    ndat_b = sum(ndat)  # /len(ndat)
                    print(nickcounter, tdat_b / ndat_b, D_pop, mjcount)  # ,

                blockstartindex = bcounter
                current_chromosome = chromosome
                prev_position = position
                continue
            elif position > (prev_position + options.block_size):
                choice_list_b = choice_list[0:blockstartindex] + choice_list[bcounter:]
                if options.Bstrat:
                    D_pop = Bstratfun(choice_list_b)
                else:
                    D_pop, nsum = Bcorrfun(choice_list_b)
                # print D_pop
                jackknife_Dp.append(D_pop)
                mjcount = bcounter - blockstartindex
                mjlist.append(mjcount)

                if options.nickverbose:
                    nickcounter += 1
                    tdat = []
                    ndat = []
                    for dat in choice_list_b:
                        tdat.append(dat[2])
                        ndat.append(dat[3])
                    tdat_b = sum(tdat)  # /len(tdat)
                    ndat_b = sum(ndat)  # /len(ndat)
                    print(nickcounter, tdat_b / ndat_b, D_pop, mjcount)  # ,

                # print jackknife_Dp
                blockstartindex = bcounter
                prev_position = position

        else:  # basepair blocks
            if chromosome != current_chromosome:

                if True:
                    sumt = sumt_main - sum(t_list_b)
                    sumn = sumn_main - sum(n_list_b)

                    if options.bootstrap != 0:
                        if len(t_list_b) > 1:
                            sumt = sum(t_list_b)
                            sumn = sum(n_list_b)
                            t_list_bootstrap.append(sumt)
                            n_list_bootstrap.append(sumn)

                    mjlist.append(len(t_list_b))

                    D_pop = sumt / sumn
                    jackknife_Dp.append(D_pop)
                    if options.verboseblocks: print(chromosome, '\t', prev_position, '\t', position, '\t',
                                                    sum(t_list_b) / sum(n_list_b), '\t', len(t_list_b), '\t',
                                                    sum(t_list_b), '\t', sum(n_list_b))
                t_list_b = []
                n_list_b = []
                current_chromosome = chromosome
                prev_position = 0
                continue
            elif position > (prev_position + options.block_size):

                if options.excludeoutliers:
                    location = str(chromosome) + '_' + str(position)
                    if location in bannedlocations:
                        t_list_b = []
                        n_list_b = []
                        prev_position = position
                        continue

                sumt = sumt_main - sum(t_list_b)
                sumn = sumn_main - sum(n_list_b)
                if options.bootstrap != 0:
                    if len(t_list_b) > 1:
                        sumt = sum(t_list_b)
                        sumn = sum(n_list_b)
                        t_list_bootstrap.append(sumt)
                        n_list_bootstrap.append(sumn)

                D_pop = sumt / sumn
                mjlist.append(len(t_list_b))

                if options.verboseblocks: print(chromosome, '\t', prev_position, '\t', position, '\t',
                                                sum(t_list_b) / sum(n_list_b), '\t', len(t_list_b), '\t',
                                                sum(t_list_b), '\t', sum(n_list_b), '\t',
                                                Dp_main - sum(t_list_b) / sum(n_list_b))
                jackknife_Dp.append(D_pop)

                t_list_b = []
                n_list_b = []
                prev_position = position

    return t_list_bootstrap, n_list_bootstrap, jackknife_Dp, mjlist
