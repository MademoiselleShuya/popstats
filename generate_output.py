# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 18:30
@Author : Shuya Zhang, Ngenie Liang
@File : generate_output.py
@Project : popstats
"""
import random
from math import sqrt
from utils import logger


def output_boostrap(options, poplabel1, poplabel2, poplabel3, poplabel4, t_list_bootstrap, n_list_bootstrap,
                    num_SNPs, targetpop1count, targetpop2count, targetpop3count, targetpop4count, Dp_main):
    if options.bootstrap != 0:
        bootDlist = []
        for x in range(0, 1000):
            pseudo_rep_t = []
            pseudo_rep_n = []
            for y in range(0, len(t_list_bootstrap)):
                pseudo_rep_t.append(random.choice(t_list_bootstrap))
                pseudo_rep_n.append(random.choice(n_list_bootstrap))
            bootD = 1.0 * sum(pseudo_rep_t) / sum(pseudo_rep_n)
            bootDlist.append(bootD)

        n_groups = float(len(t_list_bootstrap))

        from rpy2.robjects import r

        print('+'.join(poplabel1), '\t', '+'.join(poplabel2), '\t', '+'.join(poplabel3), '\t',
              '+'.join(poplabel4), '\t', float(Dp_main), '\t',
              r.quantile(bootDlist, probs=r.c(0.025, 0.975)),
              '\t', num_SNPs, '\t', int(n_groups), '\t', targetpop1count, '\t', targetpop2count, '\t',
              targetpop3count, '\t', targetpop4count)
        print('ratio %CI\t', end=' ')
        exit(0)
    else:
        logger.info('options.bootstrap is 0')


def output_block_jackknife(options, poplabel1, poplabel2, poplabel3, poplabel4, Dp_main, Dp_main2,
                           jackknife_Dp, mjlist, num_SNPs, targetpop1count, targetpop2count, targetpop3count,
                           targetpop4count, blist):
    if options.nickverbose:
        print('Total:', Dp_main2, Dp_main, '0')
        exit(0)

    if options.noweighting:
        n_groups = float(len(jackknife_Dp))
        pseudovalues = []
        for j in jackknife_Dp:
            pseudovalues.append(float((n_groups * Dp_main - (n_groups - 1.0) * j)))
        jackknife_estimate = float(sum(pseudovalues)) / n_groups

        sum_list = []
        for f in pseudovalues:
            sum_list.append(float(f - jackknife_estimate) ** 2)

        Dp_SE = sqrt(sum(sum_list) / (n_groups * (n_groups - 1.0)))

    else:
        n_groups = float(len(jackknife_Dp))

        pseudovalues = []
        for m, Dj in zip(mjlist, jackknife_Dp):
            pseudovalues.append(((num_SNPs - m) * Dj) / num_SNPs)
        jackknife_estimate = n_groups * Dp_main - float(sum(pseudovalues))

        # print jackknife_estimate
        jvallist = []
        for mj, Dj in zip(mjlist, jackknife_Dp):
            hj = num_SNPs / mj
            thetaminusj = Dj
            tj = hj * Dp_main - (hj - 1.0) * thetaminusj

            jval = ((tj - jackknife_estimate) ** 2) / (hj - 1.0)
            jvallist.append(jval)
        Dp_SE = sum(jvallist) / n_groups
        Dp_SE = sqrt(Dp_SE)

    if options.jackest:
        Dp_main = jackknife_estimate

    print('+'.join(poplabel1), '\t', '+'.join(poplabel2), '\t', '+'.join(poplabel3), '\t',
          '+'.join(poplabel4),
          '\t', float(Dp_main), '\t', Dp_SE, '\t', float(Dp_main) / float(Dp_SE), '\t', num_SNPs, '\t',
          int(n_groups), '\t', targetpop1count, '\t', targetpop2count, '\t', targetpop3count, '\t',
          targetpop4count)  # ,sum(t_list), sum(n_list)
    # exit(0)
    if options.Bcorr or options.Bstrat:
        # print '\t'.join([str(r) for r in range(0,9)])
        print('\t'.join([str(r) for r in blist]))
        print(jackknife_estimate)
