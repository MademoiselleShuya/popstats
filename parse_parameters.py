# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 15:59
@Author : Shuya Zhang, Ngenie Liang
@File : parse_parameters.py
@Project : popstats
"""
import sys
from utils import logger


def parse_input_files_path(options):
    """
    Parse parameters of input files: tfam file' path and tped file' path
    """
    options.chromosome = options.chromosome.lstrip('chr')
    if options.chromosome == 'autosomes' or options.chromosome == '0':
        options.chromosome = False
    else:
        options.chromosome = int(options.chromosome)

    if options.outfile != False:
        sys.stdout = open(options.outfile, 'w')

    if options.file:
        samples = options.file + '.tfam'
        data = options.file + '.tped'
    else:
        samples = options.tfam
        data = options.tped

    return samples, data


def parse_pops(options):
    """
    Parse parameters of populations: poplist, poplabel1, poplabel2, poplabel3, poplabel4
    """
    poplist = options.pops.split(',')

    # print a poplist
    logger.info(f'***** Poplist: {poplist} *****')

    poplabel1 = poplist[0].split('+')
    poplabel2 = poplist[1].split('+')
    poplabel3 = poplist[2].split('+')

    if options.f3:
        poplabel4 = poplabel3
    elif options.f3vanilla:
        poplabel4 = poplabel3
    else:
        poplabel4 = poplist[3].split('+')

    return poplabel1, poplabel2, poplabel3, poplabel4


def define_genetic_distance(options):
    """
    Parse parameters of genetic distance: block_size
    """
    # Use genetic distance (default 5 cM) instead of physical distance to define block size for the
    # jackknife.
    if options.morgan:
        options.block_size = 0.05
    block_size = options.block_size
    return block_size


def parse_mutation_parameters(options):
    """
    Parse parameters of mutation: mutationrate (alleles in comma-delimited list), regionchoice
    """
    # 需要补充解释
    mutationclass = []
    regionchoice = []
    if options.mutationclass:
        mutationclass = options.mutationclass.split(',')

    if options.region:
        regionchoice = [int(x) for x in options.region.split('-')]
    return mutationclass, regionchoice


def read_tfam(samples, options):
    """
    Read the tfam file
    """
    pops = []
    # popinds is a list that combines the sample name and ID in the form of col[0] + ':' + col[1]
    popinds = []
    for line in open(samples):
        col = line.split()
        # Distinguish sex
        sex = col[4]
        if options.nomales and sex == '1':
            pops.append('Ignore')
            pops.append('Ignore')
            continue
        if options.onlymales and sex == '2':
            pops.append('Ignore')
            pops.append('Ignore')
            continue
        pops.append(col[0])
        pops.append(col[0])

        # A combination of sample names and IDs
        popinds.append(col[0] + ':' + col[1])
        popinds.append(col[0] + ':' + col[1])
    return pops, popinds


def make_targetpop(pops, popinds, options, poplabel1, poplabel2, poplabel3, poplabel4):
    """
    Make a list of target populations
    """
    targetpop1 = []
    targetpop2 = []
    targetpop3 = []
    targetpop4 = []

    for i, p in enumerate(pops):
        if p in poplabel1 and (len(targetpop1) < options.maxn * 2): targetpop1.append(i)
        if p in poplabel2 and (len(targetpop2) < options.maxn * 2): targetpop2.append(i)
        if p in poplabel3 and (len(targetpop3) < options.maxn * 2): targetpop3.append(i)
        if p in poplabel4 and (len(targetpop4) < options.maxn * 2): targetpop4.append(i)

    for i, p in enumerate(popinds):
        if p in poplabel1 and (len(targetpop1) < options.maxn * 2): targetpop1.append(i)
        if p in poplabel2 and (len(targetpop2) < options.maxn * 2): targetpop2.append(i)
        if p in poplabel3 and (len(targetpop3) < options.maxn * 2): targetpop3.append(i)
        if p in poplabel4 and (len(targetpop4) < options.maxn * 2): targetpop4.append(i)

    if len(targetpop1) < 1:
        print('no', poplabel1)
        exit(0)
    if len(targetpop2) < 1:
        print('no', poplabel2)
        exit(0)
    if len(targetpop3) < 1:
        print('no', poplabel3)
        exit(0)
    if len(targetpop4) < 1:
        print('no', poplabel4)
        exit(0)
    return targetpop1, targetpop2, targetpop3, targetpop4


def make_targetpop2(options, pops):
    targetpop21 = []
    targetpop22 = []
    targetpop23 = []
    targetpop24 = []
    if options.pops2:
        poplist2 = options.pops2.split(',')
        poplabel21 = poplist2[0].split('+')
        poplabel22 = poplist2[1].split('+')
        poplabel23 = poplist2[2].split('+')
        poplabel24 = poplist2[3].split('+')

        for i, p in enumerate(pops):
            if p in poplabel21 and (len(targetpop21) < options.maxn * 2): targetpop21.append(i)
            if p in poplabel22 and (len(targetpop22) < options.maxn * 2): targetpop22.append(i)
            if p in poplabel23 and (len(targetpop23) < options.maxn * 2): targetpop23.append(i)
            if p in poplabel24 and (len(targetpop24) < options.maxn * 2): targetpop24.append(i)

    return targetpop21, targetpop22, targetpop23, targetpop24


def parse_ancestor(options, pops):
    """
    Parse parameters of ancestor: ancpop
    """
    ancpop = []
    if options.ancestor:
        anclabel = options.ancestor.split('+')
        for i, p in enumerate(pops):
            if p in anclabel: ancpop.append(i)
    return ancpop


def parse_testpop(options, pops, popinds):
    """
    Parse parameters of testpop: testpop
    """
    testpop = []
    if options.testpop:
        testlabel = options.testpop.split('+')
        for i, p in enumerate(pops):
            if p in testlabel: testpop.append(i)
        for i, p in enumerate(popinds):
            if p in testlabel: testpop.append(i)
    return testpop


def parse_ascertain123(options, pops):
    """
    Parse parameters of ascertain: ascertain1, ascertain2, ascertain3
    """

    def parse_ascertain(pops, ascertain, downsampleasc):
        """
        Parse parameters of ascertain: ascertain
        """
        ascpop = []
        ascpopcount = 0
        if ascertain:
            asclabel = ascertain.split('+')
            for i, p in enumerate(pops):
                if p in asclabel: ascpop.append(i)

            ascpopcount = len(ascpop)
            if downsampleasc:
                ascpopcount = downsampleasc
        return ascpop, ascpopcount

    ascpop, ascpopcount = parse_ascertain(pops, options.ascertain, options.downsampleasc)
    ascpop2, ascpopcount2 = parse_ascertain(pops, options.ascertain2, options.downsampleasc2)
    ascpop3, ascpopcount3 = parse_ascertain(pops, options.ascertain3, options.downsampleasc3)
    return ascpop, ascpopcount, ascpop2, ascpopcount2, ascpop3, ascpopcount3
