# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 16:14
@Author : Shuya Zhang, Ngenie Liang
@File : utils.py
@Project : popstats
"""
import sys
import logging

logger = logging.getLogger('my_log')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
logger.addHandler(ch)


def progress(x):
    """
    Print progress bar
    """
    # Assemble a string, equivalent to x + ' SNPs', f'{x} SNPs'
    out = '%s SNPs' % x  # The output
    bs = '\b' * 1000  # The backspace
    print(bs, end=' ', file=sys.stderr)
    print(out, end=' ', file=sys.stderr)


def output_log(options, poplabel2, poplabel3, poplabel4, targetpop1count, targetpop2count, t_list, n_list,
               abbalist, babalist):
    if options.autocorr:
        print('dABBA:', abbalist.count('dABBA'))
        print('ndABBA:', abbalist.count('ndABBA'))
        print('dBABA:', babalist.count('dBABA'))
        print('ndBABA:', babalist.count('ndBABA'))
        abbacount = 1.0 * abbalist.count('dABBA') / (abbalist.count('ndABBA') + abbalist.count('dABBA'))
        babacount = 1.0 * babalist.count('dBABA') / (babalist.count('ndBABA') + babalist.count('dBABA'))
        D = (1.0 * abbacount - babacount) / (babacount + abbacount)
        print(options.LD, '\t', abbacount, '\t', abbacount, D, abbalist.count('dABBA'),
              abbalist.count('ndABBA'),
              abbalist.count('ndBABA'), abbalist.count('dBABA'))
    # exit(0)

    if options.SFS or options.SFS2:
        print('+'.join(poplabel2) + '_SFS', '\t', '+'.join(poplabel4) + '=0 ' + '+'.join(poplabel3) + '=1',
              '\t',
              '+'.join(poplabel3) + '=0 ' + '+'.join(poplabel4) + '=1')
        for i in range(0, targetpop2count):
            freq1 = t_list.count(i) * 1.0
            freq2 = n_list.count(i) * 1.0
            n = 1.0 * freq1 + freq2
            # print freq1,freq2
            if freq1 == freq2:
                freq = 0.0
            else:
                freq = (1.0 * freq1 - 1.0 * freq2) / n
            # val=(abs(freq)*(1.0-abs(freq)))/n
            # print val
            # SE=sqrt(val)

            print(i, '\t', t_list.count(i), '\t', n_list.count(i), '\t', freq)  # ,'\t',SE,'\t',freq/SE
        exit(0)

    elif options.mSFS:
        total = len(t_list)
        for i in range(0, targetpop1count):
            for x in range(0, targetpop2count):
                freq1 = t_list.count((i, x)) * 1.0 / total
                freq2 = n_list.count((i, x)) * 1.0 / total
                n = 1.0 * freq1 + freq2
                if n == 0:
                    print('0', end=' ')
                    continue
                # print freq1,freq2
                freq = (1.0 * freq1 - 1.0 * freq2) / n
                # val=(abs(freq)*(1.0-abs(freq)))/n
                print(freq, end=' ')
            # SE=sqrt(val)

            # print i,x,'\t',t_list.count((i,x)),'\t',n_list.count((i,x)),'\t',freq,'\t',SE,'\t',freq/SE
            print('')

        exit(0)

