# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 15:56
@Author : Shuya Zhang, Ngenie Liang
@File : main.py
@Project : popstats
"""
from parameters import parser
from generate_output import output_boostrap, output_block_jackknife
from block_jackknife import make_dp_main, block_jackknife
from process_tped_data import process_tped_data
from parse_parameters import parse_input_files_path, parse_pops, parse_mutation_parameters, \
    read_tfam, make_targetpop, make_targetpop2, parse_ancestor, parse_testpop, parse_ascertain123
from utils import logger, output_log


# options, args = parser.parse_args()


def popstats_main(options):
    """
    main function
    """
    # can't find option pi
    if options.FST:
        logger.info(f'***** Start FST computations *****')
    elif options.FSTWC:
        logger.info(f'***** Start FSTWC computations *****')
    elif options.f2:
        logger.info(f'***** Start f2 computations *****')
    elif options.f3:
        logger.info(f'***** Start f3 computations *****')
    elif options.f3vanilla:
        logger.info(f'***** Start f3vanilla computations *****')
    elif options.f4:
        logger.info(f'***** Start f4 computations *****')
    elif options.ratio:
        logger.info(f'***** Start f4-ratio computations *****')
    elif options.f5:
        logger.info(f'***** Start f5 computations *****')
    elif options.symmetry:
        logger.info(f'***** Start symmetry computations *****')
    elif options.LD:
        logger.info(f'***** Start LD computations *****')
    elif options.ancestor:
        logger.info(f'***** Start FAB computations *****')
    else:
        logger.info(f'***** Start D computations *****')

    logger.info('***** Parsing parameters *****')
    samples, data = parse_input_files_path(options)
    poplabel1, poplabel2, poplabel3, poplabel4 = parse_pops(options)
    mutationclass, regionchoice = parse_mutation_parameters(options)
    pops, popinds = read_tfam(samples, options)
    targetpop1, targetpop2, targetpop3, targetpop4 = make_targetpop(pops, popinds, options, poplabel1,
                                                                    poplabel2,
                                                                    poplabel3, poplabel4)

    popcount = len(pops)
    targetpop1count = len(targetpop1)
    targetpop2count = len(targetpop2)
    targetpop3count = len(targetpop3)
    targetpop4count = len(targetpop4)

    targetpop21, targetpop22, targetpop23, targetpop24 = make_targetpop2(options, pops)
    ancpop = parse_ancestor(options, pops)
    testpop = parse_testpop(options, pops, popinds)
    ascpop, ascpopcount, ascpop2, ascpopcount2, ascpop3, ascpopcount3 = parse_ascertain123(options, pops)

    logger.info('***** Processing tped data *****')
    choice_list, t_list, n_list, abbalist, babalist = process_tped_data(data, options, regionchoice,
                                                                        popcount, targetpop1,
                                                                        targetpop2, targetpop3, targetpop4,
                                                                        targetpop1count, targetpop2count,
                                                                        targetpop3count,
                                                                        targetpop4count, targetpop21,
                                                                        targetpop22, targetpop23,
                                                                        targetpop24, ancpop, testpop,
                                                                        ascpop, ascpopcount, ascpop2,
                                                                        ascpopcount2, ascpop3, ascpopcount3,
                                                                        mutationclass)
    output_log(options, poplabel2, poplabel3, poplabel4, targetpop1count, targetpop2count, t_list, n_list,
               abbalist, babalist)

    # logger.info('***** Computing block jackknife *****')
    Dp_main, blist, Dp_main2, num_SNPs = make_dp_main(options, poplabel1, poplabel2, poplabel3, poplabel4,
                                                      t_list, n_list, choice_list)
    t_list_bootstrap, n_list_bootstrap, jackknife_Dp, mjlist = block_jackknife(options, t_list, n_list,
                                                                               choice_list)

    # logger.info('***** Output boostrap *****')
    output_boostrap(options, poplabel1, poplabel2, poplabel3, poplabel4, t_list_bootstrap, n_list_bootstrap,
                    num_SNPs, targetpop1count, targetpop2count, targetpop3count, targetpop4count, Dp_main)

    # logger.info('***** Output block jackknife *****')
    logger.info('***** Output *****')
    output_block_jackknife(options, poplabel1, poplabel2, poplabel3, poplabel4, Dp_main, Dp_main2,
                           jackknife_Dp, mjlist, num_SNPs, targetpop1count, targetpop2count, targetpop3count,
                           targetpop4count, blist)

    logger.info(f'***** Analysis done! *****')


if __name__ == '__main__':
    opts, args = parser.parse_args()
    popstats_main(opts)
