# -*- coding: UTF-8 -*-
"""
@Time : 2024/8/3 15:57
@Author : Shuya Zhang, Ngenie Liang
@File : parameters.py
@Project : popstats
"""
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--tfam", action="store", type="string", dest="tfam", help="tfam_file")
parser.add_option("-p", "--tped", action="store", type="string", dest="tped", help="tped_file")
parser.add_option("--file", action="store", type="string", dest="file", help="file", default=False)
parser.add_option("-c", "--chromosome", action="store", type="string", dest="chromosome", help="chromosome",
                  default='0')
parser.add_option("--pops", action="store", dest="pops", help="populations in comma-delimited list")
parser.add_option("--pops2", action="store", dest="pops2", help="populations in comma-delimited list")

parser.add_option("--informative", action="store_true", dest="informative",
                  help="only use sites where pop1,pop2 are polymorphic and pop3,pop4 are also polymorphic",
                  default=False)

parser.add_option("--ancestor", action="store", dest="ancestor", help="ancestor", default=False)
parser.add_option("--Dcorr", action="store_true", dest="Dcorr", help="Dcorr", default=False)
parser.add_option("--scramble", action="store_true", dest="scramble", help="scramble", default=False)

parser.add_option("--jackest", action="store_true", dest="jackest",
                  help="use jackknife estimate of the mean", default=False)

parser.add_option("--scrambleall", action="store_true", dest="scrambleall", help="scrambleall",
                  default=False)
parser.add_option("--onlyABBA", action="store_true", dest="onlyABBA", help="onlyABBA", default=False)
parser.add_option("--onlyBABA", action="store_true", dest="onlyBABA", help="onlyBABA", default=False)
parser.add_option("--ABBAorBABA", action="store_true", dest="ABBAorBABA", help="ABBAorBABA", default=False)

parser.add_option("--allhaploid", action="store_true", dest="allhaploid", help="allhaploid", default=False)

parser.add_option("--Bcorr", action="store_true", dest="Bcorr",
                  help="Bcorr, requires B-stat info in marker field separated by comma", default=False)

parser.add_option("--Bstrat", action="store_true", dest="Bstrat",
                  help="Bstrat, requires B-stat info in marker field separated by comma", default=False)

parser.add_option("--nickverbose", action="store_true", dest="nickverbose", help="nickverbose",
                  default=False)

parser.add_option("--Bchoice", action="store", type="string", dest="Bchoice",
                  help="Bchoice, requires B-stat info in marker field separated by comma", default=False)

parser.add_option("--popscramble", action="store_true", dest="popscramble", help="popscramble",
                  default=False)
parser.add_option("--not23", action="store_true", dest="not23", help="not23", default=False)
parser.add_option("--chromblocks", action="store_true", dest="chromblocks", help="chromblocks",
                  default=False)
parser.add_option("--numSNPblocks", action="store_true", dest="numSNPblocks", help="numSNPblocks",
                  default=False)
parser.add_option("--normal", action="store_true", dest="normal", help="normal (placeholder option)",
                  default=False)
parser.add_option("--haploidize", action="store_true", dest="haploidize", help="haploidize", default=False)
parser.add_option("--maxn", action="store", type="int", dest="maxn", help="maxn", default=1000000000)
parser.add_option("--vcf", action="store", dest="vcf", help="vcf from stdin", default=False)
parser.add_option("--testpop", action="store", dest="testpop", help="testpop for f4 ratio", default=False)

parser.add_option("--region", action="store", dest="region", help="region start-stop", default=False)

parser.add_option("--bootstrap", action="store", type="int", dest="bootstrap", help="bootstrap", default=0)

parser.add_option("--ascertain", action="store", dest="ascertain", help="ascertain", default=False)
parser.add_option("--ascertainfreq", action="store", type="int", dest="ascertainfreq", help="ascertainfreq",
                  default=-1)
parser.add_option("--downsampleasc", action="store", type="int", dest="downsampleasc", help="downsampleasc",
                  default=False)

parser.add_option("--ascertain2", action="store", dest="ascertain2", help="ascertain2", default=False)
parser.add_option("--ascertainfreq2", action="store", type="int", dest="ascertainfreq2",
                  help="ascertainfreq2", default=-1)
parser.add_option("--downsampleasc2", action="store", type="int", dest="downsampleasc2",
                  help="downsampleasc2", default=False)

parser.add_option("--ascertain3", action="store", dest="ascertain3", help="ascertain3", default=False)
parser.add_option("--ascertainfreq3", action="store", type="int", dest="ascertainfreq3",
                  help="ascertainfreq3", default=-1)
parser.add_option("--downsampleasc3", action="store", type="int", dest="downsampleasc3",
                  help="downsampleasc3", default=False)

parser.add_option("--outdiff", action="store", type="float", dest="outdiff",
                  help="minimum freqdiff between A and B in (A,B,X,Y)", default=False)
parser.add_option("--indiff", action="store", type="float", dest="indiff",
                  help="minimum freqdiff between X and Y in (A,B,X,Y)", default=False)
parser.add_option("--fixed1", action="store_true", dest="fixed1", help="fixed1", default=False)
parser.add_option("--fixed3", action="store_true", dest="fixed3", help="fixed3", default=False)
parser.add_option("--fixed4", action="store_true", dest="fixed4", help="fixed4", default=False)

parser.add_option("--outdiffexact", action="store", type="float", dest="outdiffexact",
                  help="exact (0.1 etc) freqdiff between A and B in (A,B,X,Y)", default=False)

parser.add_option("--simpleD", action="store_true", dest="simpleD", help="simpleD", default=False)
parser.add_option("--simpleDfreq", action="store_true", dest="simpleDfreq", help="simpleDfreq",
                  default=False)
parser.add_option("--simpleDtailtest", action="store_true", dest="simpleDtailtest", help="simpleDtailtest",
                  default=False)
parser.add_option("--pop2weight", action="store_true", dest="pop2weight", help="pop2weight", default=False)

parser.add_option("--clock", action="store_true", dest="clock", help="clock", default=False)

parser.add_option("--excludeoutliers", action="store_true", dest="excludeoutliers", help="excludeoutliers",
                  default=False)
parser.add_option("--Dutheilfilter", action="store_true", dest="Dutheilfilter", help="Dutheilfilter",
                  default=False)
parser.add_option("--Sriramfilter", action="store_true", dest="Sriramfilter", help="Sriramfilter",
                  default=False)

parser.add_option("--nohzcorrection", action="store_true", dest="nohzcorrection",
                  help="Skip the Hz correction for f3, but keep sample size correction. To skip sample size correction, use --f3vanilla",
                  default=False)

parser.add_option("--haploidinclade", action="store_true", dest="haploidinclade",
                  help="randomly sampled allele from each of X and Y in (A,B,X,Y)", default=False)

parser.add_option("--mutationclass", action="store", dest="mutationclass",
                  help="alleles in comma-delimited list")

parser.add_option("--onlymales", action="store", dest="onlymales", help="onlymales", default=False)
parser.add_option("--nomales", action="store", dest="nomales", help="nomales", default=False)
parser.add_option("--pop2freq", action="store", type="int", dest="pop2freq", help="pop2freq", default=False)
parser.add_option("--notransitions", action="store_true", dest="notransitions", help="notransitions",
                  default=False)
parser.add_option("--nomissing", action="store_true", dest="nomissing", help="nomissing", default=False)
parser.add_option("--f4", action="store_true", dest="f4", help="f4", default=False)
parser.add_option("--f5", action="store_true", dest="f5", help="f5", default=False)
parser.add_option("--ratio", action="store_true", dest="ratio",
                  help="ratio [pop1,pop2,testpop,pop4] / [pop1,pop2,pop3,pop4]", default=False)
parser.add_option("--f3", action="store_true", dest="f3", help="f3", default=False)
parser.add_option("--f3vanilla", action="store_true", dest="f3vanilla", help="f3vanilla", default=False)
parser.add_option("--f2", action="store_true", dest="f2", help="f2", default=False)
parser.add_option("--verboseblocks", action="store_true", dest="verboseblocks", help="verboseblocks",
                  default=False)
parser.add_option("--verbose", action="store_true", dest="verbose", help="verbose", default=False)

parser.add_option("--mutatepop2", action="store", type="float", dest="mutatepop2", help="mutatepop2",
                  default=False)
parser.add_option("--mutatepop4", action="store", type="float", dest="mutatepop4", help="mutatepop4",
                  default=False)

parser.add_option("--countDzero", action="store_true", dest="countDzero", help="countDzero", default=False)

parser.add_option("--symmetry", action="store_true", dest="symmetry", help="symmetry", default=False)

parser.add_option("--LD", action="store", type="float", dest="LD", help="LD", default=False)
parser.add_option("--SNPfreq", action="store", type="string", dest="SNPfreq", help="SNPfreq", default=False)
parser.add_option("--FST", action="store_true", dest="FST", help="FST", default=False)
parser.add_option("--FSTWC", action="store_true", dest="FSTWC", help="FSTWC", default=False)
parser.add_option("--Tdiv", action="store_true", dest="Tdiv", help="Tdiv", default=False)
parser.add_option("--autocorr", action="store_true", dest="autocorr", help="autocorr", default=False)
parser.add_option("--autocorr2", action="store_true", dest="autocorr2", help="autocorr2", default=False)
parser.add_option("--positivestat", action="store_true", dest="positivestat", help="positivestat",
                  default=False)
parser.add_option("--polymorphic", action="store_true", dest="polymorphic", help="polymorphic",
                  default=False)

parser.add_option("--Dprim", action="store_true", dest="Dprim", help="Dprim", default=False)
parser.add_option("--Dr", action="store_true", dest="Dr", help="Dr", default=False)
parser.add_option("--LDwindow", action="store", type="float", dest="LDwindow", help="LDwindow", default=5000)
parser.add_option("--hapD", action="store_true", dest="hapD", help="hapD", default=False)
parser.add_option("--twohaps", action="store_true", dest="twohaps", help="twohaps", default=False)

parser.add_option("--mincount", action="store", type="float", dest="mincount", help="mincount", default=1)

parser.add_option("--equaln", action="store_true", dest="equaln", help="equaln", default=False)

parser.add_option("--LiReich", action="store_true", dest="LiReich",
                  help="Li and Reich statistic for probability that pop2 allele is derived given heterozygote in pop1",
                  default=False)

parser.add_option("--withinfreq", action="store_true", dest="withinfreq",
                  help="compute allele frequences for the LD test (A,B),(X,Y) for each population separately",
                  default=False)
parser.add_option("--withinoutgroupsfreq", action="store_true", dest="withinoutgroupsfreq",
                  help="compute allele frequences for the LD test (A,B),(X,Y) for population A and B separately but X and Y jointly",
                  default=False)
parser.add_option("--doubletest", action="store_true", dest="doubletest",
                  help="compute the sum of the LD4 stat and the sum of the two f4 stats for each pair",
                  default=False)

parser.add_option("--morgan", action="store_true", dest="morgan",
                  help="Use genetic distance (default 5 cM) instead of physical distance to define block size for the jackknife.",
                  default=False)
parser.add_option("--outfile", action="store", type="string", dest="outfile", help="outfile", default=False)
parser.add_option("--inds", action="store_true", dest="inds", help="inds", default=False)
parser.add_option("--multi", action="store_true", dest="multi", help="multi", default=False)
parser.add_option("--nojackknife", action="store_true", dest="nojackknife", help="nojackknife",
                  default=False)
parser.add_option("--noweighting", action="store_true", dest="noweighting", help="noweighting",
                  default=False)
parser.add_option("--SFS", action="store_true", dest="SFS", help="SFS test, [ancestor,SFStarget,pop3,pop4]",
                  default=False)
parser.add_option("--mSFS", action="store_true", dest="mSFS",
                  help="mSFS test, [SFStarget1,SFStarget2,pop3,pop4]", default=False)
parser.add_option("--SFS2", action="store_true", dest="SFS2",
                  help="SFS test, [ancestor,SFStarget,pop3,pop4]", default=False)
parser.add_option("-b", "--block_size", action="store", type="float", dest="block_size", help="block_size",
                  default=5000000.0)
parser.add_option("--anc_test", action="store_true", dest="anc_test",
                  help="anc test, [ancestor,SFStarget_fixed_ancestral,pop3,pop4]", default=False)

