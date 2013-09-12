from motifs import *
from utils import *
from matplotlib import pyplot as plt
from matplotlib import rc
from copy_numbers import copy_numbers
import scikits.statsmodels.api as sm
import numpy as np

rc('text',usetex=True)

rseqs = [motif_ic(getattr(Escherichia_coli,tf)) for tf in Escherichia_coli.tfs]
rfreqs = [log2(4.6*10**6/len(getattr(Escherichia_coli,tf)))
          for tf in Escherichia_coli.tfs]

def make_plot(filename=None):
    
    plt.scatter(rfreqs,rseqs)
    plt.plot([0,20],[0,20],linestyle='--')
    plt.xlabel(r"R_{freq}")
    plt.ylabel(r"R_{seq}")
    plt.title(r"\textrm{R}_{freq}\ vs.\ \textrm{R}_{seq}\ \textrm{for}\ \textrm{27}\ \textit{E.\ coli}\ \textrm{TFs}")
    maybesave(filename)

def make_plot_with_copy_numbers(filename=None):
    r_rseqs = [motif_ic(getattr(Escherichia_coli,tf)) for tf in Escherichia_coli.tfs
               if tf in copy_numbers]
    r_rfreqs = [log2(4.6*10**6/len(getattr(Escherichia_coli,tf)))
          for tf in Escherichia_coli.tfs
                if tf in copy_numbers]
    copies = [copy_numbers[tf] for tf in Escherichia_coli.tfs if tf in copy_numbers]
    adj_rseqs = zipWith(lambda x,y:x+log2(y),r_rseqs,copies)
    plt.scatter(r_rfreqs,adj_rseqs)
    plt.plot([0,20],[0,20],linestyle='--')
    plt.xlabel(r"R_{freq}")
    plt.ylabel(r"R_{seq} + log_2(N)")
    #plt.title(r"\textrm{R}_{freq}\ vs.\ \textrm{R}_{seq}\ \textrm{for}\ \textrm{27}\ \textit{E.\ coli}\ \textrm{TFs}")
    maybesave(filename)

def plot_r_seq_against_copy_numbers(filename=None):
    r_rseqs = [motif_ic(getattr(Escherichia_coli,tf)) for tf in Escherichia_coli.tfs
               if tf in copy_numbers]
    r_rfreqs = [log2(4.6*10**6/len(getattr(Escherichia_coli,tf)))
          for tf in Escherichia_coli.tfs
                if tf in copy_numbers]
    copies = [copy_numbers[tf] for tf in Escherichia_coli.tfs if tf in copy_numbers]
    adj_rseqs = zipWith(lambda x,y:x+log2(y),r_rseqs,copies)
    plt.scatter(r_rseqs,map(log2,copies))
    plt.plot([0,20],[0,20],linestyle='--')
    plt.xlabel(r"R_{seq}")
    plt.ylabel(r"log_2(N)")
    plt.title(r"\textrm{R}_{freq}\ vs.\ \textrm{R}_{seq}\ \textrm{for}\ \textrm{27}\ \textit{E.\ coli}\ \textrm{TFs}")
    maybesave(filename)

def plot_r_freq_against_copy_numbers(filename=None):
    r_rseqs = [motif_ic(getattr(Escherichia_coli,tf)) for tf in Escherichia_coli.tfs
               if tf in copy_numbers]
    r_rfreqs = [log2(4.6*10**6/len(getattr(Escherichia_coli,tf)))
          for tf in Escherichia_coli.tfs
                if tf in copy_numbers]
    copies = [copy_numbers[tf] for tf in Escherichia_coli.tfs if tf in copy_numbers]
    adj_rseqs = zipWith(lambda x,y:x+log2(y),r_rseqs,copies)
    plt.scatter(r_rfreqs,map(log2,copies))
    plt.plot([0,20],[0,20],linestyle='--')
    plt.xlabel(r"R_{freq}")
    plt.ylabel(r"log_2(N)")
    plt.title(r"\textrm{R}_{freq}\ vs.\ \textrm{R}_{seq}\ \textrm{for}\ \textrm{27}\ \textit{E.\ coli}\ \textrm{TFs}")
    maybesave(filename)

def explain_rseq_by_rfreq_and_copy():
    r_rseqs = [motif_ic(getattr(Escherichia_coli,tf)) for tf in Escherichia_coli.tfs
               if tf in copy_numbers]
    r_rfreqs = [log2(4.6*10**6/len(getattr(Escherichia_coli,tf)))
          for tf in Escherichia_coli.tfs
                if tf in copy_numbers]
    copies = [copy_numbers[tf] for tf in Escherichia_coli.tfs if tf in copy_numbers]
    log_copies = map(log2,copies)
    X = sm.add_constant(np.column_stack((r_rfreqs,log_copies)),prepend=True)
    res = sm.OLS(r_rseqs,X).fit()
    print res.summary()
