
import numpy as np

from matplotlib import pyplot as plt
plt.style.use('norm.mplstyle')

def read_cand_data(filename):
    with open(filename) as f:
        reg  = [int(x) for x in f.readline().strip().split(',')]
        hard = [int(x) for x in f.readline().strip().split(',')]
    return np.array(reg), np.array(hard)

def read_timing_data(filename):
    ns = []
    candkeys = []

    with open(filename) as f:
        for line in f:
            q, _, _, _, nkeys = line.strip().split(',')
            ns.append((int(q)+1)//2)
            candkeys.append(float(nkeys))

    return ns, candkeys

def bin_data(ary):
    logary = np.log2(ary).astype(int)
    vals, counts = np.unique(logary, return_counts=True)

    # make sure we didn't skip any values
    if len(vals) != vals[-1]+1:
        rtn = np.zeros(vals[-1]+1, dtype=float)
        rtn[vals] = counts
    else:
        rtn = counts.astype(float)

    errs = np.sqrt(rtn)

    # normalize
    rtn /= len(ary)
    errs /= len(ary)

    return rtn, errs

def plot_candkeys(ax, ns, candkeys):
    ax.plot(ns, candkeys,
            #color='k',
            color='C2',
            marker='o',
            markersize=2,
            #linestyle='--',
            label='# cand. keys')

    ax.set_ylim(0, 9)

    ax.set_xlabel('Problem size $n$')
    ax.set_ylabel('Candidates\nchecked')

def plot_hist(ax, reg, hard):
    reg, regerr = bin_data(reg)
    hard, harderr = bin_data(hard)
    idxs = np.arange(len(reg))

    barwidth = 0.4

    ax.bar(idxs-barwidth/2, reg, barwidth, yerr=regerr, label='All')
    ax.bar(idxs+barwidth/2, hard, barwidth, yerr=harderr, label='Re-run')

    ax.set_xticks(idxs)

    ax.set_xlabel(r"$n - \mathrm{rank}(M)$")
    ax.set_ylabel("Relative frequency")

    # ax.text(0.96, 0.7, '$n = 245$',
    #         horizontalalignment='right',
    #         verticalalignment='top',
    #         transform=ax.transAxes)

    ax.legend()

def plot_data(reg, hard, ns, candkeys):

    f, (nkeyax, rankax) = plt.subplots(2, 1,
                                       figsize=(4,5.3),
                                       gridspec_kw={'height_ratios': [1, 3]})

    plot_candkeys(nkeyax, ns, candkeys)
    plot_hist(rankax, reg, hard)

    # add sub-figure labels
    nkeyax.text(0.02, 0.94, '(a)',
                horizontalalignment='left',
                verticalalignment='top',
                transform=nkeyax.transAxes)

    rankax.text(0.02, 0.94, '(b)',
                horizontalalignment='left',
                verticalalignment='top',
                transform=rankax.transAxes)

    plt.tight_layout()

    plt.savefig('candkeys.pdf')

    plt.show()

def main():
    data = read_cand_data('candkeys.csv')
    keydata = read_timing_data('timing.csv')
    plot_data(*(data+keydata))

if __name__ == '__main__':
    main()
