
import numpy as np
from matplotlib import pyplot as plt

def read_data(filename):
    with open(filename) as f:
        reg  = [int(x) for x in f.readline().strip().split(',')]
        hard = [int(x) for x in f.readline().strip().split(',')]
    return np.array(reg), np.array(hard)

def bin_data(ary):
    logary = np.log2(ary).astype(np.int)
    vals, counts = np.unique(logary, return_counts=True)

    # make sure we didn't skip any values
    if len(vals) != vals[-1]+1:
        rtn = np.zeros(vals[-1]+1, dtype=np.float)
        rtn[vals] = counts
    else:
        rtn = counts.astype(np.float)

    errs = np.sqrt(rtn)

    # normalize
    rtn /= len(ary)
    errs /= len(ary)

    return rtn, errs

def plot_data(reg, hard):
    reg, regerr = bin_data(reg)
    hard, harderr = bin_data(hard)
    idxs = np.arange(len(reg))

    barwidth = 0.4

    f, ax = plt.subplots(figsize=(3,2.5))

    plt.bar(idxs-barwidth/2, reg, barwidth, yerr=regerr, label='All')
    plt.bar(idxs+barwidth/2, hard, barwidth, yerr=harderr, label='"Hard"')

    ax.set_xticks(idxs)

    plt.xlabel(r"$n - \mathrm{rank}(M)$")
    plt.ylabel("Relative frequency")

    plt.legend()
    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
    #ncol=2, mode="expand", borderaxespad=0.)

    plt.tight_layout()

    plt.savefig('../../paper/candkeys.pdf')

    plt.show()

def main():
    data = read_data('candkeys.txt')
    plot_data(*data)

if __name__ == '__main__':
    main()
