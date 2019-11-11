
import matplotlib
from matplotlib import pyplot as plt

def read_data(filename="timing.csv"):
    ns = []
    ts = []
    iters = []
    errlows = []
    errhighs = []

    with open(filename) as f:
        for line in f:
            q, t, errlow, errhigh, niter = line.strip().split(',')
            ns.append((int(q)+1)//2)
            ts.append(float(t))
            iters.append(float(niter))
            errlows.append(float(errlow))
            errhighs.append(float(errhigh))

    return ns, ts, iters, errlows, errhighs

def plot_data(ns, ts, iters, errlows, errhighs):
    fig, (timeax, iterax) = plt.subplots(2, 1,
                                         sharex=True,
                                         figsize=(3.5, 3.5),
                                         gridspec_kw={'height_ratios' : [2.5, 1]})

    timelines = timeax.plot(ns, ts,
                            marker='o',
                            markersize=4,
                            #color='k',
                            label='Solve time')

    timeax.fill_between(ns, errlows, errhighs,
                        color='C0',
                        alpha=0.3,
                        linewidth=0,
                        #label=r'$1^{\mathrm{st}}-3^{\mathrm{rd}}$ quartile'
    )

    # fit
    c = sum(t/(n**2) for t,n in zip(ts, ns)) / len(ns)
    timelines += timeax.plot(ns, [c*n**2 for n in ns],
                             color='k',
                             linestyle='--',
                             label=r'$\propto n^2$',
                             zorder=0)

    timeax.set_xscale('log')
    timeax.set_yscale('log')

    iterlines = iterax.plot(ns, iters,
                            #color='k',
                            color='C2',
                            #linestyle='--',
                            label='# cand. keys')
    iterax.set_ylim(0, 9)

    iterax.set_xlabel('Problem size (number of qubits)')
    timeax.set_ylabel('Solve time [s]')
    iterax.set_ylabel('Keys\nchecked')

    # add more descriptive ticks
    plt.xticks([50, 100, 200, 500])
    iterax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    # add sub-figure labels
    timeax.text(0.98, 0.02, '(a)',
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=timeax.transAxes)

    iterax.text(0.98, 0.06, '(b)',
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=iterax.transAxes)

    plt.tight_layout()

    # # get all the lines in the same legend
    # lines = timelines + iterlines
    # labels = [l.get_label() for l in lines]
    # timeax.legend(lines, labels)
    timeax.legend()

    plt.savefig("../../paper/solvetime.pdf")
    plt.show()

def main():
    data = read_data()
    plot_data(*data)

if __name__ == '__main__':
    main()
