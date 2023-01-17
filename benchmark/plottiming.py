
import matplotlib
from matplotlib import pyplot as plt

def read_data(filename="timing.csv"):
    ns = []
    ts = []
    errlows = []
    errhighs = []

    with open(filename) as f:
        for line in f:
            q, t, errlow, errhigh, _ = line.strip().split(',')
            ns.append((int(q)+1)//2)
            ts.append(float(t))
            errlows.append(float(errlow))
            errhighs.append(float(errhigh))

    return ns, ts, errlows, errhighs

def plot_data(ns, ts, errlows, errhighs):
    f, ax = plt.subplots(figsize=(3.5, 3))

    ax.plot(ns, ts,
            marker='o',
            linestyle='',
            markersize=5,
            color='dodgerblue',
            mec="k",
            markeredgewidth=1,
            label='Solve time')

    ax.fill_between(ns, errlows, errhighs,
                    color='dodgerblue',
                    alpha=0.3,
                    linewidth=0,
                    #label=r'$1^{\mathrm{st}}-3^{\mathrm{rd}}$ quartile'
    )

    # fit
    c = sum(t/(n**2) for t,n in zip(ts, ns)) / len(ns)
    ax.plot(ns, [c*n**2 for n in ns],
            color='k',
            linestyle='--',
            label=r'$\propto n^2$',
            zorder=1)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Problem size $n$ (qubit count)')
    ax.set_ylabel('Solve time [s]')

    # add more descriptive ticks
    plt.xticks([50, 100, 200, 500])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    plt.tight_layout()

    ax.legend()

    plt.savefig("solvetime.pdf")
    plt.show()

def main():
    data = read_data()
    plot_data(*data)

if __name__ == '__main__':
    main()
