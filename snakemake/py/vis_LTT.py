import logging
import os
from collections import Counter, defaultdict

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import plot, xlim, xlabel, ylabel, legend, show, figure, savefig
from pastml.tree import read_tree, DATE, annotate_dates

DATE_STEP = 10

EXT_COLOR = '#4daf4a'
HIGH_COLOR = '#e41a1c'
LOW_COLOR = '#377eb8'

state2color = {'High': HIGH_COLOR, 'External': EXT_COLOR, 'Low': LOW_COLOR}


def plot_proportions_years(tree, col, ax=None, suffix='', linestyle='solid'):
    high, low, ext = Counter(), Counter(), Counter()

    year2n = defaultdict(set)
    for n in tree.traverse():
        year2n[int(getattr(n, DATE))].add(n)
    years = sorted(year2n.keys())

    for year in years:
        ns = year2n[year]
        for n in ns:
            states = getattr(n, col)
            if 'High' in states:
                high[year] += 1 / len(states)
            if 'Low' in states:
                low[year] += 1 / len(states)
            if 'External' in states:
                ext[year] += 1 / len(states)

    xs = np.array(sorted(set(high.keys()) | set(low.keys()) | set(ext.keys())))
    high = np.array([high[x] for x in xs])
    low = np.array([low[x] for x in xs])
    ext = np.array([ext[x] for x in xs])

    if ax:
        ax.plot(xs, low, color=LOW_COLOR, label='Low{}'.format(suffix), linestyle=linestyle)
        ax.plot(xs, high, color=HIGH_COLOR, label='High{}'.format(suffix), linestyle=linestyle)
        ax.plot(xs, ext, color=EXT_COLOR, label='External{}'.format(suffix), linestyle=linestyle)
        ax.set_xlim(min(xs) - .6, max(xs) + .6)
        ax.set_xlabel('Year')
        ax.set_ylabel('Number of infected individuals')
        ax.legend()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    else:
        plot(xs, low, color=LOW_COLOR, label='Low')
        plot(xs, high, color=HIGH_COLOR, label='High')
        plot(xs, ext, color=EXT_COLOR, label='External')
        xlim(min(xs) - .6, max(xs) + .6)
        xlabel('Year')
        ylabel('Number of infected individuals')
        legend()


def plot_sampling(tree, col, ax=None, suffix='', linestyle='solid'):
    state2tips = defaultdict(list)
    state2year2tips = defaultdict(lambda: defaultdict(list))

    min_year, max_year = np.inf, -np.inf
    for tip in tree:
        states = getattr(tip, col)
        date = getattr(tip, DATE)
        if len(states) == 1:
            state = next(iter(states))
            state2tips[state].append(tip.name)
            year = int(date)
            if year < min_year:
                min_year = year
            if year > max_year:
                max_year = year
            state2year2tips[state][year].append(tip.name)

    for state in ['Low', 'High', 'External']:
        year2tips = state2year2tips[state]
        xs = sorted(year2tips.keys())
        ys = [len(year2tips[_]) for _ in xs]
        for i in range(1, len(ys)):
            ys[i] += ys[i - 1]
        ax.plot(xs, ys, color=state2color[state], label=state + suffix, linestyle=linestyle)

    ax.set_xlabel('Year')
    ax.set_ylabel('Accumulated sampled cases')
    ax.legend()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(min_year - .6, max_year + .6)
    # ax.set_xticks(np.arange(min_year, max_year, step=5))
    # ax.set_yticks(np.arange(0, max(len(_) for _ in state2tips.values()) + 50, step=50))


if '__main__' == __name__:

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)
    import argparse

    data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../..', 'data')

    parser = argparse.ArgumentParser(description="Visualises PastML stats.")

    parser.add_argument('--column', default='highlow_prevalence', required=True,
                        type=str, help="the column of interest.")
    parser.add_argument('--trees', default=os.path.join(data_dir, "rep_*"),
                        type=str, help="the PASTML trees.", nargs='+')
    parser.add_argument('--labels', type=str, help="the PASTML tree labels.", nargs='+')
    parser.add_argument('--time_pdf', default=os.path.join(data_dir, "infections_time.pdf"), type=str, required=True,
                        help="the number of infected individuals in each state vs time.")
    parser.add_argument('--png', default=os.path.join(data_dir, "infections.png"), type=str, required=True,
                        help="LTT plot for the full and subsampled tree 1")

    params = parser.parse_args()

    forest = []
    for nwk in params.trees:
        tree = read_tree(nwk, columns=[params.column])
        annotate_dates([tree])
        forest.append(tree)

    # Infection plots
    with PdfPages(params.time_pdf) as pdf_pages:
            for label, tree in zip(params.labels, forest):
                fig = figure(figsize=(10, 10), dpi=100)
                plot_proportions_years(tree, params.column)
                pdf_pages.savefig(fig)

                # Done with the page
                logging.info('Analysed infections {}'.format(label))

    fig = figure(figsize=(12, 6), dpi=300)
    ax2, ax1 = fig.subplots(1, 2)
    plot_proportions_years(forest[0], params.column, ax=ax1)
    plot_proportions_years(forest[1], params.column, ax=ax1, suffix=' (subsampled)', linestyle='dashed')
    plot_sampling(forest[0], params.column, ax=ax2)
    plot_sampling(forest[1], params.column, ax=ax2, suffix=' (subsampled)', linestyle='dashed')

    x_m, x_M = min(ax1.get_xlim()[0], ax2.get_xlim()[0]), max(ax1.get_xlim()[1], ax2.get_xlim()[1])
    for ax in (ax1, ax2):
        ax.set_xlim(x_m, x_M)
    ax1.set_title('Infected individuals')
    ax2.set_title('Sampled individuals')

    savefig(params.png, bbox_inches='tight')
