import logging
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import bar, plot, xlim, xlabel, ylabel, legend, show, ylim, figure, savefig
from pastml.tree import read_tree, DATE, annotate_dates
from pastml.visualisation.cytoscape_manager import save_as_transition_html

DATE_STEP = 10

EXT_COLOR = '#4daf4a'
HIGH_COLOR = '#e41a1c'
LOW_COLOR = '#377eb8'

state2color = {'High': HIGH_COLOR, 'External': EXT_COLOR, 'Low': LOW_COLOR}


def count_transmissions(tree, mp_df, filter=lambda _: True):
    states = mp_df.columns
    from_to2count = Counter()
    for n in tree.traverse():
        if not filter(n):
            continue
        for n_state in states:
            same_state_count = 0
            for c in n.children:
                for c_state in states:
                    prob = mp_df.loc[n.name, n_state] * mp_df.loc[c.name, c_state]
                    from_to2count[(n_state, c_state)] += prob
                    if n_state == c_state:
                        same_state_count += prob
            old_v = from_to2count[(n_state, n_state)]
            from_to2count[(n_state, n_state)] -= min(mp_df.loc[n.name, n_state], same_state_count)
            if from_to2count[(n_state, n_state)] < 0:
                print(old_v, from_to2count[(n_state, n_state)], mp_df.loc[n.name, n_state], same_state_count)
    return from_to2count


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
    parser.add_argument('--mps', default=os.path.join(data_dir, "rep_*"),
                        type=str, help="the PASTML marginal probability files.", nargs='+')
    parser.add_argument('--labels', type=str, help="the PASTML tree labels.", nargs='+')
    parser.add_argument('--table', default=os.path.join(data_dir, "table.xlsx"), type=str, required=True,
                        help="Who infected whom table.")
    parser.add_argument('--out_html', default=os.path.join(data_dir, "transitions_{}.html"), type=str, required=True,
                        help="Who infected whom visualisation.")

    params = parser.parse_args()

    forest = []
    for nwk in params.trees:
        tree = read_tree(nwk, columns=[params.column])
        annotate_dates([tree])
        forest.append(tree)

    # Who infected whom
    with pd.ExcelWriter(params.table, engine='xlsxwriter') as writer:
        workbook = writer.book
        for label, tree, mp in zip(params.labels, forest, params.mps):
            mp_df = pd.read_csv(mp, sep='\t', index_col=0)
            states = mp_df.columns

            min_year, max_year = np.inf, -np.inf

            year2state_counts = defaultdict(Counter)
            state_counts = Counter()
            for tip in tree:
                year = int(getattr(tip, DATE))
                if min_year > year:
                    min_year = year
                if max_year < year:
                    max_year = year
                year //= DATE_STEP
                for state in states:
                    mp = mp_df.loc[tip.name, state]
                    state_counts[state] += mp
                    year2state_counts[year][state] += mp

            state_df = pd.DataFrame(index=states, columns=['samples'], data=[[state_counts[s]] for s in states])
            total_counts = state_df['samples'].sum()
            state_df['%'] = 100 * state_df['samples'] / total_counts
            state_df.to_excel(writer, sheet_name='{} tip states'.format(label), startrow=0, startcol=0, float_format='%.0f')

            transmission2counts = count_transmissions(tree, mp_df)

            df = pd.DataFrame(columns=states, index=states)

            for (from_state, to_state), count in transmission2counts.items():
                df.loc[from_state, to_state] = count

            total_transitions = df.sum().sum()
            for s in states:
                percentage_s = '% of {}'.format(s)
                df[percentage_s] = 0
                df.loc[percentage_s, percentage_s] = 0
            for state in states:
                df['% of {}'.format(state)] = 100 * df.loc[states, state] / df.loc[states, state].sum(skipna=True)
                df.loc['% of {}'.format(state), states] \
                    = 100 * df.loc[state, states] / df.loc[state, states].sum(skipna=True)
                df.loc[['% of {}'.format(s) for s in states], '% of {}'.format(state)] \
                    = (100 * df.loc[states, state] / total_transitions).tolist()

            print(df)

            df.to_excel(writer, sheet_name='{} transmissions'.format(label), startrow=0, startcol=0,
                        index_label='From \\ To', float_format='%.0f')

            counts = np.round(np.array(state_df['%'], dtype=float), 1)
            transitions = np.round(
                np.array(df.loc[['% of {}'.format(s) for s in states], ['% of {}'.format(s) for s in states]],
                         dtype=float), 0)
            save_as_transition_html(params.column, states, counts=counts,
                                    transitions=transitions,
                                    out_html=params.out_html.format(label, min_year, max_year),
                                    state2colour=state2color, work_dir=None,
                                    local_css_js=False, threshold=0)

            # year2count_array = {}
            # year2transmission_array = {}
            n = len(states)
            years = sorted(set((int(getattr(n, DATE))) // DATE_STEP for n in tree.traverse()))
            n_years = len(years)
            count_array = np.zeros(n * n_years, dtype=float)
            transmission_array = np.zeros((n * n_years, n * n_years), dtype=float)
            state_labels = []
            for year_i, year in enumerate(years):
                for s in states:
                    suffixed_s = '{}, {}s'.format(s, year * DATE_STEP)
                    state_labels.append(suffixed_s)
                    state2color[suffixed_s] = state2color[s]

                state_counts = year2state_counts[year]
                # total_counts = sum(state_counts.values())
                count_array[year_i * n: (year_i + 1) * n] = [state_counts[s] for s in states]
                # year2count_array[year] = np.array([state_counts[s] for s in states], dtype=float)

                y_transmission2counts = count_transmissions(tree, mp_df,
                                                            filter=lambda n: year == (int(getattr(n, DATE))) // DATE_STEP)
                a = np.zeros((n, n), dtype=float)
                for i in range(n):
                    i_state = states[i]
                    a[i, i] = y_transmission2counts[(i_state, i_state)]
                    for j in range(0, i):
                        j_state = states[j]
                        a[i, j] = y_transmission2counts[(i_state, j_state)]
                        a[j, i] = y_transmission2counts[(j_state, i_state)]
                # total_transitions = a.sum()
                transmission_array[year_i * n: (year_i + 1) * n, year_i * n: (year_i + 1) * n] = a
                # year2transmission_array[year] = a

            count_array = 100 * count_array / total_counts
            transmission_array = 100. * transmission_array / total_transitions

            counts = np.round(count_array, 1)
            transitions = np.round(transmission_array, 0)
            if np.any(transitions > 0):
                save_as_transition_html(params.column, state_labels, counts=counts,
                                        transitions=transitions,
                                        out_html=params.out_html.format(label, 'by', DATE_STEP),
                                        state2colour=state2color, work_dir=None,
                                        local_css_js=False, threshold=0)

            logging.info('Analysed who infected whom in {}'.format(label))