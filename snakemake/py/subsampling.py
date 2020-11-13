from collections import defaultdict

import numpy as np
import pandas as pd
from pastml.tree import read_tree, remove_certain_leaves, annotate_dates, DATE, DATE_CI

EXT_COLOR = '#4daf4a'
HIGH_COLOR = '#e41a1c'
LOW_COLOR = '#377eb8'

state2colour = {'High': HIGH_COLOR, 'Low': LOW_COLOR, 'External': EXT_COLOR}
column = 'highlow_prevalence'


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--tree', required=True, type=str,
                        help='Input tree with tip annotated with states and dates')
    parser.add_argument('--subtree', type=str, required=True)
    params = parser.parse_args()

    tree = read_tree(params.tree, columns=[column])
    annotate_dates([tree])
    state2tips = defaultdict(list)
    state2year2tips = defaultdict(lambda: defaultdict(list))
    tip2date = {}

    min_year, max_year = np.inf, -np.inf
    for tip in tree:
        states = getattr(tip, column)
        date = getattr(tip, DATE)
        tip2date[tip.name] = date
        if len(states) == 1:
            state = next(iter(states))
            state2tips[state].append(tip.name)
            year = int(date)
            if year < min_year:
                min_year = year
            if year > max_year:
                max_year = year
            state2year2tips[state][year].append(tip.name)
    total = len(tree)

    size = min(len(_) for _ in state2tips.values())
    min_ac_year = min(min(state2year2tips[_].keys()) for _ in ['Low', 'High'])

    print('Gonna subsample {} sequences of each state'.format(size))
    print('Sampling in AC started in {}'.format(min_ac_year))

    subsampled_ids = set()
    for state, year2tips in state2year2tips.items():
        total_tips = len(state2tips[state])
        df = pd.DataFrame(index=list(year2tips.keys()),
                          data=[[len(v)] for v in year2tips.values()],
                          columns=['n'])
        df.sort_values(by=['n'], inplace=True, ascending=True)
        left = size

        if state == 'External':
            sampled_not_before_ac = sum(df.loc[df.index >= min_ac_year, 'n'])
            if sampled_not_before_ac >= size:
                years = [y for y in df.index if y >= min_ac_year]
                df = df.loc[years, :]
            else:
                for y, tips in year2tips.items():
                    if y >= min_ac_year:
                        subsampled_ids |= set(tips)
                years = [y for y in df.index if y < min_ac_year]
                df = df.loc[years, :]
                left -= sampled_not_before_ac

        # if there are not enough sequences sampled (less than we want to keep),
        # we are going to take a bit more of other years to compensate for it
        n = len(df)
        for i, (year, row) in enumerate(df.iterrows()):
            avail = row['n']
            like_to_take = int(np.round(left / (n - i), 0))
            can_take = min(like_to_take, avail)
            if can_take:
                left -= can_take
                newly_added_tips = set(np.random.choice(year2tips[year], size=can_take, replace=False))
                subsampled_ids |= newly_added_tips

    tree = remove_certain_leaves(tree, lambda _: _.name not in subsampled_ids)
    tree.dist = 0
    tree.write(outfile=params.subtree, format=3, format_root_node=True, features=[DATE, DATE_CI])