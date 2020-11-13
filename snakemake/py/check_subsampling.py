from collections import Counter

import pandas as pd

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tab', required=True, type=str)
    parser.add_argument('--output_log', required=True, type=str)
    parser.add_argument('--column', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.input_tab, sep='\t', header=0, index_col=0)
    interesting_columns = [c for c in df.columns if params.column in c]
    df = df[interesting_columns]
    values = sorted([_ for _ in df[params.column].unique() if not pd.isna(_)])
    ids = df.index.unique()
    stats_df = pd.DataFrame(index=ids,
                            columns=['state_union', 'state_intersection', 'number_subtrees', 'full_tree_state']
                                    + values)

    def get_value_set(id, column=params.column, many_values=True):
        id_df = df.loc[id, column]
        return {_ for _ in id_df.unique() if not pd.isna(_)} \
            if many_values else ({id_df} if not pd.isna(id_df) else set())

    agreed = 0
    majority = 0
    seen = 0
    different = 0
    total = 0
    for id in ids:
        many_values = len(df.loc[id, :].shape) > 1
        ref_values = get_value_set(id, many_values=many_values)
        value_intersection = set(values)
        value_union = set()
        value_counter = Counter()
        n = 0
        for col in df.columns:
            if col != params.column:
                col_values = get_value_set(id, column=col, many_values=many_values)
                if col_values:
                    value_counter.update(col_values)
                    value_intersection &= col_values
                    value_union |= col_values
                    n += 1
        stats_df.loc[id, :] = ', '.join(sorted(value_union)), ', '.join(sorted(value_intersection) if n else []), n, \
                              ', '.join(sorted(ref_values)), *[value_counter[v] for v in values]
        is_internal_id = id.startswith('n') or id == 'root'
        if not n:
            continue
        if is_internal_id:
            total += 1
        if not value_intersection & ref_values:
            most_common_count = value_counter.most_common(1)[0][1]
            most_common_values = set(_ for _ in values if value_counter[_] == most_common_count)
            if not most_common_values & ref_values:
                if value_union & ref_values:
                    seen += 1
                else:
                    different += 1
                    print(id, n, value_counter, ref_values)
            elif is_internal_id:
                majority += 1
        elif is_internal_id:
            agreed += 1

    stats_df.to_csv(params.output_log, sep='\t', index_label='id')

    print('Out of {} internal nodes:'
          '\t{} ({:.1f}%) have a common prediction in (sub)trees'
          '\t{} ({:.1f}%) have a full tree prediction corresponding to the majority prediction in the subtrees'
          '\t{} ({:.1f}%) have a full tree prediction that was also predicted in some of the subtrees'
          '\t{} ({:.1f}%) have a different prediction in the full trees'
          .format(total,
                  agreed, 100 * agreed / total,
                  majority, 100 * majority / total,
                  seen, 100 * seen / total,
                  different, 100 * different / total,
                  ))
