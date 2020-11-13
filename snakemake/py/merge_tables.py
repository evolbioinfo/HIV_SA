import pandas as pd
from pastml.annotation import preannotate_forest
from pastml.tree import read_tree
from pastml.acr import _serialize_predicted_states, annotate_dates

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tabs', nargs='+', type=str)
    parser.add_argument('--input_names', nargs='+', type=str)
    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    parser.add_argument('--root_date', required=False, type=float, default=0)
    params = parser.parse_args()

    tree = read_tree(params.tree)
    annotate_dates([tree], root_dates=[params.root_date])
    columns = []
    for tab, name in zip(params.input_tabs, params.input_names):
        tab_df = pd.read_csv(tab, sep='\t', header=0, index_col=0)
        tab_df.columns = [name]
        columns.extend(preannotate_forest(forest=[tree], df=tab_df)[0])
    _serialize_predicted_states(columns, params.output_tab, [tree])

