import pandas as pd
from ete3 import Tree


def read_tree(nwk):
    tree = None
    for format in range(10):
        try:
            tree = Tree(nwk, format=format)
            break
        except:
            continue
    return tree


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tab', required=True, type=str)
    parser.add_argument('--input_acr', required=True, type=str)
    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tab', required=True, type=str)
    parser.add_argument('--arv_tab', required=True, type=str)
    parser.add_argument('--root_date', required=True, type=float)
    parser.add_argument('--arv', required=True, type=str)
    params = parser.parse_args()

    tree = read_tree(params.input_tree)
    for n in tree.traverse('preorder'):
        n.add_feature('date', params.root_date if n.is_root() else (getattr(n.up, 'date') + n.dist))

    df = pd.read_csv(params.input_tab, index_col=0, sep='\t')
    df = df[[params.arv]]
    df.index = df.index.map(str)
    df = df.loc[[_.name for _ in tree], :]

    acr_df = pd.read_csv(params.input_acr, index_col=0, sep='\t')
    acr_df.columns = df.columns
    acr_df.index = acr_df.index.map(str)

    df.drop(set(acr_df.index) & set(df.index), axis=0, inplace=True)

    arv_df = pd.read_csv(params.arv_tab, sep='\t')
    if params.arv in arv_df['mutation'].unique():
        arv_df = arv_df[arv_df['mutation'] == params.arv]
    else:
        arv_df = arv_df[arv_df['drug abbreviation'] == params.arv]
    drm_date = float(arv_df['year'].min())

    for n in tree.traverse():
        date = getattr(n, 'date')
        if date < drm_date:
            df.loc[n.name, params.arv] = 'sensitive'

    pd.concat([df, acr_df]).to_csv(params.output_tab, sep='\t', index_label='node id')

