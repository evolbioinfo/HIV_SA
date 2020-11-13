import pandas as pd
from ete3 import Tree, TreeNode


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

    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_forest', required=True, type=str)
    parser.add_argument('--root_date', required=True, type=float)
    parser.add_argument('--arv_tab', required=True, type=str)
    parser.add_argument('--arv', required=False, type=str, default=None)
    params = parser.parse_args()

    tree = read_tree(params.input_tree)
    arv_df = pd.read_csv(params.arv_tab, index_col=None, sep='\t')
    if params.arv in arv_df['mutation'].unique():
        arv_df = arv_df[arv_df['mutation'] == params.arv]
        message = 'First {}-provoking ARV'.format(params.arv)
    else:
        arv_df = arv_df[arv_df['drug abbreviation'] == params.arv]
        message = params.arv
    arv_year = float(arv_df['year'].min())
    print('{} was accepted in {}.'.format(message, arv_year))
    print('Root year is {}.'.format(params.root_date))

    nwks = []
    todo = [(tree, params.root_date)]
    while todo:
        node, date = todo.pop()
        if date < arv_year:
            for child in node.children:
                todo.append((child, date + child.dist))
        else:
            fake_root = TreeNode(dist=0, name='sensitive')
            fake_root.add_child(node, dist=date - arv_year)
            nwks.append(fake_root.write(format=3, format_root_node=True))
    with open(params.output_forest, 'w+') as f:
        f.write('\n'.join(nwks))
