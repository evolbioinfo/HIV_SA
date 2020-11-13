from pastml.tree import name_tree, DATE, DATE_CI, read_forest

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_tree', required=True, type=str)
    parser.add_argument('--output_tree', required=True, type=str)
    params = parser.parse_args()

    tr = read_forest(params.input_tree)[0]
    name_tree(tr)

    tr.write(outfile=params.output_tree, format_root_node=True, format=3, features=[DATE, DATE_CI])
