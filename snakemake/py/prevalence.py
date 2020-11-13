import pandas as pd


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    parser.add_argument('--subtype', required=False, type=str, default=None)
    parser.add_argument('--subtype_col', default='Sierra subtype', type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.input, sep='\t', index_col=0)
    if params.subtype and params.subtype_col:
        df = df[df[params.subtype_col] == params.subtype]
    columns = [col for col in df.columns if 'RT:' in col or 'PR:' in col or 'IN:' in col]
    df = df[columns]
    ((df == 'resistant').astype(int).sum().sort_values() / len(df)).to_csv(params.output, header=False, index=True, sep='\t')
