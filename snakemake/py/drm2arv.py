import logging
import re

import pandas as pd
import wikipedia

from sierrapy import SierraClient

QUERY = '''
drugResistance {
    gene { name },
    drugScores {
        drugClass { name },
        drug { displayAbbr, fullName },
        score,
        text,
        partialScores {
            mutations {
                text,
            },
        },
    }
}
'''


def get_date(drug):
    summary = wikipedia.page(wikipedia.search('{} HIV'.format(drug))[0]).summary
    dates = set()
    for date in re.findall(r'\s[12][901][8901]\d[^\d]', summary):
        index = summary.find(date)
        prev_sentence = summary[max(0, index - 50): index]
        if 'approv' in prev_sentence or 'sold ' in prev_sentence:
            dates.add(date[1:-1])
    return min(dates) if dates else None


if '__main__' == __name__:
    import argparse
    parser = argparse.ArgumentParser(description="Extracts SDRM drug resistance information.")
    parser.add_argument('--drms', nargs='+', type=str, help="SDRMs of interest")
    parser.add_argument('--output', required=True, type=str, help="output file in tab-delimited format.")
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    params.drm = [_.replace('PR_', 'PR:').replace('RT_', 'RT:') for _ in params.drms]

    data = []

    for dr in SierraClient().mutations_analysis(params.drms, QUERY)['drugResistance']:
        gene = dr['gene']['name']
        drug_scores = dr['drugScores']
        for ds in drug_scores:
            text = ds['text']
            if 'Resistance' not in text:
                continue
            drug_class = ds['drugClass']['name']
            drug_abbr = ds['drug']['displayAbbr'].replace('/r', '')
            drug_name = ds['drug']['fullName'].replace('/r', '')
            score = ds['score']
            partial_scores = ds['partialScores']
            for ps in partial_scores:
                for _ in ps['mutations']:
                    drm = _['text']
                    data.append(['{}:{}'.format(gene, drm), drug_class, drug_name, drug_abbr, score, text])
    df = pd.DataFrame(data, columns=['mutation', 'drug class', 'drug full name', 'drug abbreviation', 'score', 'note'])
    for drug in df['drug full name'].unique():
        df.loc[df['drug full name'] == drug, 'year'] = get_date(drug)
    df.to_csv(params.output, sep='\t', index=False)


