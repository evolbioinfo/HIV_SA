import datetime
import logging
import os
import re

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from pastml import numeric2datetime, datetime2numeric

EXTERNAL = 'External'

MEDIUM = 'Medium'

LOW = 'Low'

HIGH = 'High'

NAME = 'sequencename'
DATE = 'sampledate'
LOCATION = 'urbanrural'
PREVALENCE = 'prevalence'
HLE = 'highlow_prevalence'
HLME = 'highmedlow_prevalence'


def format_id(id):
    dot_year = '\.\d\d\d\d$'
    underscore_date = '_\d\d\d\d-\d\d-\d\d$'
    date = re.findall(dot_year, id)
    if date:
        id = re.sub(dot_year, '', id)
        date = int(date[0][1:])
    else:
        date = re.findall(underscore_date, id)
        if date:
            date = datetime.datetime.strptime(date[0][1:], '%Y-%m-%d')
            id = re.sub(underscore_date, '', id)
        else:
            date = None
    return id.replace('.', '_'), date


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    import argparse

    parser = argparse.ArgumentParser(description="Processes data files.")

    parser.add_argument('--data_in', required=True, type=str, help="the annotation file.")
    parser.add_argument('--sequences_in', required=True, type=str, help="the fasta file.")
    parser.add_argument('--sequences_out', required=True, type=str, help="the out fasta file.")
    parser.add_argument('--to_remove', required=True, type=str, help="the sequences to skip.")
    parser.add_argument('--data_out', required=True, type=str, help="the sequence annotation file.")
    parser.add_argument('--dates', required=True, type=str, help="the date annotation file.")
    params = parser.parse_args()

    # Read and fix metadata
    df = pd.io.stata.read_stata(params.data_in)
    df.index = df['id'].apply(lambda _: format_id(_)[0])
    df.drop(labels=['id'], axis=1, inplace=True)
    df[DATE] = pd.to_datetime(df[DATE], yearfirst=False, format="%d/%m/%Y")
    for cat in (HLE, HLME):
        df[cat].replace('LastVisit', EXTERNAL, inplace=True)
    df[HLME] = df[HLME].apply(lambda _: '' if _ not in {HIGH, LOW, MEDIUM, EXTERNAL} else _)
    df[HLE] = df[HLE].apply(lambda _: '' if _ not in {HIGH, LOW, EXTERNAL} else _)
    df[LOCATION].replace('.', '', inplace=True)

    df[DATE] = pd.to_datetime(df[DATE], dayfirst=True, yearfirst=False, infer_datetime_format=True)
    logging.info(df.head())

    to_remove = set()
    if params.to_remove:
        with open(params.to_remove, 'r') as f:
            to_remove = set(f.read().strip().strip('\n').split('\n'))

    ids = []
    with open(params.sequences_out, 'w+') as f:
        for rec in SeqIO.parse(params.sequences_in, 'fasta', alphabet=generic_dna):
            rec.id, date = format_id(rec.id)
            if rec.id in to_remove:
                continue
            if rec.id not in df.index:
                df.loc[rec.id, [DATE, HLE, HLME]] = date, EXTERNAL, EXTERNAL
            ids.append(rec.id)
            f.write('>{}\n{}\n'.format(rec.id, str(rec.seq).replace('_', '-')))

    df = df.loc[ids, :]
    df['lsdate'] = df[DATE].apply(lambda _: datetime2numeric(_) if isinstance(_, datetime.datetime)
        else ('b({},{})'.format(_, _ + 1) if isinstance(_, int) else 'u(2016)'))
    with open(params.dates, 'w+') as f:
        f.write('{}\n'.format(len(df)))
    df['lsdate'].to_csv(params.dates, sep='\t', header=False, mode='a')
    df[DATE] = df[DATE].apply(lambda _: _.strftime('%Y-%m-%d') if isinstance(_, datetime.datetime) else _)
    df.to_csv(params.data_out, sep='\t', index_label='id')
