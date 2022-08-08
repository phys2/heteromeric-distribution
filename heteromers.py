#!/usr/bin/env python
"""
Calculate how proteins distribute over possible homo- and
heteromeric populations based on protein abundances from
AP-MS data (affinity-purification mass-spectrometry).

Note:
  This code is not very generic (yet) and is
  still tied to the specific experimental setup.
"""
import argparse
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

#
# Experimental setup and data
#

# Targets
c1 = 0
c4 = 1
c5 = 2
targets = [c1, c4, c5]

# Combinations of protein assemblies
target_combinations = [
    [c1],
    [c4],
    [c5],
    [c1, c4],
    [c1, c5],
    [c4, c5],
    [c1, c4, c5],
]

# Series of APs with distinct order
series = [
    [c1, c4, c5],
    [c1, c5, c4],
    [c4, c1, c5],
    [c4, c5, c1],
    [c5, c1, c4],
    [c5, c4, c1],
]

# Target names (for export only)
target_names = {
    c1: 'C1',
    c4: 'C4',
    c5: 'C5',
}


def get_efficiencies(abundances: np.ndarray) -> np.ndarray:
    """
    Calculate AP efficiencies from abundances.

    :param abundances: AP abundance matrix

    :return: Matrix of 2 efficiency values per serie
    """
    eff = []
    for i, s in enumerate(series):
        # first AP in series
        t = s[0]
        k = i * len(targets) + t
        e0 = abundances[k, 0] / sum(abundances[k])

        # second AP in series
        t = s[1]
        k = i * len(targets) + t
        e1 = abundances[k, 1] / (abundances[k, 1] + abundances[k, 2])

        eff.append([e0, e1])

    return eff


def build_linear_systems(abundances: np.ndarray, colnames: [str]) -> dict[str, pd.DataFrame]:
    """
    Build input matrices for solver (system of linear equations).

    :param abundances: AP abundance matrix
    :param colnames: Column names for created DataFrames

    :return: Dict mapping each target name to its linear system dataframe
    """
    abundance_totals = np.sum(abundances, axis=1)
    efficiencies = get_efficiencies(abundances)

    results = {}
    for mti, target in enumerate(targets):
        rows = len(series) * len(targets)
        cols = len(target_combinations) + 1
        m = np.zeros((rows, cols))

        for si in range(len(series)):
            for ai in range(len(targets)):
                # first row index of current series in abundance matrix
                k = si * len(targets)
                # row index into output matrix
                i = k + ai

                # col 0 contains 'total' for current target in current series/AP
                m[i, 0] = abundances[k + mti, ai]

                # loop over homo/heteromere combinations
                for j, comb in enumerate(target_combinations):
                    # skip combination not containing current (matrix) target
                    if target not in comb:
                        continue

                    # skip combinations not containing current AP target
                    if series[si][ai] not in comb:
                        continue

                    a = abundance_totals[k + mti]

                    # factor in any previous AP's efficiency
                    for p in range(ai):
                        if series[si][p] in comb:
                            a *= 1 - efficiencies[si][p]

                    # factor in own AP's efficiency
                    if ai < 2:
                        a *= efficiencies[si][ai]

                    m[i, j + 1] = a

        name = target_names[target]
        results[name] = pd.DataFrame(data=m, columns=colnames).round().astype('int')

    return results


def save_linear_systems(systems: Mapping[str, pd.DataFrame], outfile: Path) -> None:
    """
    Save system of linear equations to a file.

    :param systems: Target names mapped to dataframe holding linear system
    :param outfile: Output file name
    """
    if outfile.suffix == '.xlsx':
        # this may raise Exception when Pandas is missing openpyxl lib
        with pd.ExcelWriter(outfile) as writer:
            for target, df in systems.items():
                df.to_excel(writer, sheet_name=target, index=False)
    else:
        sep = ',' if outfile.suffix == '.csv' else '\t'
        for target, df in systems.items():
            filename = outfile.with_stem(outfile.stem + '_' + target)
            df.to_csv(filename, index=False, sep=sep)


def solve_one(name: str, df: pd.DataFrame, verbose=False) -> pd.DataFrame:
    """
    Solve single linear equation system.

    :param name: Dataset name
    :param df: Equation matrix (sum in first column)
    :param verbose: If True, print solver output
    """
    # ignore all-zero columns
    colmask = df.sum() > 0
    df = df.loc[:, colmask]

    m = df.to_numpy()
    A = m[:, 1:]
    b = m[:, 0]
    res = lsq_linear(A, b, bounds=(0, 1), lsmr_tol='auto', verbose=0)

    if verbose:
        print('---')
        print(name)
        print(res)

    return pd.DataFrame(
        data=[res.x],  # noqa
        index=[name],
        columns=df.columns[1:]
    )


def solve_many(dataframes: Mapping[str, pd.DataFrame], verbose=False) -> pd.DataFrame:
    """Solve multiple input frames."""
    single_results = [
        solve_one(name, df, verbose=verbose)
        for name, df in dataframes.items()
    ]

    df = pd.concat(single_results, axis=0)
    # set very small values to zero
    df[df < 1E-6] = 0.0
    # sort columns
    cols = sorted(df.columns, key=lambda c: (len(c), c))
    df = df[cols]
    # normalize to 100%
    df = df.div(df.sum(axis=1), axis=0)
    df['Sum'] = df.sum(axis=1)

    return df


def main(file_abundances: Path, outfile_systems=Optional[Path],
         outfile=Optional[Path], verbose=False):
    # columns for solver input data
    colnames = ['total'] + [
        '+'.join([target_names[target] for target in comb])
        for comb in target_combinations
    ]

    abundances = np.loadtxt(file_abundances)

    systems = build_linear_systems(abundances, colnames)
    if outfile_systems:
        save_linear_systems(systems, outfile_systems)

    df_res = solve_many(systems, verbose=verbose).fillna('')
    if outfile:
        if outfile.suffix == '.xlsx':
            df_res.to_excel(outfile)
        else:
            sep = ',' if outfile.suffix == '.csv' else '\t'
            df_res.to_csv(outfile, sep=sep)
    else:
        with pd.option_context('display.float_format', '{:.1%}'.format):
            print(df_res)


def arguments():
    argp = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    argp.add_argument('abundances', help='Abundance data file')
    argp.add_argument('-s', '--out-systems', type=Path, metavar='FILE',
                      help='Save linear equations to this file or set of files.')
    argp.add_argument('-o', '--out', type=Path, metavar='FILE',
                      help='Save results to this file.')
    argp.add_argument('-v', '--verbose', action='store_true',
                      help='Print solver output.')

    argp.description = __doc__.split('.')[0]
    argp.epilog = (
        "Use '.xlsx' extension for output files to create Excel file "
        "(requires openpyxl library)."
    )

    return argp.parse_args()


if __name__ == '__main__':
    args = arguments()
    main(
        file_abundances=args.abundances,
        outfile_systems=args.out_systems,
        outfile=args.out,
        verbose=args.verbose,
    )
