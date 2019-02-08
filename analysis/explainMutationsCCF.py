import os, sys
import argparse
import warnings
import itertools
from itertools import cycle
from collections import Counter
from collections import deque
from scipy.stats import beta
from scipy.stats import binom_test



def parse_args():
    description = "Compute CCF and mutated copies from given clone and allele-specific copy-number states, proportions, and somatic mutations."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("INPUT", type=str, help="CSV file of somatic mutations.")
    parser.add_argument("-s", "--seg", required=True, type=str, help="SEG.UCN file with inferred allele and clone-specific copy-number states and proportions.")
    parser.add_argument("-t", "--tool", required=False, default='HATCHet', type=str, help="Name of the method which inferred the copy-number states and proportions.")
    parser.add_argument("-e", "--explain", action='store_true', default=False, required=False, help="Compute whether mutations are explained by the copy numbers and proportions.")
    args = parser.parse_args()

    if not os.path.isfile(args.INPUT):
        raise ValueError("Mutation file does not exist!")
    if not os.path.isfile(args.seg):
        raise ValueError("SEG.UCN file does not exist!")

    return { "snv" : args.INPUT,
             "cns" : args.seg,
             "tool" : args.tool,
             "explain" : args.explain }


def main():
    log('Reading and parsing input arguments')
    args = parse_args()

    log('Reading allele and clone-specific copy numbers in SEG.UCN format')
    cns = read_cns(args['cns'])

    log('Reading somatic mutations and inferring mutated copies')
    snv = read_snv(args['snv'], cns, args['tool'], args['explain'])

    log('Formatting and writing')
    if not args['explain']:
        print '\t'.join(['#CHR', 'POS', 'PATIENT-SAMPLE', 'TOOL', 'COV', 'COUNTS', 'ObservedVAF', 'predicted_VAF', 'Error', 'MutatedCopies', 'CCF'])
        rec = (lambda s : '\t'.join(map(str, [s['chr'], s['pos'], s['Sample'], s['tool'], s['cov'], s['counts'], s['observed_VAF'], s['predicted_VAF'], s['Err'], s['mutated_copies'], s['CCF']])))
        print '\n'.join([rec(s) for s in snv])
    else:
        print '\t'.join(['#CHR', 'POS', 'PATIENT-SAMPLE', 'TOOL', 'COV', 'COUNTS', 'ObservedVAF', 'predicted_VAF', 'Error', 'MutatedCopies', 'CCF', 'Explained'])
        rec = (lambda s : '\t'.join(map(str, [s['chr'], s['pos'], s['Sample'], s['tool'], s['cov'], s['counts'], s['observed_VAF'], s['predicted_VAF'], s['Err'], s['mutated_copies'], s['CCF'], s['Explained']])))
        print '\n'.join([rec(s) for s in snv])




def read_cns(path):
    cns = {}
    with open(path, 'r') as f:
        for line in f:
            if len(line) > 1 and line[0] != '#':
                parsed = line.strip().split()
                chrs = parsed[0]
                seg = (int(parsed[1]), int(parsed[2]))
                sample = parsed[3]
                if sample not in cns:
                    cns[sample] = {}
                if chrs not in cns[sample]:
                    cns[sample][chrs] = {}
                assert seg not in cns[sample][chrs]
                cn = [tuple(map(int, e.split('|'))) for i,e in enumerate(parsed[4:]) if i%2==0]
                u = [float(e) for i,e in enumerate(parsed[4:]) if i%2==1]
                cns[sample][chrs][seg] = [e for e in zip(cn, u) if e[1] > 0.0]
                assert cns[sample][chrs][seg][0][0] == (1, 1)
                assert 0.99 <= sum(e[1] for e in cns[sample][chrs][seg]) <= 1.01
    return cns


def read_snv(path, cns, tool, explain):
    snv = []
    with open(path, 'r') as i:
        header = [t for t in i.readline().strip().split(',')]
        for line in i:
            row = {header[x] : r for x, r in enumerate(line.strip().split(','))}
            sel = True and ''.join([l for l in row['chrom'] if l.isdigit()]) in [str(c) for c in range(1, 23)]
            sample = None
            if '{}-{}'.format(row['Patient'], row['Sample']) in cns:
                sel = sel and True
                sample = '{}-{}'.format(row['Patient'], row['Sample'])
            elif '{}_{}'.format(row['Patient'], row['Sample']) in cns:
                sel = sel and True
                sample = '{}_{}'.format(row['Patient'], row['Sample'])
            else:
                sel = False
            sel = sel and 'Somatic' == row['somatic_status']
            if sel:
                c = str(row['chrom'])
                o = int(row['position'])

                find = findseg(cns, sample, c, o)
                assert len(find) <= 2, '{}\n{}'.format(o, find)

                if len(find) == 1:
                    find = find[0][-1]
                    snv.append(record(find, tool, row, sample, explain))
                elif len(find) == 2:
                    record1 = record(find[0][-1], tool, row, sample, explain)
                    record2 = record(find[1][-1], tool, row, sample, explain)
                    if record1['Err'] <= record2['Err']:
                        snv.append(record1)
                    else:
                        snv.append(record2)
    return snv


def record(f, tool, row, sample, explain):
    best, est, ccf = estimate(f, row)
    record = {}
    record['chr'] = row['chrom']
    record['pos'] = int(row['position'])
    record['Sample'] = sample
    record['tool'] = tool
    record['predicted_VAF'] = est
    record['CCF'] = ccf
    record['observed_VAF'] = get_perc(row['tumor_var_freq'])
    record['Err'] = abs(record['observed_VAF'] - est)
    record['cov'] = int(row['tumor_reads1']) + int(row['tumor_reads2'])
    record['counts'] = '{},{}'.format(int(row['tumor_reads1']), int(row['tumor_reads2']))
    record['mutated_copies'] = ','.join(map(str, tuple(best)))
    if explain:
        record['Explained'] = isconf((int(row['tumor_reads1']), int(row['tumor_reads2'])), est, 0.05)
    return record


def findseg(data, s, c, o):
    segs = sorted(list(data[s][c]), key=(lambda x : x[0]))
    L = 0
    R = len(segs) - 1

    if not (segs[0][0] <= o <= segs[-1][1]):
        return []
    if o == segs[L][0]:
        return [(s, c, segs[L], data[s][c][segs[L]])]
    if o == segs[R][1]:
        return [(s, c, segs[R], data[s][c][segs[R]])]

    while (R - L) > 1:
        M = int(round(float(R + L) / 2.0))
        if segs[L][0] <= o < segs[M][1]:
            R = M
        else:
            L = M
    if segs[L][0] <= o < segs[L][1]:
        return [(s, c, segs[L], data[s][c][segs[L]])]
    elif segs[R][0] < o < segs[R][1]:
        return [(s, c, segs[R], data[s][c][segs[R]])]
    elif o == segs[L][1] and o == segs[R][0]:
        return [(s, c, segs[L], data[s][c][segs[L]]), (s, c, segs[R], data[s][c][segs[R]])]
    elif o == segs[L][1]:
        return [(s, c, segs[L], data[s][c][segs[L]])]
    elif o == segs[R][0]:
        return [(s, c, segs[R], data[s][c][segs[R]])]
    else:
        if not(segs[L][1] < o < segs[R][0]):
            print o
            print L
            print segs[L]
            print R
            print segs[R]
        assert segs[L][1] < o < segs[R][0]
        return []


def get_perc(s):
    if '%' in s:
        return float(s.replace('%', '')) / 100.00
    else:
        return float(s)


def estimate(inf, row):
    assert inf[0][0] == (1, 1), str(inf)
    avail = [[0]] + [range(max(i[0])+1) for i in inf[1:]]
    num = (lambda y : sum(float(e) * inf[i][1] for i,e in enumerate(y)))
    den = sum(float(sum(i[0])) * i[1] for i in inf)
    compute = (lambda y : num(y) / den)
    ests = {x : compute(x) for x in itertools.product(*avail)}
    ests = {d : ests[d] for d in ests if 0.0 <= ests[d] <= 1.0}
    best = min(ests, key=(lambda x : abs(ests[x] - get_perc(row['tumor_var_freq']))))
    purity = sum(i[1] for i in inf[1:])
    assert 0.99 <= purity + inf[0][1] <= 1.01
    mut = sum(inf[i][1] for i, e in enumerate(best) if e > 0)
    assert mut <= purity
    ccf = float(mut) / float(purity)
    return best, ests[best], ccf


def isconf((countA, countB), est, gamma):
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    [c_lower, c_upper] = beta.ppf([p_lower, p_upper], countA + 1, countB + 1)
    return c_lower <= est <= c_upper or c_lower <= (1.0 - est) <= c_upper


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def argmax(d):
    return max(d, key=(lambda x : d[x]))


def argmin(d):
    return min(d, key=(lambda x : d[x]))


def log(M):
    sys.stderr.write('#' + M + '\n')


if __name__ == '__main__':
        main()
