from argparse import ArgumentParser
import os
import copy
import gzip
import logging
from collections import OrderedDict, defaultdict, Counter
from singlecell.utils import file_transaction

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)


class readsInWell:
    def __init__(self, g, d):
        self.gene = g
        self.umi = d


def get_umi_well(in_file):
    umi_well = {}
    with gzip.open(in_file, 'rb') as in_handle:
        for l in in_handle:
            cols = l.strip().split("\t")
            read = (cols[1], cols[2])
            if read not in umi_well:
                umi_well[read] = readsInWell(cols[3], {cols[0]: int(cols[4])})
            else:
                umi_well[read].umi.update({cols[0]: int(cols[4])})
    return umi_well


def write_summary(umi_well):
    """write summary about edit distance among same read position"""
    for read in umi_well:
        umi_list = Counter(umi_well[read].umi)
        ma = calculate_matrix_distance(umi_list.keys())
        for pair in ma:
            max_umi, min_umi = pair[1], pair[0]
            if umi_list[pair[0]] > umi_list[pair[1]]:
                max_umi, min_umi = pair[0], pair[1]
            print "%s" % (" ".join(map(str, [read[0], read[1],
                                       umi_list[max_umi], umi_list[min_umi],
                                       ma[pair]])))


def merge_umis(umi_well):
    """merge umis with 1-edit-distance on the same genome position"""
    for read in umi_well:
        umi_list = Counter(umi_well[read].umi)
        logger.debug("merge_umis: popular %s" % [read[1], umi_well[read].gene])
        logger.debug("related umi to a position: %s" % umi_list.keys())
        ma = calculate_matrix_distance(umi_list.keys())
        prob = most_voted_umi(ma, umi_list)
        if prob:
            umi_well[read].umi = calc_num_umi(prob, ma, umi_list)
            #logger.debug("%s" % " ".join(map(str, [read[1], umi_well[read].gene,
            #                                        popular, max_count, umi_list[max_count],
            #                                        min_count, umi_list[min_count]])))
    umi_well[read].umi = umi_list
    return umi_well


def get_order(pairs, popular):
    """return ordered umis according to counts"""
    if pairs[0] in popular:
        return pairs[0], pairs[1]
    else:
        return pairs[1], pairs[0]


def most_voted_umi(ma, umi):
    """get the matrix distance and select the umi most vote
    to merge the rest with it"""
    if len(ma) > 0:
        iamerror, iamreal = Counter(), Counter()
        for pair in ma:
            if ma[pair] == 9:
                for p in pair:
                    if umi[p] == 1:
                        iamerror[p] += 1
                    else:
                        iamreal[p] += 1
        logger.debug("most_voted_umi: real %s" % iamreal)
        logger.debug("most_voted_umi: error %s" % iamerror)
        prob = get_error_prob(iamerror, iamreal)
        if len(prob) > 0:
            logger.debug("most_voted_umi: probabilities %s" % prob)
            prob = get_decision(prob, ma, umi)
            logger.debug("most_voted_umi: final probabilities %s" % prob)
            return prob


def get_error_prob(e, r):
    """get a probability based on e/(e+r) for each umi"""
    prob = defaultdict(float)
    for u in list(set(e.keys()).union(r.keys())):
        prob[u] = 1.0 * e[u]/(e[u]+r[u])
    return prob


def get_decision(p, ma, umi, need_decide=False):
    """decide real or error"""
    num_umi = calc_num_umi(p, ma, umi)
    logger.debug("decide: umis %s" % num_umi)
    if len([v for v in p.values() if v != 1 and v != 0]) > 0:
        need_decide = True
    logger.debug("decide: need_decide %s" % need_decide)
    while need_decide:
        unknown = [k for k, v in p.iteritems() if v < 1 or v > 0]
        logger.debug("decide: unknown %s" % p)
        p_fake = copy.deepcopy(p)
        p_fake[unknown[0]] = 1
        num_umi_fake_real = len(calc_num_umi(p_fake).keys())
        p_fake[unknown[0]] = 0
        num_umi_fake_error = len(calc_num_umi(p_fake).keys())
        if num_umi_fake_error > num_umi_fake_real:
            p[unknown[0]] = 1
        else:
            p[unknown[0]] = 0
        logger.debug("decide: assigned %s" % p)
        if len([v for v in p.values() if v != 1 and v != 0]) == 0:
            need_decide = False
        logger.debug("decide: need_decide %s" % need_decide)
    return p


def calc_num_umi(p, ma, umi):
    for pair in ma:
        if p[pair[0]] == 0 and p[pair[1]] == 1:
            umi[pair[0]] += umi[pair[1]]
            del umi[pair[1]]
            logger.debug("calc_num_umi: %s" % [pair[0], umi[pair[0]], pair[1], umi[pair[1]]])
        elif p[pair[1]] == 0 and p[pair[0]] == 1:
            umi[pair[1]] += umi[pair[0]]
            del umi[pair[0]]
            logger.debug("calc_num_umi: %s" % [pair[1], umi[pair[1]], pair[0], umi[pair[0]]])
    return umi


def calculate_matrix_distance(umis_pos):
    """find putative similar UMIs"""
    ma = {}
    for umi1 in umis_pos:
        [ma.update({(umi1, umi2): distance(umi1, umi2)}) for umi2 in umis_pos if umi1 != umi2 and (umi2, umi1) not in ma]
    logger.debug("matrix %s" % ma)
    return ma


def distance(u1, u2):
    """calculate distance"""
    score = sum([1 for nt1, nt2 in zip(u1, u2) if nt1 == nt2])
    logger.debug("distance: %s" % [u1, u2, score])
    return score


def gene_counts(umi_well):
    """summarize by gene and umi and well"""
    counts = defaultdict(Counter) 
    for read in umi_well:
        for umi in umi_well[read].umi.keys():
            counts[(umi_well[read].gene, read[1])] = umi_well[read].umi
    return counts


def write_extensive_summary(well_umi_gen, out_file):
    with file_transaction(out_file) as tx_out_file:
        with gzip.open(tx_out_file, 'wb') as out_handle:
            well_umi_gen_str = [[("\t%s\t%s\t" % (gen_well[0], gen_well[1])).join(map(str, umi)) for umi in well_umi_gen[gen_well].items()] for gen_well in well_umi_gen]
            out_handle.write("\n".join(["\n".join(item) for item in well_umi_gen_str]))
            out_handle.write("\n")


if __name__ == "__main__":
    parser = ArgumentParser(description="Get UMIs stats")
    parser.add_argument("--counts-umi", required=True,
                        help="file from umi_stats.py")
    parser.add_argument("--only-stats",
                        help="stats about edit distance among UMIs from same position",
                        action='store_true')
    args = parser.parse_args()
    umi_well = get_umi_well(args.counts_umi)
    if args.only_stats:
        write_summary(umi_well)
    else:
        umi_well = merge_umis(umi_well)
        well_umi_gen = gene_counts(umi_well)
        base, _ = os.path.splitext(args.counts_umi)
        write_extensive_summary(well_umi_gen, base + "_reduced.gz")
