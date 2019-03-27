import h5py
import argparse
from collections import defaultdict, OrderedDict
import yaml
from tqdm import tqdm

dv_dict = defaultdict(tuple)


def read_dv_candidates(dv_file):
    with open(dv_file, 'r') as tsv:
        chr_name, pos_start, pos_end, ref, alts, gt = None, None, None, None, None, None

        for i, line in enumerate(tsv):
            if (i + 1) % 3 == 1:
                chr_name, pos_start, pos_end, ref = line.rstrip().split('\t')
            if (i + 1) % 3 == 2:
                alts = line.rstrip().split('\t')
            if (i + 1) % 3 == 0:
                gt = line.rstrip().split('\t')
                if chr_name is None or pos_start is None or \
                        pos_end is None or ref is None or alts is None or gt is None:
                    print("DV CANDIDATE ERROR IN LINE: ", i)
                    print(line)
                    exit(1)
                dv_dict[(chr_name, pos_start, pos_end)] = (ref, alts, gt)
                chr_name, pos_start, pos_end, ref, alts, gt = None, None, None, None, None, None


def read_friday_candidates(friday_file):
    with open(friday_file, 'r') as tsv:
        chr_name, pos_start, pos_end, ref, alts, gt = None, None, None, None, None, None

        for i, line in enumerate(tsv):
            if (i + 1) % 3 == 1:
                chr_name, pos_start, pos_end, ref = line.rstrip().split('\t')
            if (i + 1) % 3 == 2:
                alts = line.rstrip().split('\t')
            if (i + 1) % 3 == 0:
                gt = line.rstrip().split('\t')
                if chr_name is None or pos_start is None or \
                        pos_end is None or ref is None or alts is None or gt is None:
                    print("FRIDAY CANDIDATE ERROR IN LINE: ", i)
                    print(line)
                    exit(1)
                if (chr_name, pos_start, pos_end) in dv_dict:
                    dv_ref, dv_alts, dv_gt = dv_dict[(chr_name, pos_start, pos_end)]
                    matched = True
                    if dv_ref != ref:
                        matched = False

                    for alt in alts:
                        if alt not in dv_alts:
                            matched = False

                    dv_gt_list = []
                    for gti in dv_gt:
                        if int(gti) == 0:
                            dv_gt_list.append(dv_ref)
                        else:
                            dv_gt_list.append(dv_alts[int(gti) - 1])

                    gt_list = []
                    for gti in gt:
                        if int(gti) == 0:
                            gt_list.append(ref)
                        else:
                            gt_list.append(alts[int(gti) - 1])

                    dv_gt_list = sorted(dv_gt_list)
                    gt_list = sorted(gt_list)

                    if dv_gt_list != gt_list:
                        matched = False

                    if not matched:
                        print("CANDIDATES DID NOT MATCH IN: ", ref, pos_start, pos_end)
                        print("DV CANDIDATES:\t", dv_ref, dv_alts, dv_gt)
                        print("FR CANDIDATES:\t", ref, alts, gt)
                else:
                    print("NO CANDIDATE FOUND BY DV IN: ", chr_name, pos_start, pos_end, ref, alts, gt)

                chr_name, pos_start, pos_end, ref, alts, gt = None, None, None, None, None, None

if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dv_candidates",
        type=str,
        required=True,
        help="DV candidate list."
    )
    parser.add_argument(
        "--friday_candidates",
        type=str,
        required=True,
        help="friday candidate list."
    )

    FLAGS, unparsed = parser.parse_known_args()
    read_dv_candidates(FLAGS.dv_candidates)
    read_friday_candidates(FLAGS.friday_candidates)

