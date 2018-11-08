import math
import argparse
import torch
import numpy as np


def get_strand_color(is_rev):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    is_rev = int(math.ceil(is_rev))
    if is_rev == 254:
        return 'R'
    if is_rev == 240:
        return '1'
    elif is_rev == 70:
        return '0'
    else:
        return ' '


def get_alt_type(alt_type_color):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    alt_type_color = int(math.ceil(alt_type_color))
    if alt_type_color == 0:
        return ' '
    elif alt_type_color == 5:
        return '0'
    elif alt_type_color == 240:
        return '1'
    elif alt_type_color == 125:
        return '2'
    elif alt_type_color == 254:
        return 'R'


def get_base_from_color(base_color):
    color = int(math.ceil(base_color))
    global_base_color_reverse = {250: 'A', 30: 'C', 180: 'G', 100: 'T', 5: 'N', 10: '.', 20: '*'}
    if color in global_base_color_reverse:
        return global_base_color_reverse[color]
    else:
        return ' '
    # 'A': 25.0, 'C': 75.0, 'G': 125.0, 'T': 175.0, '*': 225.0


def get_quality_by_color(quality):
    """
    Get a color spectrum given mapping quality
    :param map_quality: value of mapping quality
    :return:
    """
    quality = int(math.ceil(quality))
    color = math.floor(((quality / 254) * 9))
    if color == 0:
        return ' '
    return str(color)


def analyze_tensor(image):
    # base_color, base_quality_color, map_qual_color, strand_color, alt_color
    img_c, img_w, img_h = image.size()
    image = np.array(image.data * 254)
    img_h = 50
    # label_c = label.size(0)
    # print()
    # for i in range(label_c):
    #     print(label[i].item(), end='')
    # print()
    print("BASE CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
                print(get_base_from_color(image[0][j][i]), end='')
        print()

    print("BASE QUALITY CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[1][j][i]), end='')
        print()
    print("MAPPING QUALITY CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[2][j][i]), end='')
        print()
    print("STRAND DIRECTION CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_strand_color(image[3][j][i]), end='')
        print()
    print("ALT FREQ CHANNEL:")
    for i in range(img_h):
        for j in range(img_w):
            print(get_alt_type(image[4][j][i]), end='')
        print()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tensor_file",
        type=str,
        required=True,
        help="Training data description csv file."
    )
    parser.add_argument(
        "--label_file",
        type=str,
        required=True,
        help="Training data description csv file."
    )
    FLAGS, unparsed = parser.parse_known_args()
    if FLAGS.label_file:
        labels = torch.load(FLAGS.label_file)
        for i in range(0, labels.size(0)):
            print(labels[i].item(), end='')
        print()

    if FLAGS.tensor_file:
        image = torch.load(FLAGS.tensor_file)
        analyze_tensor(image)
