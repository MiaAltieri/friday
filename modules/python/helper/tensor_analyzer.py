import math
def get_strand_color(is_rev):
    """
    Get color for forward and reverse reads
    :param is_rev: True if read is reversed
    :return:
    """
    is_rev = is_rev.item()
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
    alt_type_color = alt_type_color.item()
    if alt_type_color == 0:
        return ' '
    elif alt_type_color == 240:
        return '1'
    elif alt_type_color == 125:
        return '2'
    elif alt_type_color == 254:
        return 'R'


def get_base_from_color(base_color):
    color = base_color.item()
    global_base_color_reverse = {250: 'A', 30: 'C', 180: 'G', 100: 'T', 0: ' ', 10: 'N'}
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
    quality = quality.item()
    color = math.floor(((quality / 254) * 9))
    if color == 0:
        return ' '
    return str(color)


def analyze_tensor(image, label):
    # base_color, base_quality_color, map_qual_color, strand_color, alt_color
    img_c, img_w, img_h = image.size()
    img_h = 30
    label_c = label.size(0)
    print()
    for i in range(label_c):
        print(label[i].item(), end='')
    print()
    for i in range(img_h):
        for j in range(img_w):
                print(get_base_from_color(image[0][j][i]), end='')
        print()

    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[1][j][i]), end='')
        print()

    for i in range(img_h):
        for j in range(img_w):
            print(get_quality_by_color(image[2][j][i]), end='')
        print()

    for i in range(img_h):
        for j in range(img_w):
            print(get_strand_color(image[3][j][i]), end='')
        print()

    for i in range(img_h):
        for j in range(img_w):
            print(get_alt_type(image[4][j][i]), end='')
        print()