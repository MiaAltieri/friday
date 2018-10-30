class ActiveRegionOptions(object):
    MIN_REGION_SIZE = 80
    MAX_REGION_SIZE = 1000
    REGION_EXPANSION = 25
    MIN_MAPPING_QUALITY = 20
    MIN_BASE_QUALITY = 20
    # the linear regression model is used inside C++ code


class CandidateFinderOptions(object):
    # base and map quality
    MIN_BASE_QUALITY = 14
    MIN_MAP_QUALITY = 15
    SAFE_BASES = 50


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 20

    MAX_READS_IN_REGION = 1500
    RANDOM_SEED = 2719747673


class DeBruijnGraphOptions(object):
    MIN_K = 10
    MAX_K = 100
    STEP_K = 1
    MIN_EDGE_SUPPORT = 2
    MAX_ALLOWED_PATHS = 256

    # base and map quality
    MIN_BASE_QUALITY = 15
    MIN_MAP_QUALITY = 20
