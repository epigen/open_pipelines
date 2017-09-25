""" Utility functions for working with a Sample instance """

from ..const import CHIP_COMPARE_COLUMN


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def determine_comparison(sample, attr=CHIP_COMPARE_COLUMN):
    """
    Return the name of a sample to compare to the one given.

    This relationship is encoded in a Sample's YAML config file, either
    directly or by way of a named column for a particular row (Sample) in a
    sheet of sample annotations, as in the looper project.

    :param Sample sample: the Sample instance for which to determine name of
        comparison Sample
    :param str attr: name of the attribute on the given Sample that should
        store the name of the comparison Sample
    :return str: name of comparison Sample
    :raise AttributeError: if the given Sample lacks the attribute that
        allegedly stores the name of the comparison Sample
    """
    return getattr(sample, attr)
