""" Exceptional pipeline conditions """

__author__ = "Vince Reuter"
__email__ ="vreuter@virginia.edu"


__all__ = ["InvalidFiletypeException"]



class InvalidFiletypeException(Exception):
    """ A pipeline will likely have a restricted domain of input types. """
    pass
