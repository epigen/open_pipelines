""" Utility functions for working with a sample annotation sheet """

from functools import partial
import os
import pandas as pd
from ..const import CHIP_COMPARE_COLUMN, CHIP_MARK_COLUMN

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


_DELIMITER_BY_EXTENSION = {".csv": ",", ".txt": "\t", ".tsv": "\t"}



def write_annotations_values(sheet,
        delimiter=None, output_path=None, **update_by_column):
    """
    Use or derive values for column(s) for an annotations sheet.

    The values are either specified directly in the given keyword arguments,
    or they're inferred from existing data in the annotations sheet, using
    functions in the given keyword arguments. Each key in the variable
    keyword arguments should correspond to an existing annotations field to
    update or one that's intended to be added. A value that's a single string
    is applied for all samples while a collection of raw strings is used as a
    sequence of individual values, one per sample, and is applied as such.
    A value that's a function is interpreted as a strategy with which to
    derive values for the column (key) to which it's mapped. Each such
    function should accept the original, parsed annotations table (DataFrame)
    as an argument.

    :param str sheet: path to file with sample annotations data
    :param str delimiter: separator between values and field names in the
        annotations sheet file
    :param str output_path: filepath to which to write the result;
        this could be the same or different than the input file. If different,
        the original file is maintained and the new one is written. If not
        provided, the original file is replaced with the result.
    :param Mapping[str, str | Iterable[str] | function] update_by_column:
        content (or way to generate it) for new or existing column in the
        given table.
    """
    pass



def write_chip_marks(sample_sheet, marks, mark_col_name=CHIP_MARK_COLUMN,
        replace_marks=True, delimiter=None,
        output_path=None, replace_existing_file=False):
    """
    Add ChIP mark(s) to a sample annotations sheet.
    
    Since this is intended to have package-external use, it's assumed that 
    there's some user familiarity with the data at hand and that therefore it's 
    been definitely determined that marks needed to be added updated, so the 
    default behavior is to replace existing marks. The default behavior is
    more conservative with respect to output, though, returning a parsed and
    modified data table rather than writing it to disk or going even further
    and overwriting the existing file. Those results can be achieved with
    provision of a valid argument to one of the appropriate corresponding
    parameters, though.

    :param str | pandas.core.frame.DataFrame sample_sheet: path to annotations
        sheet file or DataFrame as would be parsed from one
    :param function | str | Iterable[str] marks: if a single string, a single
        ChIP mark to apply as the value across all of the sheet's samples. If
        instead a collection of raw text labels, then the length must match
        the number of samples, as the value is interpreted as a sequence of
        mark names, with the order of the sequence matching that of the
        Samples. If a function, it should accept the DataFrame with the
        Sample annotations as input and generate a single label or a
        vector of mark labels to assign for the samples.
    :param bool replace_marks: whether to overwrite existing ChIP mark names
        if they're present
    :param str delimiter: delimiter used to separate values; this is only
        relevant if the argument to sample_sheet is a filepath (i.e., it's
        meaningless if given an already-parsed table of annotations)
    :param output_path: path to output file, optional; if not provided and
        the flag to replace an existing file is not set, a DataFrame is
        returned rather than anything being written to disk.
    :param bool replace_existing_file: whether to overwrite the given
        sample_sheet (relevant if that's a filepath)
    :return NoneType | pandas.core.frame.DataFrame: null if output is written
        to disk, otherwise the DataFrame representing sample annotations
    """

    # Allow filepath as input.
    if isinstance(sample_sheet, str):
        # Infer delimiter if needed.
        if delimiter is None:
            _, ext = os.path.splitext(sample_sheet)
            try:
                delimiter = _DELIMITER_BY_EXTENSION[ext.lower()]
            except KeyError:
                raise ValueError(
                    "No delimiter provided and it cannot be inferred from "
                    "extension '{}' (file '{}')".format(ext, sample_sheet))

        # Parse the annotations file.
        sample_sheet = pd.read_csv(sample_sheet, sep=delimiter)

        if replace_existing_file:
            # Set the output path to the input path to simplify the later
            # conditional on output path for whether to write to disk.
            output_path = sample_sheet

    # Restrict input to filepath or DataFrame.
    elif not isinstance(sample_sheet, pd.DataFrame):
        raise TypeError("Sample sheet must be a filepath or a pandas "
                        "DataFrame: {}".format(type(sample_sheet)))

    # Do the actual mark addition.
    sample_sheet = _add_chip_marks(sample_sheet,
        marks, mark_col_name=mark_col_name, replace=replace_marks)

    # Return the annotations sheet as DataFrame or write to disk.
    if output_path is None:
        return sample_sheet
    sample_sheet.to_csv(output_path, index=False, sep=delimiter or ",")



def _add_chip_marks(sample_sheet, marks,
        mark_col_name=CHIP_MARK_COLUMN, replace=False):
    """
    Add ChIP mark names to sample annotations table.

    To the given table of sample annotations (DataFrame), add ChIP mark(s)
    that are either directly and explicitly specified or that are inferred
    with a given callable, to which the annotations table is passed as an
    argument.

    """

    if mark_col_name in sample_sheet.columns and not replace:
        print("Column '{}' already exists, doing nothing".format(mark_col_name))

    else:
        if isinstance(marks, str):
            # Recycle single raw mark name for each Sample.
            pass
        elif hasattr(marks, "__call__"):
            marks = marks(sample_sheet)
        elif len(sample_sheet) != len(marks):
            raise ValueError("{} rows and {} marks".format(
                len(sample_sheet), len(marks)))
        sample_sheet[mark_col_name] = marks

    return sample_sheet



def _update_table(table, values, column, replace=False):
    if column in table.columns and not replace:
        print("Column '{}' already exists, doing nothing".format(column))
    else:
        if isinstance(values, str):
            # Recycle single raw mark name for each Sample.
            pass
        elif hasattr(values, "__call__"):
            values = values(table)
        elif len(table) != len(values):
            raise ValueError("{} rows and {} values".format(
                len(table), len(values)))
        table[column] = values
    return table


add_chip_marks = partial(_update_table(column=CHIP_MARK_COLUMN))
designate_control = partial(_update_table, column=CHIP_COMPARE_COLUMN)

