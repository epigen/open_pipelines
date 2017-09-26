""" Utility functions for working with a sample annotation sheet """

from functools import partial
import os
import sys
if sys.version_info < (3, 3):
    from collections import Sized
else:
    from collections.abc import Sized
import pandas as pd
from ..const import CHIP_COMPARE_COLUMN, CHIP_MARK_COLUMN


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


_DELIMITER_BY_EXTENSION = {".csv": ",", ".txt": "\t", ".tsv": "\t"}



def write_annotations_values(path_sheet_file,
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
    function should accept the a parsed annotations table, (DataFrame)
    as an argument.

    :param str path_sheet_file: path to file with sample annotations data
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
    
    if not isinstance(path_sheet_file, str):
        raise TypeError("Sample sheet must be path to file with tabular data.")

    # Infer delimiter if needed.
    if delimiter is None:
        _, ext = os.path.splitext(path_sheet_file)
        try:
            delimiter = _DELIMITER_BY_EXTENSION[ext.lower()]
        except KeyError:
            raise ValueError(
                "No delimiter provided and it cannot be inferred from "
                "extension '{}' (file '{}')".format(ext, path_sheet_file))

    # Parse the annotations file.
    table = pd.read_csv(path_sheet_file, sep=delimiter)

    try:
        update_by_column = update_by_column.items()
    except AttributeError:
        # Interpret the input as already formatted as a collection of pairs.
        pass

    # For each column, finalize the values to update.
    update_by_column = map(
        lambda vals_col_pair: _format_values_for_table(
            values=vals_col_pair[0], table=table), update_by_column)
    for column, values in update_by_column:
        table[column] = values
        
    output_path = output_path or path_sheet_file
    table.to_csv(output_path, index=False, sep=delimiter)



def _format_values_for_table(values, table):
    """
    Do some validation and finalization of values for table insertion.

    :param Iterable | function values: A single value for all rows of a table,
        a collection of individual values, or a function with which to infer
        values when given the table as an argument; if a collection, the size
        must match the number of rows in the table; if a function, it must
        take the table as an argument.
    :param pandas.core.frame.DataFrame table: the sample annotations table
        for which the values are destined
    :return object | Iterable[object]: either a single value to use for each row
        in the table, or a collection of such values
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError("Table must be a DataFrame; got {}".format(type(table)))
    if isinstance(values, Sized):
        if len(values) != len(table):
            raise ValueError(
                "{} values for {} rows".format(len(values), len(table)))
    elif hasattr(values, "__call__"):
        values = values(table)
    return values



def _update_table(table, values, column, replace=False):
    """
    Update a data table (e.g., annotations sheet) with values for a column.

    Single value is supported (set for each row), as is a collection of
    atomic values (must be exactly one for each row). A function with which to
    derive values when passed the table itself as an argument is also valid
    as an argument for the values.

    :param pandas.core.frame.DataFrame table:
    :param object | Iterable[object] | function values: value for each row,
        or a way to derive such values when given the table
    :param str column: name of the column to set or update
    :param bool replace: whether to replace the column indicated if it
        already exists
    :return pandas.core.frame.DataFrame: the updated table
    """
    if column in table.columns and not replace:
        print("Column '{}' already exists, doing nothing".format(column))
    else:
        values = _format_values_for_table(values=values, table=table)
        table[column] = values
    return table



# Specific versions of more generic function (i.e., partially-parameterized)
add_chip_marks = partial(_update_table, column=CHIP_MARK_COLUMN)
designate_control = partial(_update_table, column=CHIP_COMPARE_COLUMN)
