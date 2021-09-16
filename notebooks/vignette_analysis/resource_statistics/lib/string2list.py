"""
Tools for handling info in mitochondria tables
"""


def string2list(stringlist,targtype=int):
    """
    In DataFrames that have been saved to .csv files,
    lists are sometimes converted to string types.
    In mitochondrial dataframes this is often the case.
    This tool converts a list with the contents of a
    DataFrame column (with string-type elements) to a
    list with list-type elements. The content of each
    list element is the defined "target type".
    @stringlist: list of string-type list elements
    @targtype: numeric type of each list element's
    items. Default is int.
    """
    temp = ''
    if stringlist.find(',') != -1:
        temp = stringlist.strip('][').split(',')
    else:
        temp = stringlist.strip('][').split()
    return [targtype(q) for q in temp if q != '']





