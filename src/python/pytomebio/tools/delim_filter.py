import os
import attr
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
from typing import Union


@attr.s(frozen=True)
class ColumnFilter:
    """A filter for a column in a delimited data file.

    Attributes:
        key: the column name to filter
        value: the value (integer, float, or string)
    """

    key: str = attr.ib()
    value: Union[str, int, float] = attr.ib()

    def fail_min(self, datum: Dict[str, Union[str, int, float]]) -> bool:
        """Returns true the given data falls below this filter's value."""
        assert isinstance(datum[self.key], type(self.value))
        return datum[self.key] < self.value  # type: ignore

    def fail_max(self, datum: Dict[str, Union[str, int, float]]) -> bool:
        """Returns true the given data falls above this filter's value."""
        assert isinstance(datum[self.key], type(self.value))
        return datum[self.key] > self.value  # type: ignore

    def __str__(self) -> str:
        return f"{self.key}:{self.value}"

    @classmethod
    def build(cls, keyval: str) -> "ColumnFilter":
        """Builds a ColumnFilter from a key-value tuple that's colon-delimited."""
        (key, value) = keyval.split(":", maxsplit=1)
        return ColumnFilter(key=key, value=convert_type(value))


def convert_type(value: str) -> Union[str, int, float]:
    """Converts a string value to a numeric type, if possible.

    If there's a period in the value, attempts to convert to a float,
    otherwise to an integer.  If the numeric conversion is unsuccessful,
    return the original string.

    Args:
        value: the value to convert
    """
    try:
        return float(value) if "." in value else int(value)
    except ValueError:
        return value


def delim_filter(
    *,
    in_txt: Path,
    sort_by: Optional[List[str]] = None,
    min_filter: Optional[List[ColumnFilter]] = None,
    max_filter: Optional[List[ColumnFilter]] = None,
    filter_any: bool = False,
    reverse_sort: bool = False,
    prepend_file_name: Optional[str] = None,
    delimiter: str = "\t",
    header: bool = True,
) -> None:
    """Filter a delimited data file.

    Does not support quotes in the columns.

    Args:
        in_txt: the path to the tab-delimited file to filter.
        sort_by: sort by the given column.  May be specified multiple times.
        min_filter: filters the given column name to have at least the minimum value.  Should
                    be specified as the column name and minimum value, colon delimited (e.g.
                    `<column>:<value>`).  May be specified multiple times.
        max_filter: filters the given column name to have at least the maximum value.  Should
                    be specified as the column name and maximum value, colon delimited (e.g.
                    `<column>:<value>`).  May be specified multiple times.
        filter_any: true to filter a row if any of the filters fail, false to filter a row if
                    all of the filters fail.
        reverse_sort: sort descending (defaults to ascending)
        prepend_file_name: prepend the file name in the first column of the output.  This will
                           trim the file name up to the first period.  The value should be the
                           name of the column to use.
        delimiter: the input delimiter
        header: include the header in the output.
    """
    first_col_value: str = in_txt.name.split(".", maxsplit=1)[0]
    with in_txt.open("r") as reader:
        in_header = reader.readline().rstrip(os.linesep).split(delimiter)
        assert sort_by is None or all(
            col in in_header for col in sort_by
        ), "--sort-by not in the header: " + ", ".join(
            col for col in sort_by if col not in in_header
        )
        assert min_filter is None or all(
            f.key in in_header for f in min_filter
        ), "--min-filter not in the header: " + ", ".join(
            f.key for f in min_filter if f.key not in in_header
        )
        assert max_filter is None or all(
            f.key in in_header for f in max_filter
        ), "--max-filter not in the header: " + ", ".join(
            f.key for f in max_filter if f.key not in in_header
        )
        data = []
        for line_number, line in enumerate(reader, 2):
            fields = [convert_type(value) for value in line.rstrip(os.linesep).split(delimiter)]
            assert len(in_header) == len(
                fields
            ), f"# of fields mismatched the header on line {line_number}"
            datum: Dict[str, Union[str, int, float]] = {}
            if prepend_file_name is not None:
                datum = {prepend_file_name: first_col_value}
            datum.update(dict(zip(in_header, fields)))

            if filter_any:
                if min_filter is not None and any(f.fail_min(datum) for f in min_filter):
                    continue
                if max_filter is not None and any(f.fail_max(datum) for f in max_filter):
                    continue
            else:
                if min_filter is not None and all(f.fail_min(datum) for f in min_filter):
                    continue
                if max_filter is not None and all(f.fail_max(datum) for f in max_filter):
                    continue

            data.append(datum)

        if sort_by is not None:
            data = sorted(data, key=lambda datum: tuple(datum[col] for col in sort_by))
            if reverse_sort:
                data = data[::-1]

        if header:
            if prepend_file_name is not None:
                in_header = [prepend_file_name] + in_header
            print(delimiter.join(in_header))
        for datum in data:
            print(delimiter.join(str(value) for value in datum.values()))
