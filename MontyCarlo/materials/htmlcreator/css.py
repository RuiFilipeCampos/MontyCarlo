from copy import deepcopy

_DEFAULT_CSS = {
    'body': {
        'font-family': 'sans-serif',
        'max-width': '80%',
        'margin': 'auto',
        'margin-bottom': '80px',
    },
    'div.pandas-dataframe': {
        'overflow': 'auto',
    },
    'table.dataframe': {
        'border-collapse': 'collapse',
        'border': 'none',
    },
    'table.dataframe td': {
        'background-color': '#fffef4',
        'font-size': '14px',
        'text-align': 'center',
        'white-space': 'nowrap',
        'margin': '0',
        'padding-top': '0.4em',
        'padding-bottom': '0.4em',
        'padding-left': '0.5em',
        'padding-right': '0.5em',
        'border': '1px solid #d7d7d7',
    },
    'table.dataframe th': {
        'background-color': '#f9f9f9',
        'font-size': '12px',
        'text-align': 'center',
        'white-space': 'nowrap',
        'padding-left': '1em',
        'padding-right': '1em',
        'padding-top': '0.5em',
        'padding-bottom': '0.5em',
        'border': '1px solid #e0e0e0',
    },
    'table.dataframe tr:nth-child(even) td': {
        'background': '#fffaf0',
    },
    'table.dataframe tr:nth-child(odd) td': {
        'background': '#fffff0',
    },
    'table.dataframe tr:hover td': {
        'background-color': '#ddeeff',
    },
    'table.dataframe tbody th': {
        'font-weight': 'normal',
    },
}


class CSS(dict):
    """CSS container."""

    def __init__(self):
        self.update(deepcopy(_DEFAULT_CSS))

    def __str__(self) -> str:
        lines = []
        for selector, declaration in self.items():
            lines.append(f'{selector} {{ ')
            for property, value in declaration.items():
                lines.append(f'{property}: {value}; ')
            lines.append('} ')
        lines[-1] = lines[-1].rstrip()
        return ''.join(lines)
