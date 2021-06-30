import base64
import io
import pathlib
from typing import Optional, Union

import bs4
import numpy as np
import pandas as pd
import PIL
import plotly
from bs4 import BeautifulSoup
from PIL import Image

from .css import CSS

_DEFAULT_HTML = (
    '<!DOCTYPE html>'
    '<html lang="en">'
    '<head><meta charset="UTF-8"></head>'
    '<body></body>'
    '</html>'
)

HEADER_LEVEL_VALUES = {'h1', 'h2', 'h3', 'h4', 'h5', 'h6'}
TEXT_ALIGN_VALUES = {'left', 'center', 'right', 'justify', 'inherit', 'start', 'end'}


class HTMLDocument:
    """HTML Document class."""

    def __init__(self) -> None:
        self.soup = BeautifulSoup(_DEFAULT_HTML, 'html.parser')
        self.css = CSS()

    def add_header(
        self,
        header: str,
        level: str = 'h2',
        align: str = 'left',
    ) -> None:
        """Add header."""
        assert level in HEADER_LEVEL_VALUES, level
        assert align in TEXT_ALIGN_VALUES, align
        style = f'text-align: {align};'
        header_tag = self.soup.new_tag(level, style=style)
        header_tag.string = header
        self.soup.body.append(header_tag)

    def add_paragraph(
        self,
        text: str,
        size: str = '16px',
        indent: str = '0',
        align: str = 'left',
    ) -> None:
        """Add text paragraph."""
        assert align in TEXT_ALIGN_VALUES, align
        style = f'font-size:{size}; text-indent: {indent}; text-align: {align};'
        p_tag = self.soup.new_tag('p', style=style)
        p_tag.string = text
        self.soup.body.append(p_tag)

    def add_line_break(self) -> None:
        """Add line break."""
        br_tag = self.soup.new_tag('br')
        self.soup.body.append(br_tag)

    def add_page_break(self) -> None:
        """Add page break."""
        style = 'page-break-after: always;'
        pb_tag = self.soup.new_tag('p', style=style)
        self.soup.body.append(pb_tag)

    def add_image(
        self,
        image: Union[np.ndarray, PIL.Image.Image, pathlib.Path],
        title: Optional[str] = None,
        height: Optional[Union[int, str]] = None,
        width: Optional[Union[int, str]] = None,
        pixelated: bool = False,
    ) -> None:
        """Embed image."""
        image_encoded_str = self._encode_image(image)
        image_src = f'data:image/png;base64, {image_encoded_str}'
        self._add_image_tag(
            src=image_src,
            title=title,
            height=height,
            width=width,
            pixelated=pixelated,
        )

    def add_image_link(
        self,
        image_link: Union[pathlib.Path, str],
        title: Optional[str] = None,
        height: Optional[Union[int, str]] = None,
        width: Optional[Union[int, str]] = None,
        pixelated: bool = False,
    ) -> None:
        """Add image link (Path or URL)."""
        if isinstance(image_link, pathlib.Path):
            image_src = str(image_link)
        elif isinstance(image_link, str):
            image_src = image_link
        else:
            raise TypeError(
                f'image_link is of type {type(image_link)}, '
                f'but it should be {pathlib.Path} or {str}.'
            )
        self._add_image_tag(
            src=image_src,
            title=title,
            height=height,
            width=width,
            pixelated=pixelated,
        )

    def add_table(self, df: pd.DataFrame) -> None:
        """Embed pandas DataFrame."""
        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f'df is of type {type(df)}, but it should be of type {pd.DataFrame}.'
            )
        div_tag = self.soup.new_tag('div', attrs={'class': 'pandas-dataframe'})
        table = BeautifulSoup(str(df.to_html()), 'html.parser')
        self._simplify_double_thead_tr(table)
        div_tag.append(table)
        self.soup.body.append(div_tag)

    def add_plotly_figure(self, fig: plotly.graph_objs.Figure) -> None:
        """Add plotly figure."""
        if not isinstance(fig, plotly.graph_objs.Figure):
            raise TypeError(
                f'fig is of type {type(fig)}, '
                f'but it should be {plotly.graph_objs.Figure}.'
            )
        plotly_figure_html = plotly.io.to_html(
            fig=fig,
            full_html=False,
            include_plotlyjs='cdn',
        )
        div_tag = self.soup.new_tag('div', attrs={'class': 'plotly-figure'})
        plotly_figure = BeautifulSoup(str(plotly_figure_html), 'html.parser')
        div_tag.append(plotly_figure)
        self.soup.body.append(div_tag)

    def set_title(self, title: str) -> None:
        """Set document title."""
        for child in list(self.soup.head.children):
            if child.name == 'title':
                child.decompose()
        title_tag = self.soup.new_tag('title')
        title_tag.string = str(title)
        self.soup.head.append(title_tag)

    def write(self, filepath: str) -> None:
        """Save document to filepath."""
        self._set_style()
        with io.open(str(filepath), 'w', encoding='utf8') as f:
            f.write(str(self.soup.prettify()))

    def _add_image_tag(
        self,
        src: str,
        title: Optional[str] = None,
        height: Optional[Union[int, str]] = None,
        width: Optional[Union[int, str]] = None,
        pixelated: bool = False,
    ) -> None:
        """Add image tag."""
        style = 'border:1px solid #021a40; margin: 3px 3px;'
        if pixelated:
            style += ' image-rendering: pixelated;'
        img_tag = self.soup.new_tag('img', src=src, style=style)
        if title:
            img_tag.attrs['title'] = str(title)
        if height:
            img_tag.attrs['height'] = str(height)
        if width:
            img_tag.attrs['width'] = str(width)
        self.soup.body.append(img_tag)

    def _encode_image(
        self,
        image: Union[np.ndarray, PIL.Image.Image, pathlib.Path],
    ) -> str:
        """Encode image to base64 string."""
        if isinstance(image, np.ndarray):
            if image.dtype != np.uint8:
                raise RuntimeError(
                    f'image.dtype is {image.dtype}, but it should be uint8.'
                )
            if not (image.ndim == 2 or image.ndim == 3):
                raise RuntimeError(
                    f'image.ndim is {image.ndim}, but it should be 2 or 3.'
                )
            buff = io.BytesIO()
            Image.fromarray(image).save(buff, format='PNG')
            encoded = base64.b64encode(buff.getvalue())
        elif isinstance(image, PIL.Image.Image):
            buff = io.BytesIO()
            image.save(buff, format='PNG')
            encoded = base64.b64encode(buff.getvalue())
        elif isinstance(image, pathlib.Path):
            encoded = base64.b64encode(open(str(image), 'rb').read())
        else:
            raise TypeError(
                f'image is of type {type(image)}, but it should be one of: '
                f'{np.ndarray}, {PIL.Image.Image} or {pathlib.Path}.'
            )
        image_encoded_str = encoded.decode('utf-8')
        return image_encoded_str

    def _set_style(self) -> None:
        """Set document style."""
        for child in list(self.soup.head.children):
            if child.name == 'style':
                child.decompose()
        style_tag = self.soup.new_tag('style')
        style_tag.string = str(self.css)
        self.soup.head.append(style_tag)

    def _simplify_double_thead_tr(self, table: bs4.BeautifulSoup) -> None:
        """Simplify a double header row in HTML table generated by pandas."""
        thead = table.find('thead')
        tr_list = thead.find_all('tr')
        if len(tr_list) != 2:
            return
        tr1, tr2 = tr_list
        th_list1 = tr1.find_all('th')
        th_list2 = tr2.find_all('th')
        if len(th_list1) != len(th_list2):
            return
        for th1, th2 in zip(th_list1, th_list2):
            if th1.string and th2.string:
                return
        for th1, th2 in zip(th_list1, th_list2):
            if not th1.string:
                th1.string = th2.string
        tr2.decompose()
