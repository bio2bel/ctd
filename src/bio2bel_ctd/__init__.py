# -*- coding: utf-8 -*-

from . import enrich, manager
from .enrich import *
from .manager import *

__all__ = [
    manager.__all__ +
    enrich.__all__
]

__version__ = '0.0.1-dev'

__title__ = 'bio2bel_ctd'
__description__ = "A package for converting CTD to BEL"
__url__ = 'https://github.com/bio2bel/ctd'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'charles.hoyt@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Charles Tapley Hoyt'
