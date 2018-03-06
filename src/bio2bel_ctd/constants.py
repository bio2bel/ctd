# -*- coding: utf-8 -*-

from bio2bel.utils import get_connection, get_data_dir

MODULE_NAME = 'ctd'

DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)
