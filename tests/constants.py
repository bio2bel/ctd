# -*- coding: utf-8 -*-

import logging
import os

from bio2bel.testing import make_temporary_cache_class_mixin
from bio2bel_ctd import Manager

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))

TemporaryDatabaseMixin = make_temporary_cache_class_mixin(Manager)


class PopulatedDatabaseMixin(TemporaryDatabaseMixin):
    @classmethod
    def setUpClass(cls):
        super(PopulatedDatabaseMixin, cls).setUpClass()
        cls.manager.populate() # FIXME use test data
