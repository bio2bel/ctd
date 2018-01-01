# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import logging
import os
import tempfile
import unittest

from bio2bel_ctd import Manager

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))


class TemporaryConnectionMixin(unittest.TestCase):
    """Creates a :class:`unittest.TestCase` that has a persistent file for use with SQLite during testing."""

    @classmethod
    def setUpClass(cls):
        """Creates a temporary file to use as a persistent database throughout tests in this class. Subclasses of
        :class:`TemporaryCacheClsMixin` can extend :func:`TemporaryCacheClsMixin.setUpClass` to populate the database
        """
        cls.fd, cls.path = tempfile.mkstemp()
        cls.connection = 'sqlite:///' + cls.path
        log.info('test database at %s', cls.connection)

    @classmethod
    def tearDownClass(cls):
        """Closes the connection to the database and removes the files created for it"""
        os.close(cls.fd)
        os.remove(cls.path)


class TemporaryDatabaseMixin(TemporaryConnectionMixin):
    @classmethod
    def setUpClass(cls):
        super(TemporaryDatabaseMixin, cls).setUpClass()
        cls.manager = Manager(connection=cls.connection)

    @classmethod
    def tearDownClass(cls):
        # add deletion of database
        super(TemporaryDatabaseMixin, cls).tearDownClass()


class PopulatedDatabaseMixin(TemporaryDatabaseMixin):
    @classmethod
    def setUpClass(cls):
        super(PopulatedDatabaseMixin, cls).setUpClass()
        cls.manager.populate()
