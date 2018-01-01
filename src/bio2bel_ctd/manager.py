# -*- coding: utf-8 -*-

__all__ = [
    'Manager'
]


class Manager(object):
    def __init__(self, connection=None):
        pass

    def populate(self):
        """Populates the database"""

    def enrich_chemicals(self, graph):
        """Finds chemicals that can be mapped and enriched with the CTD

        :param pybel.BELGraph graph:
        """
        raise NotImplementedError

    @staticmethod
    def ensure(connection=None):
        """Checks and allows for a Manager to be passed to the function.

        :param connection: can be either a already build manager or a connection string to build a manager with.
        """
        if connection is None or isinstance(connection, str):
            return Manager(connection=connection)

        if isinstance(connection, Manager):
            return connection

        raise TypeError
