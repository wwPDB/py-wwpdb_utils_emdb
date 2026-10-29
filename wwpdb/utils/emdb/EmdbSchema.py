##
# File:    EmdbSchema.py
# Date:    29-March-2021
"""
Wrapper for access to the EMDB schema distributed by :mod:`emdb_schemas`.
"""

from contextlib import contextmanager

try:
    from importlib.resources import as_file
    from importlib.resources import files
except ImportError:  # Python 3.8
    from importlib_resources import as_file
    from importlib_resources import files


class EmdbSchema(object):
    @staticmethod
    def _schema_resource():
        return files("emdb_schemas").joinpath("emdb.xsd")

    def getSchemaPath(self):
        """Return the schema path for a conventional, unpacked installation.

        New callers that require a filesystem path should use
        :meth:`getSchemaPathContext`, which also supports zip-based imports.
        """
        return str(self._schema_resource())

    @contextmanager
    def getSchemaPathContext(self):
        """Yield a filesystem path to the packaged EMDB XSD.

        ``importlib.resources.as_file`` keeps a temporary file available for
        the duration of the context when the dependency is installed from a
        zip-based distribution.
        """
        with as_file(self._schema_resource()) as schema_path:
            yield str(schema_path)
