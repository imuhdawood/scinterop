"""Custom exception hierarchy for scinterop.

All exceptions inherit from :class:`ScinteropError`. Use this base
class to catch any scinterop error generically.

Exception tree::

    ScinteropError
    ├── DetectionError
    ├── FormatAdapterError
    │   ├── MtxAdapterError
    │   ├── H5adAdapterError
    │   └── RdsAdapterError
    ├── ValidationError
    ├── ExecError
    │   ├── RExecError
    │   └── PythonExecError
    ├── CacheError
    └── ProvenanceError
"""

import logging

logger = logging.getLogger("scinterop")


class ScinteropError(Exception):
    """Base exception for all scinterop errors."""


class DetectionError(ScinteropError):
    """Raised when format detection fails on a file or directory."""


class FormatAdapterError(ScinteropError):
    """Base exception for format-specific adapter failures."""


class MtxAdapterError(FormatAdapterError):
    """Raised on MTX read/write failures."""


class H5adAdapterError(FormatAdapterError):
    """Raised on H5AD read/write failures."""


class RdsAdapterError(FormatAdapterError):
    """Raised on RDS read/write failures."""


class ValidationError(ScinteropError):
    """Raised when a CanonicalObject fails structural validation."""


class ExecError(ScinteropError):
    """Base exception for external script execution failures."""


class RExecError(ExecError):
    """Raised when an R script fails or times out."""


class PythonExecError(ExecError):
    """Raised when a Python script fails or times out."""


class CacheError(ScinteropError):
    """Raised on scratch directory or temp file failures."""


class ProvenanceError(ScinteropError):
    """Raised on provenance log read/write failures."""
