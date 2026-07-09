import logging

logger = logging.getLogger("scinterop")


class ScinteropError(Exception):
    pass


class DetectionError(ScinteropError):
    pass


class FormatAdapterError(ScinteropError):
    pass


class MtxAdapterError(FormatAdapterError):
    pass


class H5adAdapterError(FormatAdapterError):
    pass


class RdsAdapterError(FormatAdapterError):
    pass


class ValidationError(ScinteropError):
    pass


class ExecError(ScinteropError):
    pass


class RExecError(ExecError):
    pass


class PythonExecError(ExecError):
    pass


class CacheError(ScinteropError):
    pass


class ProvenanceError(ScinteropError):
    pass
