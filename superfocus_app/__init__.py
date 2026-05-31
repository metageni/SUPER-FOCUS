# -*- coding: utf-8 -*-

try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
    try:
        version = _pkg_version("superfocus")
    except PackageNotFoundError:
        version = "1.8"
except ImportError:
    version = "1.8"
