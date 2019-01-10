# -*- coding: utf-8 -*-

import pkg_resources

try:
    version = pkg_resources.require("superfocus")[0].version
except:
    raise ValueError('Cannot find version number')
