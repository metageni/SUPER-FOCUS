# Python Coding Standards

## Shebang & docstring
- First line: `#!/usr/bin/env python3` (NOT `/usr/bin/python3` — breaks conda envs)
- Second: module-level docstring describing the file

## Import order (strictly by line length within each group)

```python
import os          # group 1: native import, shortest first
import csv
import random
import logging

from pathlib import Path          # group 5: native from...import
from collections import defaultdict

import numpy as np                # group 3: third-party import

from scipy.optimize import nnls   # group 7: third-party from...import

from superfocus_app import version          # group 11: local from...import
from superfocus_app.pipeline import (       # multi-symbol: parenthesised, one per line
    run,
)
```

Groups: 1=native import, 3=third-party import, 5=native from...import,
7=third-party from...import, 9=local import, 11=local from...import.
Blank line between every group.

## Docstrings (PEP 257)
All functions, classes, methods must have docstrings with Args/Returns/Raises.

## Logging
Use `logging` module. Never `print()` in core logic.
`logging.basicConfig()` only in `cli.py` (the entry point).

## CLI separation
- `pipeline.py` = core logic only, no CLI code
- `cli.py` = CLI only (click), no business logic

## Trailing newline
Every Python file ends with a blank line.
