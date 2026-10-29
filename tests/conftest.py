"""Every test runs in a throw-away CB_HOME, so nothing lands in the user's own workspace.

`server.common` reads CB_HOME once, when it is first imported. A test that sets it later (test_core's PDF test does)
is too late once another test has imported the module, and the items it saves show up in the app's Recent list.
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ["CB_HOME"] = tempfile.mkdtemp(prefix="cb_test_")
