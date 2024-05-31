import warnings
from contextlib import contextmanager

@contextmanager
def suppress_warnings():
    warnings.filterwarnings("ignore")
    try:
        yield
    finally:
        warnings.resetwarnings()