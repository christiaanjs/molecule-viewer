import pathlib

def _set_wd(wd):
    global _wd
    _wd = pathlib.Path(wd)
    _wd.mkdir(parents=True, exist_ok=True)

_set_wd("molecules")
