from enum import Enum

class ConformationMethod(Enum):
    MMFF94 = "MMFF94"
    UFF = "UFF"
    ETKDG = "ETKDG"

_default_conformation_method = ConformationMethod.MMFF94

def _set_default_conformation_method(method):
    global _default_conformation_method
    _default_conformation_method = method