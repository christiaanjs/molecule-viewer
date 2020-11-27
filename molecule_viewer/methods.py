from enum import Enum

class ConformationMethod(Enum):
    """An enum for 3D conformer generation methods.

    MMMFF94: Halgren, T. A. (1996). Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94.
    UFF: Rappé, A. K., Casewit, C. J., Colwell, K. S., Goddard III, W. A., & Skiff, W. M. (1992). UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations.
    ETKDG: Riniker, S., & Landrum, G. A. (2015). Better informed distance geometry: using what we know to improve conformation generation.

    """
    MMFF94 = "MMFF94"
    UFF = "UFF"
    ETKDG = "ETKDG"

_default_conformation_method = ConformationMethod.MMFF94

def _set_default_conformation_method(method):
    global _default_conformation_method
    _default_conformation_method = method