"""fitstability: peak-fit stability diagnostics for XPS.

Implements the multi-start box-plot test of Moeini et al., J. Electron Spectrosc. Relat.
Phenom. 250 (2021) 147094, with the Monte Carlo and uniqueness-plot comparators it discusses.
"""
from .experiments import (MonteCarloResult, MultiStartResult, UniquenessResult, monte_carlo,
                          multistart, random_start, scannable_parameters, uniqueness,
                          unscannable_reason, UNIQUENESS_METHOD_NOTE)
from .casa import CasaFitError, CasaModel, casa_models, parse_fit
from .library import OriginFit, ReferenceLibrary, ReferenceSpectrum, calibrate_origins
from .fit import fit, fit_nominal
from .io import (ProtocolSet, load_protocol, load_protocols, protocol_to_cfg, read_spectrum, save_protocol,
                 sniff_table, synthetic_spectrum)
from .model import FitProblem, FitResult, PeakSpec, Protocol, StartPoint
from .vamas import Block, VamasError, VamasFile, read_vamas, write_vamas
from .verdict import Verdict, assess, summary_table

__all__ = [
    "FitProblem", "FitResult", "PeakSpec", "Protocol", "StartPoint", "fit", "fit_nominal",
    "multistart", "monte_carlo", "uniqueness", "random_start",
    "MultiStartResult", "MonteCarloResult", "UniquenessResult",
    "Verdict", "assess", "summary_table",
    "read_spectrum", "sniff_table", "save_protocol", "protocol_to_cfg", "load_protocol", "load_protocols", "ProtocolSet", "read_vamas", "write_vamas",
    "CasaFitError", "CasaModel", "casa_models", "parse_fit",
    "ReferenceLibrary", "ReferenceSpectrum", "OriginFit", "calibrate_origins",
    "VamasFile", "Block", "VamasError", "scannable_parameters", "unscannable_reason", "UNIQUENESS_METHOD_NOTE", "synthetic_spectrum",
]
__version__ = "0.1.0"
