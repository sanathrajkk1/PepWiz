"""PepWiz — GUI-based peptide MS/MS analysis and visualization tool."""

# Version
__version__ = "0.1.0"

# Public API (lightweight import surface)
from .match_engine import (
    calc_fragments,
    generate_theoretical_by,
    legacy_summary_from_spectrum,
    nearest_match,
    ion_meta,
    ppm_error,
    PROTON,
    WATER,
)

from .mzml_utils import (
    open_reader,
    precursor_mz_from_spec,
    list_precursors_with_counts,
    average_spectrum,
)

from .msconvert_utils import (
    find_msconvert,
    run_msconvert,
)

# Optional: visualization & legacy writer if you want them importable too
try:
    from .visualize import export_fragment_image, export_annotated_spectrum
    from .io_legacy import write_legacy_out
except Exception:
    # These modules might be moved/renamed later; don't break imports if missing.
    pass

__all__ = [
    "calc_fragments", "generate_theoretical_by", "legacy_summary_from_spectrum",
    "nearest_match", "ion_meta", "ppm_error", "PROTON", "WATER",
    "open_reader", "precursor_mz_from_spec", "list_precursors_with_counts", "average_spectrum",
    "find_msconvert", "run_msconvert",
    "export_fragment_image", "export_annotated_spectrum", "write_legacy_out",
    "__version__",
]
