"""Public Python interface for CellBLASTer."""

from .CellBlaster import (
    DATABASE_CONFIG,
    DEFAULT_NONCODING_RNA_KEYWORDS,
    CellBlaster,
)

__version__ = "1.0.0"

__all__ = [
    "CellBlaster",
    "DATABASE_CONFIG",
    "DEFAULT_NONCODING_RNA_KEYWORDS",
    "__version__",
]
