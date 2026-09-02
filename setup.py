from pathlib import Path

from setuptools import find_packages, setup


PROJECT_ROOT = Path(__file__).parent
README_PATH = PROJECT_ROOT / "README.md"

DESCRIPTION = (
    "CellBLASTer: a universal plant scRNA-seq annotation tool inspired by "
    "cellular BLAST strategies"
)


setup(
    name="CellBLASTer",
    version="1.0.0",
    author="Lin Du",
    author_email="3051065449@qq.com",
    description=DESCRIPTION,
    long_description=(
        README_PATH.read_text(encoding="utf-8")
        if README_PATH.is_file()
        else DESCRIPTION
    ),
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(include=["CellBLASTer", "CellBLASTer.*"]),
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.23",
        "pandas>=1.5",
        "requests>=2.28",
        "matplotlib>=3.6",
        "seaborn>=0.12",
        "scipy>=1.9",
        "scikit-learn>=1.2",
        "numba>=0.57",
        "scanpy>=1.9",
    ],
    entry_points={
        "console_scripts": [
            "cellblaster=CellBLASTer.CellBlaster:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        "single-cell RNA-seq",
        "cell-type annotation",
        "plant",
        "orthogroup",
        "bioinformatics",
    ],
    zip_safe=False,
)
