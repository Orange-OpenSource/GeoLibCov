# Software Name: DataCov
# Version: 0.1.0
# SPDX-FileCopyrightText: Copyright (c) 2023 Orange
# SPDX-License-Identifier: BSD-3-Clause
#
# This software is distributed under the BSD 3-Clause "New" or "Revised" License,
# see the "LICENSE.txt" file for more details.
#
# Author: Danny Qiu <danny.qiu@orange.com>

from setuptools import setup, find_packages

setup(
    name='DataCov',
    version="0.1.0",
    description="A small package for cell topography generation and coverage modelling",
    author='Danny Qiu',
    author_email='danny.qiu@orange.com',
    license='BSD 3-Clause',
    classifiers=[  
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Researchers",
        "Programming Language :: Python :: 3 :: Only",
    ],
    packages=[
        'datacov',
    ],
    install_requires=[
        'numpy', 
        'geopandas',
        'pandas',
        'matplotlib',
        'folium',
        'mapclassify',
        'tqdm',
        'ipykernel',
        'shapely',
        'scipy',
    ]
)