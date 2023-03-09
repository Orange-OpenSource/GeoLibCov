from setuptools import setup, find_packages

setup(
    name='DataCov',
    version="0.0.1",
    description="A small package for cell topography generation and coverage modelling",
    author='Danny Qiu',
    author_email='danny.qiu@orange.com',
    license='To be defined.',
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
    ]
)