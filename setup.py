from setuptools import setup, find_packages

setup(
    name='PhyloProPy',
    version='0.1',
    packages=find_packages(),
    description='Read, write and modify PhyloProfile files in Python',
    author='Felix Langschied',
    author_email='langschied@bio.uni-frankfurt.de',
    install_requires=[
        'ete3',
        'pandas',
        'scikit-learn',
        'plotly',
        'numpy',
        'seaborn',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': ["phyloSNE = PhyloProPy.standalone_tsne:main"],
    },
)
