import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", 'r') as fh:
    requirements = [line.strip() for line in fh.readlines()]


setuptools.setup(
    name="BondGraphTools",
    version="0.4.6",
    author="Pete Cudmore",
    author_email="peter.cudmore@uqconnect.edu.au",
    description="Bond Graph Modelling Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BondGraphTools/BondGraphTools",
    classifiers=(
        "Intended Audience :: Science/Research",
        "Intended Audience :: End Users/Desktop",
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ),
    keywords="modelling control engineering",
    packages=['BondGraphTools'],
    package_dir={'BondGraphTools': 'BondGraphTools'},
    package_data={'BondGraphTools': ['components/*.json']},
    extras_require={
        'docs': [
            'sphinx >= 1.7',
            'sphinx_rtd_theme']},
    install_requires=requirements
#     include_dirs = [r'C:\Program Files\SUNDIALS\include'],
#     library_dirs = [r'C:\Program Files\SUNDIALS\lib'],
#     libraries = ['sundials_cvode', 'sundials_nvecserial', ...]  # 视你使用模块而定
)
