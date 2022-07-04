import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fispact_dedx", # Replace with your own username
    version="0.0.1",
    author="Ross Allen",
    author_email="r.allen.4@pgt.bham.ac.uk",
    description="Tool for charged particle activations accounting for energy loss",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Allentro/FISPACT_energy_loss",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib', 
        'mplhep']
        #'https://github.com/Allentro/pysrim.git']
)
