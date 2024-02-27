import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Cyclic_peptides",
    version="0.0.1",
    author="Ruben Muñoz",
    author_email="ruben.21.12.m@gmail.com",
    description="A Python package to create, prepare, and simulate cyclic peptides",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
    'matplotlib',
    'rosetta'
    ],
)
