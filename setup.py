import setuptools

def get_readme():
    with open("README.md", "rt", encoding="utf-8") as fh:
        return fh.read()

def get_requirements():
    with open("requirements.txt", "rt", encoding="utf-8") as fh:
        return [line.strip() for line in fh.readlines()]

def get_version():
    with open("celltypist/__init__.py", "rt", encoding="utf-8") as fh:
        for line in fh.readlines():
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1].strip()
    raise RuntimeError("Unable to find version string in celltypist/__init__.py")

setuptools.setup(
    name="celltypist",
    version=get_version(),
    author="Chuan Xu",
    author_email="cx1@sanger.ac.uk",
    description="A tool for semi-automatic cell type classification",
    long_description=get_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Teichlab/celltypist",
    packages=setuptools.find_packages(),
    install_requires=get_requirements(),
    include_package_data=True,
    entry_points={
        'console_scripts': ['celltypist=celltypist.command_line:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 5 - Production/Stable",
    ],
    python_requires='>=3.6',
)
