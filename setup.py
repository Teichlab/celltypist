import setuptools

with open("README.md", "rt", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "rt", encoding="utf-8") as fh:
    install_requirements = [line.strip() for line in fh.readlines()]

setuptools.setup(
    name="celltypist-dev",
    version="0.1.10",
    author="Tomas Pires de Carvalho Gomes",
    author_email="tpcg@sanger.ac.uk",
    description="A tool for semi-automatic cell type annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cellgeni/celltypist",
    packages=setuptools.find_packages(),
    install_requires=install_requirements,
    include_package_data=True,
    entry_points={
        'console_scripts': ['celltypist=celltypist.command_line:main'],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    python_requires='>=3.6',
)
