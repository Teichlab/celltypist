#pip install sphinx sphinx_rtd_theme
cd docs/
sphinx-apidoc --force --no-toc --module-first -o source ../celltypist
make clean
make html
