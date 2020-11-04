#!/usr/bin/env bash
set -euo pipefail


echo -e "[*] Deleting things we don't want to include..."
echo -e "\t- env/"
rm -rf env
echo -e "\t- dist/"
rm -rf dist/
echo -e "\t- celltypist_dev.egg-info/"
rm -rf celltypist_dev.egg-info/
echo -e "\t- build/"
rm -rf build/
echo -e "\t- celltypis/"
echo -e "\t\t- __pycache__"
rm -rf celltypist/__pycache__
echo -e "\t\t- data/models/"
rm -rf celltypist/data/models/*


#PACKAGE_NAME=celltypist-dev
#echo -e "[*] Uninstall current $PACKAGE_NAME (if installed)..."
#pip uninstall $PACKAGE_NAME

echo -e "[*] Run tests..."
echo -e "\t- TODO: make real tests"
#python setup.py test

echo -e "[*] Build..."
python setup.py sdist bdist_wheel

echo -e "[*] Find wheel..."
#WHEEL_FILE="$PWD/$(find dist/*whl)"
#echo -e "[*] $WHEEL_FILE"

echo -e "[*] Install wheel..."
#pip install $WHEEL_FILE
twine upload dist/*

echo -e "[*] Done!"