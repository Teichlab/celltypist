#!/usr/bin/env bash
set -euo pipefail


echo -e "[*] Deleting current build from local..."
echo -e "\t- dist/"
rm -rf dist/
echo -e "\t- celltypist_dev.egg-info/"
rm -rf celltypist_dev.egg-info/
echo -e "\t- build/"
rm -rf build/

PACKAGE_NAME=celltypist-dev
echo -e "[*] Uninstall current $PACKAGE_NAME (if installed)..."
pip uninstall $PACKAGE_NAME

echo -e "[*] Copying package data to test folder..."
rm -rf tests/data
cp -r celltypist/data tests/data

echo -e "[*] Run tests..."
python setup.py test

echo -e "[*] Build..."
#python setup.py sdist bdist_wheel

echo -e "[*] Find wheel..."
# WHEEL_FILE="$PWD/$(find dist/*whl)"
# echo -e "[*] $WHEEL_FILE"

echo -e "[*] Install wheel..."
#pip install $WHEEL_FILE
#twine upload dist/*

echo -e "[*] Done!"