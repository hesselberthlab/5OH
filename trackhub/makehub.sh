#! /usr/bin/env bash

rst="description.rst"
html="description.html"
rst2html.py $rst $html

url="http://amc-sandbox.ucdenver.edu/~jhessel/hubs/5OH/hub.txt"

cp -r * ~/public_html/hubs/5OH/

hubCheck $url

echo "hub is at $url"
