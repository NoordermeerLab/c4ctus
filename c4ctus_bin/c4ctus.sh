#!/bin/sh

A=`dirname $0`
value=`cat "${A}/help.txt"`
echo "$value" | fold -w 80 -s
