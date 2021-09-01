#!/usr/bin/env python

"""
===========================
Get Restriction site from Restriction enzyme
===========================
"""


from Bio.Restriction.Restriction import RestrictionBatch
import sys


RestrEnzyme = sys.argv[1]

batch = RestrictionBatch()
batch.add(RestrEnzyme)
enzyme = batch.get(RestrEnzyme)
print enzyme.site
