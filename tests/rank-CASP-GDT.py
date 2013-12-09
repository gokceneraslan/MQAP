#!/usr/bin/env python

## MQAP
# Copyright (C) 2013 Gokcen Eraslan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import urllib2
import csv

import pandas as pd

RESULT_URI = 'http://www.predictioncenter.org/download_area/CASP10/results_LGA_sda/%s.SUMMARY.lga_sda.txt'

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--target', help='Name of the target',
                              required=True)
parser.add_argument('-o', '--output', help='Output file in CSV format',
                              required=True)

args = parser.parse_args()
resultcontent = urllib2.urlopen(RESULT_URI % args.target)
#resultcontent = open('T0687.SUMMARY.lga_sda.txt')

resultdf = pd.read_csv(resultcontent, delim_whitespace=True)
df = pd.DataFrame([x.split('.')[0] for x in resultdf.NAME.values], columns=['Predictions'])
df['Scores'] = resultdf.GDT_TS.values
df.sort(columns=['Scores', 'Predictions'], ascending=[0, 1], inplace=True)
df.to_csv(args.output, index=False, quoting=csv.QUOTE_NONNUMERIC)
