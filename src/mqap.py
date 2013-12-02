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

from __future__ import print_function

import argparse
import sys

from generateTraining import generateFeatures
from train import train
from predict import predict

# we are using argparse module to implement argument parsing
parser = argparse.ArgumentParser()

# argparse can be used to implement git-style subcommand based arguments
subparsers = parser.add_subparsers(title='Subcommands',
                                   description='This program requires a '
             'mandatory subcommand (like git). See "%s <command> -h" to '
             'read about a specific subcommand.' % sys.argv[0],
             dest='subcommand')

# generatetraining subcommand ############################
generatetraining = subparsers.add_parser('generatetraining',
                                         help='Generates training set')

generatetraining.add_argument('-t', '--target', help='Target structure in '
                              'PDB format.', required=True)

generatetraining.add_argument('-m', '--models', help='Model structures in PDB '
                              'format. Multiple files can be separated by '
                              'comma.', required=True)

generatetraining.add_argument('-o', '--output',
                              help='Name of the output file in CSV format.',
                              required=True)

generatetraining.add_argument('-d', '--distance',
                              help='Maximum distance between target and model'
                              ' atoms to determine class labels.'
                              ' Default value is 3.5 angstrom.', default=3.5)

# train subcommand ############################
trainCommand = subparsers.add_parser('train', help='Trains a RandomForest '
                                     'based on given training set and '
                                     'produces a Pickle file to be used in '
                                     'prediction step.')

trainCommand.add_argument('-t', '--training',
                          help='Training set in CSV format. Class labels must'
                          ' be at the last column.', required=True)

trainCommand.add_argument('-o', '--output', help='RandomForest output file in '
                          'pickle format.', required=True)

trainCommand.add_argument('-n', '--ntrees', help='Number of trees in the '
                          'forest. Default: 100', default=100)

# predict subcommand ############################
predictCommand = subparsers.add_parser('predict', help='Predicts the quality '
                                       'of given model file in PDB format.')

predictCommand.add_argument('-r', '--randomforest',
                            help='RandomForest input file in pickle format.',
                            required=True)
predictCommand.add_argument('-m', '--model',
                            help='Model file to be predicted in PDB format',
                            required=True)
predictCommand.add_argument('-o', '--output',
                            help='Model quality output file in CSV format.',
                            required=True)

args = parser.parse_args()


if args.subcommand == 'generatetraining':
    generateFeatures(args.models, target=args.target, output=args.output,
                     distance=float(args.distance))

elif args.subcommand == 'train':
    train(args.training, int(args.ntrees), args.output)

elif args.subcommand == 'predict':
    predict(args.randomforest, args.model, args.output)

else:
    sys.exit('Invalid subcommand.')
