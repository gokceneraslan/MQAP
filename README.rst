=========================================================================================
Fast and accurate method for assessment of model quality for protein structure prediction
=========================================================================================

MQAP is a classifier-based quality assessment method utilizing local structural features of protein models (protein secondary structure, solvent accessibility, spatial neighborhood), using already available bioinformatic tools (STRIDE, PSIPRED, netSurfP, NACCESS) and machine learning methods (RandomForestClassifier class from scikit-learn framework).

Command line interface:
=======================

.. code-block:: console

    $ ./mqap -h

    usage: mqap.py [-h] {generatetraining,train,predict} ...

    optional arguments:
      -h, --help            show this help message and exit

    Subcommands:
      This program requires a mandatory subcommand (like git). See "./mqap.py
      <command> -h" to read about a specific subcommand.

      {generatetraining,train,predict}
        generatetraining    Generates training set
        train               Trains a RandomForest based on given training set and
                            produces a Pickle file to be used in prediction step.
        predict             Predicts the quality of given model file in PDB
                            format.

First argument is a mandatory sub-command, similar to that of git or
subversion CLI. :code:`generatetraining` command produces a CSV file out of a target 
PDB file and its predictions which are also in PDB format.

.. code-block:: console

    $ ./mqap generatetraining -h

    usage: mqap.py generatetraining [-h] -t TARGET -m MODELS -o OUTPUT
                                    [-d DISTANCE]

    optional arguments:
      -h, --help            show this help message and exit
      -t TARGET, --target TARGET
                            Target structure in PDB format.
      -m MODELS, --models MODELS
                            Model structures in PDB format. Multiple files can be
                            separated by comma.
      -o OUTPUT, --output OUTPUT
                            Name of the output file in CSV format.
      -d DISTANCE, --distance DISTANCE
                            Maximum distance between target and model atoms to
                            determine class labels. Default value is 3.5 angstrom.

Multiple model files can be specified using a comma-separated list of file names. 
This command first aligns all model structures onto the target structure and
computes local structure features such as surface accessibility and secondary
structure of residues. Output is saved into a CSV so that :code:`train` command can
parse it. Class labels are positive if residues of the aligned model file is
close enough to those of target file. Maximum distance is specified by
:code:`--distance/-d` parameter which has a default value of 3.5 Angstrom.

:code:`train` command parses given input CSV training file and trains a RandomForest
using scikit-learn framework.


.. code-block:: console

    $ ./mqap train -h

    usage: mqap.py train [-h] -t TRAINING -o OUTPUT [-n NTREES]

    optional arguments:
      -h, --help            show this help message and exit
      -t TRAINING, --training TRAINING
                            Training set in CSV format. Class labels must be at
                            the last column.
      -o OUTPUT, --output OUTPUT
                            RandomForest output file in pickle format.
      -n NTREES, --ntrees NTREES
                            Number of trees in the forest. Default: 100

Trained RandomForest instance is saved into the output file in pickle format
so that the next step (:code:`predict` command) can use it to predict assessment
scores of a new model file in PDB format. Number of trees in the random forest
can also be specified using :code:`--ntrees/-n` parameter.

:code:`predict` command requires a new model file of which assesment scores will be
computed. Another required parameter is the random forest file in pickle
format. Output is the structural features as well as the predicted class
labels.

.. code-block:: console

    $ ./mqap predict -h

    usage: mqap.py predict [-h] -r RANDOMFOREST -m MODEL -o OUTPUT

    optional arguments:
      -h, --help            show this help message and exit
      -r RANDOMFOREST, --randomforest RANDOMFOREST
                            RandomForest input file in pickle format.
      -m MODEL, --model MODEL
                            Model file to be predicted in PDB format
      -o OUTPUT, --output OUTPUT
                            Model quality output file in CSV format.

Requirements:
=============

MQAP requires following Python packages:

- `pandas <http://pandas.pydata.org/>`_ (for generating DataFrames and loading from/saving to CSV files conveniently)
- `scikit-learn <http://scikit-learn.org/>`_ (for classifying using RandomForestClassifier class)
- numpy (for usual matrix operations)
- `ProDy <http://www.csb.pitt.edu/prody/>`_ (for PDB parsing, alingment and manipulation)


