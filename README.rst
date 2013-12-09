=========================================================================================
Fast and accurate method for assessment of model quality for protein structure prediction
=========================================================================================

MQAP is a classifier-based quality assessment method utilizing local structural features of protein models (protein secondary structure, solvent accessibility, spatial neighborhood), using already available bioinformatic tools (STRIDE, PSIPRED, netSurfP, NACCESS) and machine learning methods (RandomForestClassifier class from scikit-learn framework).

Command line interface:
=======================

.. code-block:: console

    $ ./mqap -h

    usage: mqap [-h] {generatetraining,train,predict,rank,cv} ...

    optional arguments:
      -h, --help            show this help message and exit

    Subcommands:
      This program requires a mandatory subcommand (like git). See "./mqap
      <command> -h" to read about a specific subcommand.

      {generatetraining,train,predict,rank,cv}
        generatetraining    Generates training set
        train               Trains a RandomForest based on given training set and
                            produces a Joblib file to be used in prediction step.
        predict             Predicts the quality of given model file in PDB
                            format.
        rank                Rank given prediction outputs in CSV format.
        cv                  Perform a k-fold cross validation using given CSV
                            file.

First argument is a mandatory sub-command, similar to that of git or
subversion CLI. :code:`generatetraining` command produces a CSV file out of a target 
PDB file and its predictions which are also in PDB format.

.. code-block:: console

    $ ./mqap generatetraining -h

    usage: mqap generatetraining [-h] -t TARGET -m MODEL [MODEL ...] -o OUTPUT
                                 [-d DISTANCE]

    optional arguments:
      -h, --help            show this help message and exit
      -t TARGET, --target TARGET
                            Target structure in PDB format.
      -m MODEL [MODEL ...], --model MODEL [MODEL ...]
                            Model structure(s) in PDB format. Multiple files can
                            be given.
      -o OUTPUT, --output OUTPUT
                            Name of the output file in CSV format.
      -d DISTANCE, --distance DISTANCE
                            Maximum distance between target and model atoms to
                            determine class labels. Default value is 3.5 angstrom.


Multiple model files can be specified using a list of file names. 
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
                            RandomForest output file in Joblib format.
      -n NTREES, --ntrees NTREES
                            Number of trees in the forest. Default: 100

Trained RandomForest instance is saved into the output file in Joblib format
so that the next step (:code:`predict` command) can use it to predict assessment
scores of a new model file in PDB format. Number of trees in the random forest
can also be specified using :code:`--ntrees/-n` parameter.

:code:`predict` command requires a new model file of which assesment scores will be
computed. Another required parameter is the random forest file in Joblib
format. Output is the structural features as well as the predicted class
labels.

.. code-block:: console

    $ ./mqap predict -h

    usage: mqap predict [-h] -r RANDOMFOREST -m MODEL [MODEL ...]
                        [-o [OUTPUT [OUTPUT ...]]]

    optional arguments:
      -h, --help            show this help message and exit
      -r RANDOMFOREST, --randomforest RANDOMFOREST
                            RandomForest input file in Joblib format.
      -m MODEL [MODEL ...], --model MODEL [MODEL ...]
                            Model file(s) to be predicted in PDB format
      -o [OUTPUT [OUTPUT ...]], --output [OUTPUT [OUTPUT ...]]
                            Model quality output file(s) in CSV format.If not
                            supplied, output file name will be based on model file
                            name.

:code:`rank` command computes a quality score per model file and sorts them
accordingly. Output is written to the supplied file.

.. code-block:: console

    $ ./mqap rank -h

    usage: mqap rank [-h] -p PREDICTION [PREDICTION ...] -o OUTPUT

    optional arguments:
      -h, --help            show this help message and exit
      -p PREDICTION [PREDICTION ...], --predictions PREDICTION [PREDICTION ...]
                            Prediction output file(s) produced by predict command.
      -o OUTPUT, --output OUTPUT
                            Ranking output file in CSV format.

:code:`cv` command performs a k-fold cross-validation on a given data set.
Data set is specified in CSV format. First column represents target file
whereas the rest of the columns represent corresponding model files for the
target given in the first column. Output is saved into given file. Also,
training sets, random forest files and prediction files can be saved into
given directories. Additionally, cutoff distance and ntrees arguments can be
supplied for training steps. Finally, previously saved trainingset and random
forest files can be reused in next runs through :code:`-u, --reuse` argument.

.. code-block:: console

    $ ./mqap cv -h

    usage: mqap cv [-h] -i INPUT [-k FOLD] -o OUTPUT [-t TRAININGDIR]
                   [-r RANDOMFORESTDIR] [-p PREDICTIONDIR] [-d DISTANCE]
                   [-n NTREES] [-u]

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            In the input file, first column represents target file
                            while the rest represents model files.
      -k FOLD, --fold FOLD  Fold size of k-fold CV. Default:5
      -o OUTPUT, --output OUTPUT
                            Cross-validation output file in CSV format.
      -t TRAININGDIR, --trainingdir TRAININGDIR
                            If specified, all training CSVs will be saved into
                            this dir.
      -r RANDOMFORESTDIR, --randomforestdir RANDOMFORESTDIR
                            If specified, all trained RFs will be saved into this
                            dir in Joblib format.
      -p PREDICTIONDIR, --predictiondir PREDICTIONDIR
                            If specified, all predictions will be saved into this
                            dir.
      -d DISTANCE, --distance DISTANCE
                            Maximum distance between target and model atoms to
                            determine class labels. Default value is 3.5 angstrom.
      -n NTREES, --ntrees NTREES
                            Number of trees in the forest. Default: 100
      -u, --reuse           If enabled, previously saved training/randomforest
                            files are reused. Default: Disabled

Requirements:
=============

MQAP requires following Python packages and command line tools:

- `pandas <http://pandas.pydata.org/>`_ (for generating DataFrames and loading from/saving to CSV files conveniently)
- `scikit-learn <http://scikit-learn.org/>`_ (for classifying using RandomForestClassifier class)
- numpy (for usual matrix operations)
- `ProDy <http://www.csb.pitt.edu/prody/>`_ (for PDB parsing, alignment and manipulation)
- `Joblib <http://pythonhosted.org/joblib/>`_ (for dumping/loading
  RandomForest files)
- `DSSP <http://swift.cmbi.ru.nl/gv/dssp/>`_ executable (for computing surface
  accessibility and secondary structure) which must be in :code:`$PATH` with name :code:`dssp`
- `STRIDE <http://webclu.bio.wzw.tum.de/stride/>`_ executable (for computing surface
  accessibility and secondary structure) which must be in :code:`$PATH` with name :code:`stride`
- `NetSurfP 1.0 <http://www.cbs.dtu.dk/cgi-bin/sw_request?netsurfp>`_
  executable (for computing surface exposure from sequence information) must
  be in :code:`$PATH` with name :code:`netsurfp`
- `PSIPRED 3.3 <http://bioinf.cs.ucl.ac.uk/software_downloads/>`_ executable
  must be in :code:`$PATH` with name :code:`runpsipred`

