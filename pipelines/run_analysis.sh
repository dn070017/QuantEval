#!/usr/bin/env bash

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

$SCRIPTDIR/download_dataset.sh
$SCRIPTDIR/prepare_reference.sh
$SCRIPTDIR/read_simulation.sh
$SCRIPTDIR/preprocessing.sh
$SCRIPTDIR/assembly.sh
$SCRIPTDIR/quantification.sh
$SCRIPTDIR/postprocessing.sh
$SCRIPTDIR/QuantEval.sh 
