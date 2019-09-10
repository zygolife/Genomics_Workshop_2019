#!/usr/bin/bash
#SBATCH -N 1 -n 1 -p short --mem 1gb  --out logs/init.log

mkdir -p input logs
pushd input
ln -s /bigdata/stajichlab/shared/tmp/Afum/DRR022920/* .
ln -s /bigdata/stajichlab/shared/tmp/Afum/DRR022919/* .
popd
