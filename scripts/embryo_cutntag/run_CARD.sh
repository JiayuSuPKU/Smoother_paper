#!/bin/bash
datadir=$1
resdir=$2
assay_resolution=$3
sudo chmod -R 777 $datadir
sudo chmod -R 777 $resdir
sudo -u rstudio Rscript CARD_deconv.R $1 $2 $3