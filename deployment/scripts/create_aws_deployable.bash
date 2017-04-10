#!/bin/bash

if [ $# -eq 0 ]
  then
    echo ''
    echo "Invalid input parameters"
    echo "Usage: ./create_aws_deployable.bash NAME_OF_DESTINATION_FOLDER_TO_CREATE"
    echo ''
    exit 1
fi

DEPLOYABLE_FOLDER=$1

mkdir $DEPLOYABLE_FOLDER
mkdir $DEPLOYABLE_FOLDER/submissions
mkdir $DEPLOYABLE_FOLDER/.ebextensions

cp *.py $DEPLOYABLE_FOLDER/
cp deployment/config.txt $DEPLOYABLE_FOLDER/
cp requirements.txt $DEPLOYABLE_FOLDER/
cp deployment/config.txt $DEPLOYABLE_FOLDER/
cp -rf deployment/ebextensions/*.config $DEPLOYABLE_FOLDER/.ebextensions/
cp deployment/scripts/deploy_to_amazon.bash $DEPLOYABLE_FOLDER/

cp -rf static $DEPLOYABLE_FOLDER/
cp -rf templates $DEPLOYABLE_FOLDER/
cp -rf transdbs $DEPLOYABLE_FOLDER/
cp -rf cava $DEPLOYABLE_FOLDER/
