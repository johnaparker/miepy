#!/bin/bash

curl https://refractiveindex.info/download/database/rii-database-2017-09-05.zip -o miepy/materials/database.zip
unzip miepy/materials/database.zip -d miepy/materials
rm miepy/materials/database.zip
