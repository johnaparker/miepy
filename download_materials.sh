#!/bin/bash

wget https://refractiveindex.info/download/database/rii-database-2017-09-05.zip -O miepy/materials/database.zip
unzip miepy/materials/database.zip -d miepy/materials
rm miepy/materials/database.zip
