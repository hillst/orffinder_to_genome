#!/bin/bash
tar -xjf testfiles.tar.bz2
../orf_to_genome.py -r 13-10-25_9509_annotation.gff3 -b OrfPredictor.pep -o output.gff3
