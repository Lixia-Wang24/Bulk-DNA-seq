#!/bin/bash

#fasterq-dump -e 8 --split-files SRR14514145.sra -O ./test1

fasterq-dump -e 8 --split-files --include-technical SRR14514145.sra -O ./test2
