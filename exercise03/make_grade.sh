#!/usr/bin/env bash
# File       : make_grade.sh
# Description: Generate your exercise grade
# Copyright 2020 ETH Zurich. All Rights Reserved.
#
# EXAMPLE:
# python grade.py \
#     --question1 10 \ # my scored points
#     --comment1 'Add a comment to justify your score' \ # optional
#     --comment1 'Add a comment if you had problems somewhere' \ # optional
#     --question2 50 \ # my scored points
#
# FOR HELP:
# python grade.py --help
#
# The script generates a grade.txt file. Submit your grade on Moodle:
# https://moodle-app2.let.ethz.ch/course/view.php?id=13666

python3 grade.py \
    --question1 60 \
    --comment1 'Perfect' \
    --question2 37 \
    --comment2 'No explanation given for the different convergence of the eigenvalues (oscillation in the weight updates): -3p.' \
