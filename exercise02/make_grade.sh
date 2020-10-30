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
    --question1 35 \
    --comment1 'Perfect' \
    --question2 20 \
    --comment2 'Solution is different from master solution but still correct.' \
    --question3 18 \
    --comment3 '-5 for missing the need of a barrier after do_work in task 3)a). Bonus of 3 points for finding an issue that was not caught by the master solution.' \
