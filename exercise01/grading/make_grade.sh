#!/usr/bin/env bash
# File       : make_grade.sh
# Description: Generate your exercise grade
# Copyright 2020 ETH Zurich. All Rights Reserved.
#
# EXAMPLE:
# python grade.py \
#     --question1 0 \
#     --comment1 'Add a comment to justify your score' \
#     --question2 0 \
#     --question2 0 \
#
# FOR HELP:
# python grade.py --help
#
# The script generates a grade.txt file. Submit your grade on Moodle:
# https://moodle-app2.let.ethz.ch/course/view.php?id=13666

# Note: --question1 and --question2 are not graded
python3 grade.py \
    --question1 5 \
    --comment1 'Made assumption that "more" is additional and not relative, but justified it in the assumptions. Therefore, I still award the full 0.5 points for that part.' \
    --question2 10 \
    --comment2 'Plot axis are suboptimal, change of x-axis would have been better, but still captures the essentials of the exercise.' \
    --question3 9.75 \
    --comment3 'Subtracted 0.25 for erronous behaviour of the random permutation code. Implementation seems correct. Argumentation is also correct.' \
