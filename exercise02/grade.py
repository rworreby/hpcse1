#!/usr/bin/env python
# File       : grade.py
# Description: Generate grading submission file
# Copyright 2020 ETH Zurich. All Rights Reserved.
question_conf = {
        'Name' : 'Exercise 2',
        'Questions' : {
            'Question 1': {'Total Points': 35},
            'Question 2': {'Total Points': 20},
            'Question 3': {'Total Points': 20}
            }
        }

import argparse
import datetime
import sys

def parse_args():
    parser = argparse.ArgumentParser()

    for i in range(1, len(question_conf['Questions'])+1, 1):
        parser.add_argument('-q{:d}'.format(i),'--question{:d}'.format(i),
                type=int, default=0,
                help='Scored points for Question {:d}'.format(i))
        parser.add_argument('-c{:d}'.format(i),'--comment{:d}'.format(i),
                type=str, action='append', nargs='*',
                help='Comments for Question {:d} (you can add multiple comments)'.format(i))

    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parse_args()

    grade = lambda s,m: 3.0 + (6.0-3.0) * float(s)/m
    summary = {}
    score = 0
    maxpoints = 0
    header = '{name:s}: {date:s}\n'.format(
        name = question_conf['Name'], date = str(datetime.datetime.now()))
    width = len(header.rstrip())
    summary[0] = [header]
    for i in range(1, len(question_conf['Questions'])+1, 1):
        content = []
        qscore  = args['question{:d}'.format(i)]
        qmax    = question_conf['Questions']['Question {:d}'.format(i)]['Total Points']
        qscore  = max(0 , min(qscore, qmax))
        content.append( 'Question {id:d}: {score:d}/{max:d}\n'.format(
            id = i, score = qscore, max = qmax)
            )
        comments = args['comment{:d}'.format(i)]
        if comments is not None:
            for j,comment in enumerate([s for x in comments for s in x]):
                content.append( ' -Comment {id:d}: {issue:s}\n'.format(
                    id = j+1, issue = comment.strip())
                    )
        for line in content:
            width = width if len(line.rstrip())<width else len(line.rstrip())
        score += qscore
        maxpoints += qmax
        summary[i] = content
    assert maxpoints > 0
    with open('grade.txt', 'w') as out:
        out.write(width*'*'+'\n')
        for lines in summary.values():
            for line in lines:
                out.write(line)
            out.write(width*'*'+'\n')
        out.write('Grade: {:.2f}'.format(grade(score, maxpoints)))
