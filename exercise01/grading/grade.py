#!/usr/bin/env python
# File       : grade.py
# Description: Generate grading submission file
# Copyright 2020 ETH Zurich. All Rights Reserved.
question_conf = {
  'Name' : 'Exercise 1',
  'Questions' : {
    'Question 1': {'Total Points': 5},
    'Question 2': {'Total Points': 10},
    'Question 3': {'Total Points': 10},
    }
}

import argparse
import datetime
import sys

def parse_args():
  parser = argparse.ArgumentParser()

  for i in range(1, len(question_conf['Questions'])+1, 1):

    parser.add_argument( f'-q{i}',f'--question{i}',
                         type=float, default=0, help=f'Scored points for Question {i}')

    parser.add_argument( f'-c{i}',f'--comment{i}',
                         type=str, action='append', nargs='*',
                         help=f'Comments for Question {i} (you can add multiple comments)')

  return vars(parser.parse_args())

if __name__ == "__main__":
  args = parse_args()

  grade = lambda s,m: 3.0 + (6.0-3.0) * float(s)/m

  summary = {}
  score = 0
  maxpoints = 0

  name = question_conf['Name']
  date = str(datetime.datetime.now())
  header = f'{name}: {date} \n'
  width = len(header.rstrip())
  summary[0] = [header]


  for i in range(1, len(question_conf['Questions'])+1, 1):
    content = []
    qscore  = args[f'question{i}']
    qmax    = question_conf['Questions'][f'Question {i}']['Total Points']
    qscore  = max(0 , min(qscore, qmax))
    content.append( f'Question {i}: {qscore}/{qmax}\n')

    comments = args[f'comment{i}']
    if comments is not None:
      for j,comment in enumerate([s for x in comments for s in x]):
        content.append( f' -Comment {i}-{j+1}: {comment.strip()}\n')
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
      out.write('Grade: {:.2f}\n'.format(grade(score, maxpoints)))
