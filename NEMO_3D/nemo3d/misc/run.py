#!/usr/local/bin/python

import threading
import popen2

script = "run.pbs"
command = "/home/packages/OpenPBS/bin/qsub ./%s" % script
(ignored, out) = popen2.popen2(command)
num = out.readline().split('.')[0]
file_name = "/home/fabiano/%s.o%s" % (script, num)
is_finished = 0
timer = threading.Event()

def check_finished():
    is_finished = "--- Done. ---" in open(file_name, 'r').readlines()
    
while not is_finished:
    print "Waiting..."
    timer.wait(30.0)
