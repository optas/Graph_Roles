'''
Created on Feb 16, 2015

Contact: pachlioptas@gmail.com

Copyright notice: 
Copyright (c) 2015, Panagiotis Achlioptas
You are free to use, change, or redistribute this code in any way you want for non-commercial purposes only.
'''

import time
import math
import smtplib
import numpy as np

def p(message=""):
    '''
    shortcut for raw_input()
    '''
    return raw_input(message)


def is_increasing(l):
    return all(l[i] <= l[i+1] for i in xrange(len(l) - 1))


def is_decreasing(l):
    return all(l[i] >= l[i+1] for i in xrange(len(l) - 1))


def is_monotone(l):
    n = len(l) - 1
    if l[0] < l[n]:   
        return is_increasing(l)
    elif l[0] > l[n]:     
        return is_decreasing(l)
    else:
        return all(l[i] == l[i+1] for i in xrange(len(l) - 1))

     
def average_time(kRuns, testFunction, *args, **kwargs):
    """Computes the average time of execution of -testFunction- in -kRuns- executions.
    Input: -kRuns- number of times that we will execute the -testFunction- 
    """
        
    t = 0
    for i in xrange(kRuns):
        t0 = time()
        testFunction(*args, **kwargs)
        t1 = time()
        t += (t1 - t0)
    return t / float(kRuns)


def print_first_k(iterable, k):    
    '''
    print the first k 'sub'-objects stored in each object of an iterable.
    '''
    assert(k>0)
    i = 0    
    for item in iterable:
        print item,
        i += 1
        if i == k:
            print ""
            i = 0

def print_fancy(inString, decoration="-", delim=None):
    emptySpace = 3
    inLen = len(inString) + 2*emptySpace

    frameLen = max(40, inLen)
    diff = frameLen - inLen
    if diff > 0:
        center =  int(math.ceil(diff/float(2)))
    else:
        center = 0

    print "\n", decoration * frameLen
    print decoration * frameLen
    if center != 0:
        print "*" *(center-emptySpace), " "*emptySpace, inString, " "*emptySpace, "*" * (center-emptySpace)
    else:
        print inString
    print decoration * frameLen
    print decoration * frameLen , "\n"


def sendemail(from_addr, to_addr_list, cc_addr_list, subject, message, login, password, smtpserver='smtp.gmail.com:587'):
    header  = 'From: %s\n' % from_addr
    header += 'To: %s\n' % ','.join(to_addr_list)
    header += 'Cc: %s\n' % ','.join(cc_addr_list)
    header += 'Subject: %s\n\n' % subject
    message = header + message
 
    server = smtplib.SMTP(smtpserver)
    server.starttls()
    server.login(login,password)
    problems = server.sendmail(from_addr, to_addr_list, message)
    server.quit()


def emailPanos(message="Hi love.", subject="Panos emails you."):   
    sendemail(['eukolomnimonitos@gmail.com'], ['optas@stanford.edu'], cc_addr_list=[],
              subject=subject, message=message,
              login="eukolomnimonitos", password="Liquid_Sun",
              smtpserver='smtp.gmail.com:587')
    

def relative_to_worst(observations):
    '''
    observations = [0, 1, 2]
    In [36]: relative_to_worst(observations)
    Out[36]: [0, 0.0, 1.0]
    '''
    obs = np.array(sorted(observations))
    minNonZero  = np.where(obs!=0)[0][0]
    minNonZero  = float(obs[minNonZero])

    res = []
    for i in xrange(len(observations)):
        if observations[i] == 0 and minNonZero > 0: #alternatively there exist a negative minimum 
            res.append(0)
        else:
            res.append(abs(observations[i] - minNonZero)/abs(minNonZero))
    return res
