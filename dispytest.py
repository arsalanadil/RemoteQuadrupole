#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:55:32 2017

@author: arsalanadil
"""
import time
def compute(n):
   time.sleep(n)
   return n

if __name__ == '__main__':
    import dispy, logging
    nodes = ['10.1.1.254']
    cluster = dispy.JobCluster(compute, nodes=nodes, ext_ip_addr='127.0.0.1', port=51348, loglevel=logging.DEBUG)
    job = cluster.submit(2)
    print('result: %s' % job())