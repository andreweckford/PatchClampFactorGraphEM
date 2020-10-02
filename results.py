#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 14:03:39 2020

@author: andreweckford
"""

from resultHandling.resultHandler import Results


def main():
    r = Results()
    r.parseResults()
    print(r.kpErrors())
    print(r.emErrors())

if __name__ == "__main__":
    main()