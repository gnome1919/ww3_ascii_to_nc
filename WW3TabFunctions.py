# -*- coding: utf-8 -*-

# This file contains functions which are used in ProjectMain.py(w)
# Copyright 2017 Farrokh A. Ghavanini <ghavanini@gmail.com>

#
#

from numpy import array

def extractParam(filename, first_line, last_line, param_no):
    """Extract a parameter values and save it in a list"""
    with open(filename, 'r') as inputfile:
        list_container = []
        for i, line in enumerate(inputfile, 1):
            if first_line <= i <= last_line:
                list_container.append(float(line.split()[param_no]))
    return list_container


def extractLine(filename, line_no):
    """Extracts a specific line of an ascii file"""
    with open(filename, 'r') as inputfile:
        for i, line in enumerate(inputfile, 1):
            if i == line_no:
                break
        header = line
    return header


def timestepCounter(filename):
    """Finds number of timesteps by counting 'Time' string in tab file."""
    with open(filename, 'r') as inputfile:
        count = 0
        for line in inputfile:
            if 'Time' in line:
                count += 1
    return count


def blockIdentifier(filename):
    """Identifies the block of data for each timestep in tab file"""
    try:        
        with open(filename, 'r') as inputfile:
            newline_counter = 0
            for i, line in enumerate(inputfile):
                if line.lstrip().startswith("Time"):
                    first_line = i + 5
                if line.lstrip() == "":
                    newline_counter = newline_counter + 1
                if newline_counter == 2:
                    last_line = i - 1
                    break
        return first_line, last_line
    except FileNotFoundError:
        return 0, 0
    except UnboundLocalError:
        return 0, 1
        

def listToArray(my_list):
    """Taking a list and save it in a numpy array"""
    my_array = array(my_list)
    return my_array
