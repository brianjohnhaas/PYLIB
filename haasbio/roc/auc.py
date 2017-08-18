#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import unittest
import logging

logger = logging.getLogger(__name__)


def auc(ordered_pts_list):


    pts = ordered_pts_list[:]
    
    auc = 0

    ptA = pts.pop(0)
    
    while pts:
        ptB = pts.pop(0)

        x1 = ptA[0]
        y1 = ptA[1]

        assert(x1 >= 0)
        assert(y1 >= 0)

        x2 = ptB[0]
        y2 = ptB[1]

        assert(x2 >= 0)
        assert(y2 >= 0)

        delta_x = abs(x2-x1)
        delta_y = abs(y2-y1)

        min_y = min(y1, y2)

        base_trap_area = min_y * delta_x
        top_triangle_area = delta_y * delta_x / 2

        trap_area = base_trap_area + top_triangle_area

        if x2 > x1:
            auc += trap_area
        else:
            auc -= trap_area

        logger.debug("auc between pts {} and {}: {}".format(ptA, ptB, trap_area))
        
        # set up for next pt
        ptA = ptB

    logger.debug("auc total area: {}".format(auc))
    
    return auc




def test_auc_calc():
    logger.setLevel(logging.DEBUG)

    pts = [ (0,0), (1,1), (2,2) ]

    auc(pts)


if __name__ == '__main__':
    test_auc_calc()
