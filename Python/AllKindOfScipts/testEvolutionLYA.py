# -*- coding: utf-8 -*-
#
# created by Hélion du Mas des Bourboux
# < helion.du-mas-des-bourboux@cea.fr >

import myTools
from const_delta import *

import subprocess
import time
import sys
import numpy
import matplotlib.pyplot as plt
from iminuit import Minuit


b=1.6
def a(z):
	a0 = -7.04712
	a1 =  2.88486
	a2 = -0.300941
	return numpy.power(10., a0 + a1*z +a2*z*z)

xxx = numpy.arange(0., 6., 0.01)
yyy = a(xxx)

plt.plot(xxx,yyy)
myTools.deal_with_plot(False,False,False)
plt.show()


z_all = [2.38492854474]
dataOverSimu = [1.2]
z            = [3.16418109371]
dataOverSimu = [0.826889908268]
