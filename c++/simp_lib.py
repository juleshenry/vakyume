import ctypes
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
handle = ctypes.CDLL(dir_path + "/libTest.so")

handle.my_Function.argtypes = [ctypes.c_int]


def my_function(num):
    return handle.my_function(num)
