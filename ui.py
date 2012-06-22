from __future__ import print_function

import sys



def _do_print(x, newline, color, file):
	msg = "\033[00;%im%s\033[00m%s" % (color, str(x), newline and "\n" or "")
	file.write(msg)
def print_msg(x = "", newline = True, file = sys.stdout):
	file.write(x + ("\n" if newline else ""))
def print_info(x = "", newline = True, file = sys.stdout):
	_do_print(x, newline, 32, file)
def print_err(x = "", newline = True, file = sys.stdout):
	_do_print(x, newline, 31, file)
def print_warn(x = "", newline = True, file = sys.stdout):
	_do_print(x, newline, 33, file)
def print_disabled(x = "", newline = True, file = sys.stdout):
	_do_print(x, newline, 30, file)
