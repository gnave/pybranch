# This file is to perform grep using python code 

def grep_open(grep_key, open_file):
	read_lines = open(open_file, 'r').readlines()
	return grep_lines(grep_key, read_lines)

def grep_lines(grep_key, read_lines):
	return_lines = []
	for line in read_lines:
		if grep_key in line:
			return_lines.append(line)
	return return_lines