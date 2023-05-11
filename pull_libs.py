#!/usr/bin/env python3

import argparse
import subprocess

parser = argparse.ArgumentParser(description='Pulls all or some submodules from external repositories.')
parser.add_argument('--libs', nargs='+', help='The libraries to pull.', choices=('ccpp',), default=[])
args = parser.parse_args()

def pull(lib_path, resursive=False):
	if resursive:
		cmd = f'git submodule update --resursive --init {lib_path}'
	else:
		cmd = f'git submodule update --init {lib_path}'
	print(f'==> {cmd}')
	res = subprocess.run(cmd, shell=True, check=True)
	
pull('lib/container')
pull('lib/datetime')
pull('lib/fiona')
pull('lib/flogger')
pull('lib/string')

if 'ccpp' in args.libs:
	pull('src/physics/ccpp', resursive=True)
