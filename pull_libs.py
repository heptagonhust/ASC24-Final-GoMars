#!/usr/bin/env python3

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Pulls all or some submodules from external repositories.')
parser.add_argument('--libs', nargs='+', help='The libraries to pull.', choices=('noahmp','rrtmgp'), default=[])
args = parser.parse_args()

def run(cmd):
	print(f'==> {cmd}')
	res = subprocess.run(cmd, shell=True, check=True)

def pull(lib_path, recursive=False):
	if recursive:
		cmd = f'git submodule update --recursive --init {lib_path}'
	else:
		cmd = f'git submodule update --init {lib_path}'
	run(cmd)
	
pull('lib/container')
pull('lib/datetime')
pull('lib/fiona')
pull('lib/flogger')
pull('lib/string')

project_root = os.getcwd()

if 'noahmp' in args.libs:
	pull('lib/noahmp')
	os.chdir('lib/noahmp')
	run('git apply ../noahmp.diff')
	os.chdir(project_root)

if 'rrtmgp' in args.libs:
	pull('lib/rrtmgp')
	os.chdir('lib/rrtmgp')
	run('git apply ../rrtmgp.diff')
	os.chdir(project_root)
