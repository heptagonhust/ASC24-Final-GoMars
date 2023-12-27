#!/usr/bin/env python3

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Pulls necessary static data from external repositories.')
parser.add_argument('-p', '--planet', help='Select the planet to run', choices=('earth', 'mars'), required=True)
args = parser.parse_args()

def run(cmd):
	print(f'==> {cmd}')
	res = subprocess.run(cmd, shell=True, check=True)

gmcore_root = os.path.dirname(os.path.realpath(__file__))
data_root = f'{gmcore_root}/data'
if not os.path.isdir(data_root):
	os.makedirs(data_root)
os.chdir(data_root)
if os.path.isdir(args.planet):
	os.remove(args.planet)

if args.planet == 'earth':
	pass
elif args.planet == 'mars':
	run('git clone https://gitee.com/dongli85/gomars-data mars')
