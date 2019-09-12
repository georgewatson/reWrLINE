#! /usr/bin/env python3

"""
        __        __    _     ___ _   _ _____
 _ __ __\ \      / / __| |   |_ _| \ | | ____|
| '__/ _ \ \ /\ / / '__| |    | ||  \| |  _|
| | |  __/\ V  V /| |  | |___ | || |\  | |___
|_|  \___| \_/\_/ |_|  |_____|___|_| \_|_____|

reWrLINE: A reimplementation of WrLINE

(c) 2019 George D. Watson, University of York
https://georgewatson.me

Based on WrLINE
by Thana Sutthibutpong, Sarah Harris, and Agnes Noy.
Please cite
Sutthibutpong T, Harris S A and Noy A 2015 J. Chem. Theory Comput. 11 2768-75
https://doi.org/10.1021/acs.jctc.5b00035

This is the test suite.
It is recommended that you run this before using this implementation of the
software.
Tests ensure consistency with WrLINE, not correctness.
"""

import caxislib
import filecmp
import numpy as np
from termcolor import colored

print(__doc__)
print("---\n")

tolerance = 0.1

name = 'test'
num_bp = 336
num_steps = 8

a1 = np.array([330.0, 50.0, 0.0])
a2 = np.array([330.0, 40.0, 10.0])
b1 = np.array([330.0, 45.0, 20.0])
b2 = np.array([335.0, 40.0, 15.0])
z = np.array([5.0, 0.0, -2.0])

strand_a, strand_b, midpoints = caxislib.read(name, num_bp, num_steps)
helix_axis = caxislib.helix_axis(num_bp, num_steps, midpoints, strand_a)
full_twist = caxislib.full_twist(name, num_bp, num_steps, strand_a, strand_b,
                                 helix_axis)
caxis = caxislib.caxis(name, num_bp, num_steps, midpoints, full_twist)

tests = {
    "rotate_to_x": [np.linalg.norm(caxislib.rotate_to_x(a1)), 1.7320508075689],
    "rotate_to_z": [np.linalg.norm(caxislib.rotate_to_z(a1)), 1.7320508075689],
    "twist 1\t": [caxislib.twist(a1, a2, b1, b2, z), 71.997556736987],
    "twist 2\t": [caxislib.twist(a1, b1, a2, b2, z), -15.06995574963],
    "twist 90deg": [caxislib.twist(np.array([0, 0, 0]),
                                   np.array([0, 1, 0]),
                                   np.array([0, 0, 0]),
                                   np.array([1, 0, 0]),
                                   np.array([0, 0, 1])), -90],
    "read strand_a": [sum(sum(sum(strand_a))), 1057248.34],
    "read strand_b": [sum(sum(sum(strand_b))), 1057277.65],
    "read midpoints": [sum(sum(sum(midpoints))), 1057262.995],
    "helix_axis": [sum(sum(sum(helix_axis))), 1057262.995],
    "full_twist": [sum(sum(full_twist)), 83325.202562320],
    "caxis\t": [sum(sum(sum(caxis))), 1057251.5419271],
    "C.mdcrd\t": [filecmp.cmp(f'{name}/C.mdcrd', f'{name}/C.mdcrd.original'),
                  True],
    "tw.ser\t": [filecmp.cmp(f'{name}/tw.ser', f'{name}/tw.ser.original'),
                 True]
}

pass_text = colored('[PASS]', 'green')
fail_text = colored('[FAIL]', 'red')

failures = []

print()
for test, assertion in tests.items():
    print(test, end="\t")
    try:
        if abs(assertion[1] - assertion[0]) < tolerance:
            print(pass_text)
        else:
            failures.append(test)
            print(fail_text)
            print(f"\tExpected\t{assertion[1]}")
            print(f"\tGot\t\t{assertion[0]}")
    except Exception as e:
        failures += 1
        print(fail_text)
        print("\tThe following exception was raised:")
        print(f"\t{e}")

print()
if failures:
    print(f"Total failures:\t{len(failures)}")
    for test in failures:
        print(f"\t{test}")
else:
    print("All tests passed!")
