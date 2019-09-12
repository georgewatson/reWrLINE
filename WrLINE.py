#! /usr/bin/env python3

import sys
import os
import writhe
import caxislib

name = sys.argv[1]
top = sys.argv[2]
traj = sys.argv[3]
num_bp = int(sys.argv[4])
num_steps = int(sys.argv[5])

# Strip trajectory to get C1' coordinates
os.system(' '.join(['bash stripC.sh', name, top, traj]))

print(f"Processing {name}")
print("Reading files & initialising arrays")
strand_a, strand_b, midpoints = caxislib.read(name, num_bp, num_steps)

print("Calculating first-order helical axis")
helix_axis = caxislib.helix_axis(num_bp, num_steps, midpoints, strand_a)

print("Calculating twist")
twist = caxislib.full_twist(name, num_bp, num_steps, strand_a, strand_b,
                            helix_axis)

print("Calculating helical axis")
caxis = caxislib.caxis(name, num_bp, num_steps, midpoints, twist)

print("Calculating register angles")
sinreg = caxislib.sinreg(name, num_bp, num_steps, midpoints, caxis)

print("Writing output .xyz and .3col files")
caxislib.make_files(name, num_bp, num_steps, midpoints, caxis)

wr = writhe.main(name, num_bp, num_steps)

print(f"Job {name} done!")
