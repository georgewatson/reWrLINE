#!/bin/sh

name=$1
top=$2
traj=$3

mkdir -p $name

cpptraj <<EOF
parm $top
trajin $traj
strip !(@C1')
trajout $name/C.mdcrd
EOF
