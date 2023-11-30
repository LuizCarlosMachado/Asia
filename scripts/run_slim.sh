#!/usr/bin/env bash

#id=`echo $1`
id=`echo $((1 + $RANDOM % 1000000))`
sel=`echo $1`
daf=`echo $2`
output=sim_sel${sel}_daf${daf}_id${id}.ts

while [ ! -f "$output" ]
do
  slim -d "Tgen=5000" -d "sel_coeff=${sel}" -d "OutTree='sim_sel${sel}_daf${daf}_id${id}.ts'" -d "end_freq=${daf}" simulate.slim
done


