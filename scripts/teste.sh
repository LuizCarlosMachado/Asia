#!/usr/bin/env bash

#random=`echo $((1 + $RANDOM % 1000000))`
#id=`echo $1`
#daf=`echo $2`
#sel=`echo $3`
#output = sim_s${sel}_${random}.ts

while getopts u:a:f: flag
do
    case "${flag}" in
        i) id=`echo $1;;
        d) daf=`echo $2;;
        s) sel=`echo $3;;
    esac
done

slim -d "Tgen=5000" -d "sel_coeff=${sel}" -d "OutTree='sim_s${sel}_${id}.ts'" -d "end_freq=${daf}" simulate.slim




