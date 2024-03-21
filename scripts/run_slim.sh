#!/usr/bin/env bash

# This script runs the SLiM simulation with the given selection coefficient and derived allele frequency

# Generate a random ID for this simulation run to ensure unique output filenames.
# This ID is a random number between 1 and 1,000,000.
id=`echo $((1 + $RANDOM % 1000000))`

# Store the first command-line argument as the selection coefficient.
sel=`echo $1`

# Store the second command-line argument as the derived allele frequency.
daf=`echo $2`

# Define the output filename using the selection coefficient, derived allele frequency, and the random ID.
output=sim_sel${sel}_daf${daf}_id${id}.ts

# Loop until the output file is successfully created.
while [ ! -f "$output" ]
do
  # Run the SLiM simulation with the specified parameters.
  # -d option is used to define variable values within the SLiM script.
  # Tgen=5000 defines the number of generations.
  # sel_coeff is the selection coefficient passed to the simulation.
  # OutTree specifies the name of the output tree sequence file.
  # end_freq is the desired end frequency of the derived allele.
  # simulate.slim is the SLiM model script to be executed.
  slim -d "Tgen=5000" -d "sel_coeff=${sel}" -d "OutTree='sim_sel${sel}_daf${daf}_id${id}.ts'" -d "end_freq=${daf}" simulate.slim
done

