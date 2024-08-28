#!/bin/bash

# Enable globstar for recursive file matching
shopt -s globstar

# Clear or create the energies.temp file
: > energies_raw.temp

# Loop over each .out file in subdirectories
for file in GaussOutputFiles/*.out
do 
    
    # Get name of a compound
    base_filename=$(basename "$file" .out)
    # Output filename without extension followed by the result of formatting
    echo -n "$base_filename " >> energies_raw.temp
    
    # Get corrected Gibbs Free value 
    $thermochange/formatters/formatted_energy_outputter.sh "$file" "$1" "$2" >> energies_raw.temp
done

awk '{print $1,$5}' energies_raw.temp > energies.temp

# Process the data with the Python script
python3 -W ignore "$(dirname "$0")/calculating_G_for_microkinetics.py" energies.temp $1 $2

# Clean up temporary file
rm energies.temp energies_raw.temp
