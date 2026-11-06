#!/bin/bash
TOTAL_ATOMS=2048
NEW_N=8

find . -type f -name "brawl.inp" | while read -r file; do
    echo "Processing $file"

    # Update n_1, n_2, n_3
    sed -i "s/^n_1=.*/n_1=$NEW_N/" "$file"
    sed -i "s/^n_2=.*/n_2=$NEW_N/" "$file"
    sed -i "s/^n_3=.*/n_3=$NEW_N/" "$file"

    # Read concentrations (keep spaces)
    conc=$(grep "^!species_concentrations" "$file" | awk -F= '{gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}')

    # Compute integer species numbers from concentrations
    new_nums=$(python3 - <<EOF
import math
conc = [float(x) for x in "$conc".split()]
total = $TOTAL_ATOMS
floors = [math.floor(c*total) for c in conc]
remainder = total - sum(floors)
frac = [(c*total - math.floor(c*total), i) for i, c in enumerate(conc)]
for _, i in sorted(frac, reverse=True)[:remainder]:
    floors[i] += 1
print(" ".join(str(x) for x in floors))
EOF
)

    # Replace species_numbers line
    sed -i "s/^species_numbers.*/species_numbers = $new_nums/" "$file"
done

