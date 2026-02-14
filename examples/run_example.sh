#!/bin/bash

# Create example directories
mkdir -p formula-notes
mkdir -p chapters_output

# Create a sample page
cat <<EOF > formula-notes/page_sample.txt
Chapter 15: Heat Exchangers
15-1 Heat Transfer Rate
Q = U * A * delta_T
where:
Q := heat transfer rate
U := overall heat transfer coefficient
A := heat transfer area
delta_T := temperature difference
EOF

echo "Running vakyume_begin OCR tool..."
python3 ../vakyume_begin.py --input formula-notes --output chapters_output

echo "Results in examples/chapters_output:"
ls chapters_output
cat chapters_output/page_sample.py
