#!/bin/bash

mkdir -p jpg_output

for pdf in *.pdf; do
    name="${pdf%.pdf}"
    pdftoppm -jpeg -r 300 "$pdf" "jpg_output/${name}"
    echo "Converted: $pdf"
done

echo "All done! Files saved to jpg_output/"
