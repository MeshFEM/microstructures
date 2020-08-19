#!/usr/bin/env zsh

./RemoveUnsupported $@

outDir=$(dirname $1)
./SliceSupporter $1 $outDir
# copy over the support slices
for i in $outDir/out_*.png; do
    convert $i -set colorspace RGB -depth 8 ${i%png}bmp
    rm $i;
done

# Copy over the modified object slices
for i in $@; do
     convert ${i%.bmp}.cleaned.png -set colorspace RGB -depth 8 $outDir/out_$((${$(basename $i .bmp)##*_}+28)).bmp
done

# Delete the original slices.
rm $@
rm $outDir/*.png
