# Note: this only needs to be done when a new type of support
# structure/foundation is desired.
for i in {0001..0100}; do
    convert /Users/jpanetta/Research/microstructures/isosurface_inflator/test_unzip/test_unzip\\unopt_layout_50xy_50z_$i.bmp -crop 20x20+823+691 +repage support_$i.png;
done
for i in {0001..200}; do
    convert /Users/jpanetta/Research/microstructures/isosurface_inflator/test_unzip/test_unzip\\unopt_layout_50xy_50z_$i.bmp -crop 157x157+488+275 + repage corner_$i.png;
done
