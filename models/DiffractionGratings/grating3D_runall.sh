GMSHDIR="/home/demesy/programs/gmsh-4.4.1-Linux64/bin"
for t in bisin checker halfellipsoid hole pyramid torus 2Dlamellar solarcell conv skew
do
    $GMSHDIR/gmsh grating3D.pro -setstring test_case $t -
    mv res3D res3D_$t
done

for t in bisin checker halfellipsoid hole pyramid torus 2Dlamellar solarcell conv skew
do
    $GMSHDIR/gmsh grating3D.geo -setstring test_case $t res3D_$t/*.pos &
    echo "____________\n"$t
    python grating3D_postplot.py res3D_$t
done