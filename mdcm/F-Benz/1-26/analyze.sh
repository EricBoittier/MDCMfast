
ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin
HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=$HOME/mdcm_fast/notebooks/dc.xyz
PCUBE=/pchem-data/meuwly/boittier/home/ref/optimised-cube/fluro/esp_fluro-ben-fixed.cube
DCUBE=/pchem-data/meuwly/boittier/home/ref/optimised-cube/fluro/esp_fluro-ben-fixed.cube
MAXATMCHG=60
NTRY=1
QXYZ=dc.xyz
FRAMES=frames.txt
MODEL=model.mdcm
INITIALXYZ=dc.xyz
OUTPUT=26_charges_refined.xyz
OUTPUTCUBE=26charges.cube
NAME=$INITIALXYZ
$BINDIR/pcubefit.x -generate -xyz $INITIALXYZ -esp $PCUBE -dens $DCUBE -v > ${NAME}.out
# Examine quality of fitted charges by comparing newly fitted model and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $OUTPUTCUBE -dens $DCUBE > analyze-cube-${NAME}.log
    