
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
OUTPUT=48_charges_refined.xyz
OUTPUTCUBE=48charges.cube
NAME=$INITIALXYZ

    $BINDIR/pcubefit.x -xyz $INITIALXYZ $INITIALXYZ -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v
    