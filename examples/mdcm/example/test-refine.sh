
ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin

HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=$HOME/mdcm_fast/notebooks/dc.xyz
PCUBE=$HOME/acem/ESP.cube
DCUBE=$HOME/acem/Density.cube
MAXATMCHG=60
NTRY=1

$BINDIR/pcubefit.x -xyz $INITIALXYZ $INITIALXYZ -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v 
