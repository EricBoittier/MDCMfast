ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin

HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=dc.xyz
PCUBE=$HOME/ref/optimised-cube/chloro/esp_chloro-ben.cube
DCUBE=$HOME/ref/optimised-cube/chloro/esp_chloro-ben.cube
MAXATMCHG=60
NTRY=1

$BINDIR/pcubefit.x -xyz $INITIALXYZ $INITIALXYZ -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v 
