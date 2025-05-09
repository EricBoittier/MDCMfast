
ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin

HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=$HOME/mdcm_fast/notebooks/dc.xyz
PCUBE=$HOME/acem/ESP.cube
DCUBE=$HOME/acem/Density.cube
MAXATMCHG=60
NTRY=1

QXYZ=27_charges_refined.xyz
FRAMES=frames.txt
OUTPUT=test.mdcm

$BINDIR/comb-xyz-to-dcm.pl $QXYZ $PCUBE $FRAMES $OUTPUT


