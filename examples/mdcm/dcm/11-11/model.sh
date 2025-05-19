
ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin
HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=$HOME/mdcm_fast/notebooks/dc.xyz
PCUBE=/pchem-data/meuwly/boittier/home/dcm/esp-dcm.cube
DCUBE=/pchem-data/meuwly/boittier/home/dcm/esp-dcm.cube
MAXATMCHG=60
NTRY=1
QXYZ=dc.xyz
FRAMES=frames.txt
MODEL=model.mdcm
INITIALXYZ=dc.xyz
OUTPUT=11_charges_refined.xyz
OUTPUTCUBE=11charges.cube
NAME=$INITIALXYZ
$BINDIR/comb-xyz-to-dcm.pl $QXYZ $PCUBE $FRAMES $MODEL