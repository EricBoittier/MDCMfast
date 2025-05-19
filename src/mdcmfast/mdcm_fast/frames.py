from io import StringIO
from pathlib import Path

import ase
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import rdkit
from dcmnet.plotting_3d import plot_3d_models, plot_3d_molecule
from IPython.display import HTML, Image, display
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDetermineBonds
from rdkit.Chem.Draw import IPythonConsole

frame_str = "{a1} {a2} {a3} BO ! atom indices involved in frame {frame} ! {frameid}"
charge_str = (
    "{nchg} 0      ! no. chgs and polarizabilities for atom {atomNumber} ({idx})"
)
blank_lines = "     0.000000        0.000000        0.000000       0.00000"


header_script = """
ROOT=/pchem-data/meuwly/boittier/home/MDCM
BINDIR=$ROOT/bin
HOME=/pchem-data/meuwly/boittier/home
INITIALXYZ=$HOME/mdcm_fast/notebooks/dc.xyz
PCUBE={CUBE}
DCUBE={CUBE}
MAXATMCHG=60
NTRY=1
QXYZ=dc.xyz
FRAMES=frames.txt
MODEL=model.mdcm
INITIALXYZ=dc.xyz
OUTPUT={NCHG}_charges_refined.xyz
OUTPUTCUBE={NCHG}charges.cube
NAME=$INITIALXYZ
"""


def mol_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            is_aromatic=atom.GetIsAromatic(),
            atom_symbol=atom.GetSymbol(),
        )

    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType()
        )

    return G


def draw_mol_nx(Gmol, labels=None, ax=None):
    Gmol_atom = nx.get_node_attributes(Gmol, "atom_symbol")
    Gmol_atom = {k: f"{k}.{v}" for k, v in list(Gmol_atom.items())}
    color_map = {"C": "orange", "O": "red", "N": "blue", "S": "yellow", "H": "grey"}

    _colors = []
    for idx in Gmol.nodes():
        if Gmol.nodes[idx]["atom_symbol"] in color_map:
            _colors.append(color_map[Gmol.nodes[idx]["atom_symbol"]])

    return nx.draw(
        Gmol,
        pos=nx.kamada_kawai_layout(Gmol),
        labels=Gmol_atom if labels is None else labels,
        with_labels=True,
        node_color=_colors,
        node_size=800,
        ax=ax,
    )


def map_fp_bits_to_atoms(mol, fp):
    bit_info = {}
    AllChem.GetMorganFingerprint(mol, radius=1, bitInfo=bit_info)
    return {bit: list(atoms) for bit, atoms in bit_info.items()}


def fingerprints(molecule):
    AllChem.Compute2DCoords(molecule)
    Chem.SanitizeMol(molecule)
    # molecule = Chem.RemoveHs(molecule)
    fp = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=1)
    atom_maps = map_fp_bits_to_atoms(molecule, fp)
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=0, bitInfo=bi)
    bits = list(bi.keys())
    pngs = []
    for bit in bits:
        hits = bi[bit]
        hit_atoms = set()
        hit_bonds = set()
        for atom_idx, radius in hits:
            if radius == 0:
                hit_atoms.add(atom_idx)
            elif radius == 1:
                bond = molecule.GetBondBetweenAtoms(atom_idx, hits[0][0])
                if bond:
                    hit_bonds.add(bond.GetIdx())
        # Draw the molecule with highlighted atoms and bonds
        d = Draw.MolDraw2DCairo(
            300, 300
        )  # Replace with appropriate drawer for your environment
        d.drawOptions().prepareMolsBeforeDrawing = False
        Draw.PrepareAndDrawMolecule(
            d, molecule, highlightAtoms=hit_atoms, highlightBonds=hit_bonds
        )
        d.FinishDrawing()
        from IPython.display import Image

        png = d.GetDrawingText()
        pngs.append(png)

    return pngs, bits, bi


def make_bit_dict(bi):
    bi_dict = {}
    for k, b in bi.items():
        for t in b:
            bi_dict[t[0]] = k
    return bi_dict


# Function to find all paths of length 3 (three nodes)
def find_paths_of_length_three(G, no_Hs=False):
    paths = []
    atomnames = nx.get_node_attributes(G, "atom_symbol")
    for source in G.nodes():
        for target in G.nodes():
            if source != target:
                # if no_Hs: and atomnames[target] != "H" and atomnames[source] != "H":
                for path in nx.all_simple_paths(
                    G, source=source, target=target, cutoff=2
                ):
                    if len(path) == 3:
                        paths.append(path)
    return paths


def make_frames(graphMol, triplets, frametypes, frame_priority):
    frames = []
    frame_types = []
    atomnames = nx.get_node_attributes(graphMol, "atom_symbol")
    heavy_nodes = [_ for _ in range(len(atomnames)) if atomnames[_] != "H"]
    hydrogens = [_ for _ in range(len(atomnames)) if atomnames[_] == "H"]
    while len(heavy_nodes) != 0:
        #  add frames by frame priority
        for pf in frame_priority:
            for i, atypframe in enumerate(frametypes):
                if atypframe == pf:
                    frames.append(triplets[i])
                    frame_types.append(pf)
                    # remove them from the nodes list
                    for a in triplets[i]:
                        if a in heavy_nodes:
                            heavy_nodes.remove(a)
    return frames, frame_types, hydrogens


def make_frame_file(resName, frames, frame_types, hydrogens):
    print("1 0          ! no. residue types defined here")
    print("")
    print(f"{resName} !residue name")
    print(f"{len(frames)} ! no. axis system frames")
    count_charges = 0
    atoms_added = []
    for i, (frame, frametype) in enumerate(zip(frames, frame_types)):
        if not np.all([a in atoms_added for a in frame]):
            print(
                frame_str.format(
                    a1=frame[0] + 1,
                    a2=frame[1] + 1,
                    a3=frame[2] + 1,
                    frame=i + 1,
                    frameid=frametype,
                )
            )
            for j, a in enumerate(frame):
                if a not in atoms_added:
                    print(charge_str.format(nchg=2, atomNumber=a + 1, idx=j))
                    print(blank_lines)
                    print(blank_lines)
                    atoms_added.append(a)
                    count_charges += 2
                else:
                    print(charge_str.format(nchg=0, atomNumber=a + 1, idx=j))
    frameCount = i
    for h in hydrogens:
        print(
            frame_str.format(
                a1=h + 1, a2=0, a3=0, frame=frameCount + 1, frameid="(hydrogen)"
            )
        )
        print(charge_str.format(nchg=1, atomNumber=h + 1, idx=0))
        print(blank_lines)
        atoms_added.append(h)
        frameCount += i
        count_charges += 1

    print("!")
    print("!ncghs", count_charges)
    atoms_added.sort()
    print("!", atoms_added)
    print("!natoms", len(atoms_added))


def batch_to_xyzblock(batch):
    block = ""
    xyz_str = ""
    count = 0
    for elem, xyz in zip(batch["Z"], batch["R"]):
        if elem != 0:
            count += 1
            xyz_str += "{} {}\n".format(
                ase.data.chemical_symbols[elem], " ".join([f"{_:.3f}" for _ in xyz])
            )
    block += f"{count}\n\n"
    block += xyz_str
    return block


def batch_to_mol(batch):
    xyz_block = batch_to_xyzblock(batch)

    mol = Chem.MolFromXYZBlock(xyz_block)
    rdDetermineBonds.DetermineBonds(mol, charge=0)
    # Annotate each atom with its index
    for atom in mol.GetAtoms():
        atom.SetProp("atomNote", str(atom.GetIdx() + 1))

    # Draw the molecule with annotated atom numbers
    # img = Draw.MolToImage(mol, kekulize=True)
    # img.show()

    return mol


def make_frames_txt(mol, resid):
    mols = [mol]
    Gmols = [mol_to_nx(mol) for mol in mols]

    pictures = []
    all_bits = []
    bi_dicts = []
    for mol in mols:
        pngs, bits, bi = fingerprints(mol)
        pictures.append(pngs)
        all_bits.append(bits)
        bi_dicts.append(make_bit_dict(bi))

    print(bi_dicts)
    names = [resid]
    formal_charges = [rdkit.Chem.rdmolops.GetFormalCharge(m) for m in mols]
    legends = [f"{n} ({fc})" for fc, n in zip(formal_charges, names)]

    paths3 = []
    results = []

    for i, G in enumerate(Gmols):
        # Find paths of length 3 in our sample graph
        paths_length_three = find_paths_of_length_three(Gmols[i])
        res = [tuple([bi_dicts[i][x] for x in _]) for _ in paths_length_three]
        paths3.append(paths_length_three)
        results.append(res)

    #  find all possible frame combinations and sort by frequency
    all_frame_types = pd.DataFrame(results).to_numpy().flatten()
    unique_frames = pd.DataFrame(all_frame_types).value_counts()
    unique_frames = list(unique_frames.index)
    frame_priority = [unique_frames[0][0]]
    for frame in unique_frames[1:]:
        #  check if the palindrome is not there
        if frame[0][-1::-1] not in frame_priority:
            frame_priority.append(frame[0])

    i = 0
    frames, frame_types, hydrogens = make_frames(
        Gmols[i], paths3[i], results[i], frame_priority
    )
    resName = legends[i].split()[0]

    frame_txt = f"{resid}\n"
    mol = mols[0]
    # bit dictionary reversed
    rev_bit_dict = {v: k for k, v in bi_dicts[i].items()}
    print(frame_types)
    # calculate the atomic weight of the frame
    weights = [
        np.sum([mol.GetAtomWithIdx(rev_bit_dict[a]).GetMass() for a in frame])
        for frame in frame_types
    ]
    print(weights)
    # sort frames by weight
    print(frames)
    frames = [x for _, x in sorted(zip(weights, frames), key=lambda pair: -pair[0])]
    print(frames)

    visited = []
    for _ in frames:
        if (_[0] not in visited) or (_[1] not in visited) or (_[2] not in visited):
            s = f"{_[0]+1} {_[1]+1} {_[2]+1} BO\n"
            frame_txt += s
            for x in _:
                if x not in visited:
                    visited.append(x)

    return frame_txt


def make_files_for_mdcm(n_dcm, mono, dipo, batch, path, cube):
    vis_out = plot_3d_models(mono, dipo, np.array([mono.flatten()]), batch, 1)
    V1, V2, V3, mol, dcmol, combined = vis_out
    ase.io.write(path / "atoms.xyz", mol)
    n_charges = n_dcm
    qs = mono.flatten()[:n_charges]
    qrs = dipo.reshape(-1, 3)[:n_charges]

    header = f"""{n_charges}
    s                      x[A]                      y[A]                      z[A]                      q[e]
    """

    file_text = ""
    file_text += header

    for xyz, q in zip(qrs, qs):
        x, y, z = xyz
        if q > 0:
            t = "O"
        else:
            t = "N"
        s = f"{t}\t{x:3.10f}\t{y:3.10f}\t{z:3.10f}\t{q:3.10f}\n"
        file_text += s

    with open(path / "dc.xyz", "w") as f:
        f.writelines(file_text)

    mol = batch_to_mol(batch)
    frames = make_frames_txt(mol, "ACE")

    with open(path / "frames.txt", "w") as f:
        f.writelines(frames)

    refine_script = """
    $BINDIR/pcubefit.x -xyz $INITIALXYZ $INITIALXYZ -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v
    """
    model_script = """$BINDIR/comb-xyz-to-dcm.pl $QXYZ $PCUBE $FRAMES $MODEL"""
    analysis_script = """$BINDIR/pcubefit.x -generate -xyz $INITIALXYZ -esp $PCUBE -dens $DCUBE -v > ${NAME}.out
# Examine quality of fitted charges by comparing newly fitted model and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $OUTPUTCUBE -dens $DCUBE > analyze-cube-${NAME}.log
    """

    # refine file
    with open(path / "refine.sh", "w") as f:
        f.writelines(header_script.format(CUBE=cube, QXYZ="atoms.xyz", NCHG=n_charges))
        f.writelines(refine_script)
    # mdcm model script
    with open(path / "model.sh", "w") as f:
        f.writelines(header_script.format(CUBE=cube, QXYZ="atoms.xyz", NCHG=n_charges))
        f.writelines(model_script)
    # analysis script
    with open(path / "analyze.sh", "w") as f:
        f.writelines(header_script.format(CUBE=cube, QXYZ="atoms.xyz", NCHG=n_charges))
        f.writelines(analysis_script)

    return mol
