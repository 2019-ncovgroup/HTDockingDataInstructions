import glob
import sys
import numpy as np
from openeye import oedocking, oechem
from pymol import cmd


def create_receptor(output_name, pdb_file, bp_min, bp_max):
    check_oeb = output_name
    proteinStructure = oechem.OEGraphMol()
    ifs = oechem.oemolistream(pdb_file)
    ifs.SetFormat(oechem.OEFormat_PDB)
    oechem.OEReadMolecule(ifs, proteinStructure)

    box = oedocking.OEBox(*bp_max, *bp_min)

    receptor = oechem.OEGraphMol()
    s = oedocking.OEMakeReceptor(receptor, proteinStructure, box)
    oedocking.OEWriteReceptorFile(receptor, check_oeb)


def minN(x, y):
    if x is None and y is None:
        return None
    if x is None:
        return y
    if y is None:
        return x
    if x <= y:
        return x
    return y


def maxN(x, y):
    if x is None and y is None:
        return None
    if x is None:
        return y
    if y is None:
        return x
    if y >= x:
        return y
    return x


def get_min_max(fname, inc=5):
    data = {}
    with open(fname, 'r') as f:
        next(f)
        next(f)
        next(f)
        for line in f:
            line = line.split()
            if line[0] == 'TER':
                break

            pocket_id = int(line[4])
            if pocket_id in data:
                data[pocket_id].append((float(line[5]), float(line[6]), float(line[7])))
            else:
                data[pocket_id] = [(float(line[5]), float(line[6]), float(line[7]))]

    res = {}
    for pocket, data in data.items():
        x_min, x_max = None, None
        y_min, y_max = None, None
        z_min, z_max = None, None

        for coord in data:
            x_min = minN(x_min, coord[0])
            x_max = maxN(x_max, coord[0])
            y_min = minN(y_min, coord[1])
            y_max = maxN(y_max, coord[1])
            z_min = minN(z_min, coord[2])
            z_max = maxN(z_max, coord[2])

        x_min -= inc
        x_max += inc
        y_min -= inc
        y_max += inc
        z_min += inc
        z_max += inc

        res[pocket] = (x_min, x_max, y_min, y_max, z_min, z_max)

    return res


def read_fppocket_file(fname):
    pockets = []
    pocket_count = 1
    with open(fname, 'r') as f:
        for line in f:
            if "Pocket" in line:
                for line in f:
                    if "Volume" in line:
                        line = float(line.strip().split('\t')[-1])
                        pockets.append((pocket_count, line))
                        pocket_count += 1
                        break

    return pockets


directories = glob.glob(sys.argv[1] + "*/")

with open('out.csv', 'w') as fout:
    for struct in directories:
        pockets = read_fppocket_file(glob.glob(struct + "*.txt")[0])
        protein_name = struct.split("/")[-2].split("_")[0]
        pocket_data = get_min_max(glob.glob(struct + "*.pqr")[0])
        pdb_file = glob.glob(struct + "*.pdb")[0]

        cmd.reinitialize()
        cmd.load(pdb_file)
        cmd.do("remove not polymer")
        cmd.do("save {}".format(pdb_file))

        mean = np.mean(list(zip(*pockets))[1])
        std = np.std(list(zip(*pockets))[1])
        print(mean, std)
        for pocket_id, volume in pockets:
            print(pocket_id, protein_name)
            if volume > (mean + std):
                # min_coords, max_coords =
                print(pocket_id)
                print(pocket_data[pocket_id])
                create_receptor("out/" + protein_name + "_pocket" + str(pocket_id) + "_receptor.oeb", pdb_file,
                                [pocket_data[pocket_id][0], pocket_data[pocket_id][2], pocket_data[pocket_id][4]],
                                [pocket_data[pocket_id][1], pocket_data[pocket_id][3], pocket_data[pocket_id][5]])

                fout.write(
                    "{},{},{},{},{},{},{},{},{},{}\n".format(protein_name + "_pocket" + str(pocket_id), protein_name,
                                                             pocket_id, glob.glob(struct + "*.pdb")[0],
                                                             pocket_data[pocket_id][0], pocket_data[pocket_id][1],
                                                             pocket_data[pocket_id][2], pocket_data[pocket_id][3],
                                                             pocket_data[pocket_id][4], pocket_data[pocket_id][5]))
