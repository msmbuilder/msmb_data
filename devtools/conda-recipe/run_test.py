import sys
import glob
import re
import mdtraj

DIRNAME = "{}/share/msmb_data".format(sys.prefix)


def get_trajs(dataname):
    return sorted(
        glob.glob("{}/trajectory*.*".format(dataname)),
        key=lambda x: int(re.match(r'trajectory\-(\d+)\.\w\w\w', x).group(1))
    )


def get_top(dataname):
    topname = {
        'met_enkephalin': '1plx',
        'alanine_dipeptide': 'ala2',
        'fs_peptide': 'fs-peptide',

    }[dataname]
    return ("{home}/{dataname}/{topname}.pdb"
            .format(home=DIRNAME, dataname=dataname, topname=topname))


def get(dataname):
    return get_top(dataname), get_trajs(dataname)


def test(dataname):
    print("Testing {} ...".format(dataname))
    top = mdtraj.load(get_top(dataname))
    for tfn in get_trajs(dataname):
        mdtraj.load(tfn, top=top)


test('met_enkephalin')
test('fs_peptide')
test('alanine_dipeptide')
