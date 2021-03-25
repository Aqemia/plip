from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction
import logging


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex

def test_bayer_reference():
    # tests detection of aromatic LH on Bayer reference compound
    complex_ = characterize_complex('aqemia/Pose_Compose_Ref.pdb', 'LIG:A:1')
    interactions = complex_.interaction_sets['LIG:A:1']
    assert(len(interactions.all_hbonds_ldon)==3)
    assert(len(interactions.all_hbonds_pdon)==2)
    residues_ldon = set([hbond.resnr for hbond in interactions.all_hbonds_ldon])
    residues_pdon = set([hbond.resnr for hbond in interactions.all_hbonds_pdon])
    assert(residues_ldon== {43,79,81})
    assert(residues_pdon=={141, 81})
    print('Bayer test for aromatic LH passed !')