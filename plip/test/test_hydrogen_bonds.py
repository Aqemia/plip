import unittest

from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction


def characterize_complex(pdb_file: str, binding_site_id: str, as_string: bool = False) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file, as_string=as_string)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex.interaction_sets[binding_site_id]


def strip_hydrogens(pdb_file: str) -> str:
    """Return the ATOM/HETATM records of a PDB file without hydrogens (and without CONECT records)."""
    with open(pdb_file) as handle:
        return ''.join(line for line in handle
                       if line.startswith(('ATOM', 'HETATM', 'TER', 'END')) and line[76:78] != ' H')


class HydrogenBondTestCase(unittest.TestCase):

    def setUp(self):
        self._nohydro = config.NOHYDRO

    def tearDown(self):
        config.NOHYDRO = self._nohydro

    def test_4dst_nondeterministic_protonation(self):
        config.NOHYDRO = False
        for i in range(0, 10):
            interactions = characterize_complex('./pdb/4dst.pdb', 'GCP:A:202')
            all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
            self.assertTrue(len(all_hbonds) == 16 or len(all_hbonds) == 17)

    def test_4dst_deterministic_protonation(self):
        config.NOHYDRO = True
        for i in range(0, 10):
            interactions = characterize_complex('./pdb/4dst_protonated.pdb', 'GCP:A:202')
            all_hbonds = interactions.hbonds_ldon + interactions.hbonds_pdon
            self.assertTrue(len(all_hbonds) == 16)

    def test_no_protonation(self):
        # Input without hydrogens: ligand donors are only found if PLIP protonates
        pdb_string = strip_hydrogens('./pdb/1x0n_state_1.pdb')
        config.NOHYDRO = True
        interactions1 = characterize_complex(pdb_string, 'DTF:A:174', as_string=True)
        self.assertEqual(len(interactions1.hbonds_ldon), 0)
        config.NOHYDRO = False
        interactions2 = characterize_complex(pdb_string, 'DTF:A:174', as_string=True)
        self.assertGreater(len(interactions2.hbonds_ldon), 0)

    def test_no_protonation_explicit_hydrogens(self):
        # Hydrogens present in the input (PyMOL-style TER record without serial number) are used as-is
        config.NOHYDRO = True
        interactions = characterize_complex('./pdb/1x0n_state_1.pdb', 'DTF:A:174')
        self.assertEqual({hbond.resnr for hbond in interactions.hbonds_ldon}, {107, 109})
