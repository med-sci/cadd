from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw


UNWANTED = Path('cadd', 'ADME', 'unwanted.txt')


class UnwantedSubs:
    def __init__(self, smiles):
        self.unwanted_subs = self._unwanted
        self._smiles = smiles

    @property
    def _unwanted(self):
        unwanted = {}
        with open(UNWANTED, 'r') as f:
            for line in f.readlines():
                unwanted[line.split()[0]] = line.split()[1]
        f.close()
        return unwanted

    def get_unwanted(self):
        matches = []
        clean = []
        for smile in self._smiles:
            mol = Chem.MolFromSmiles(smile)
            clear = True
            for name, smarts in self.unwanted_subs.items():
                if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                    matches.append(
                        {
                            'smiles': smile,
                            'name': name,
                            'smarts': smarts
                        }
                    )
                    clear = False
            if clear:
                clean.append(smile)
        return matches, clean

    @staticmethod
    def visualize(u_subs):
        matches = []
        mols = []
        names = []
        for i in u_subs:
            mol = Chem.MolFromSmiles(i['smiles'])
            matches.append(mol.GetSubstructMatch(Chem.MolFromSmarts(i['smarts'])))
            names.append(i['name'])
            mols.append(mol)
        return Draw.MolsToGridImage(mols, highlightAtomLists=matches, legends=names)




