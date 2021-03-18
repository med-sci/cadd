from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import Chem


class PAINS:
    def __init__(self, smiles):
        self._smiles = smiles
        self._pains = []

    def get_pains(self):
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        matches = []
        pains = []
        for s in self._smiles:
            mol = Chem.MolFromSmiles(s)
            entry = catalog.GetFirstMatch(mol)
            if entry is not None:
                matches.append({
                    'smiles': s,
                    'pains': entry.GetDescription()
                })
                pains.append(s)
        self._pains = pains
        if len(matches) == 0:
            return 'No pains found'
        else:
            return matches

    def exclude(self):
        self._smiles = [s for s in self._smiles if s not in self._pains]
        return self._smiles
