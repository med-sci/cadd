from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
import pandas as pd


class ADME:
    def __init__(self, smiles):
        self.smiles = smiles
        self.props = self._get_lipinski_props()
        self.lipinski_props = self._props

    def _get_lipinski_props(self):
        mols = [Chem.MolFromSmiles(smile) for smile in self.smiles]
        props_df = pd.DataFrame(
            {
                "Smile": self.smiles,
                "Molecular_weight": [Descriptors.ExactMolWt(m) for m in mols],
                "Acceptors": [Lipinski.NumHAcceptors(m) for m in mols],
                "Donors": [Lipinski.NumHDonors(m) for m in mols],
                "LogP": [Crippen.MolLogP(m) for m in mols],
            }
        )
        return props_df

    def _check_props(self):
        fulfillment = []
        for i in range(self.props.shape[0]):
            s = sum(
                [
                    self.props["Molecular_weight"][i] <= 500,
                    self.props["Acceptors"][i] <= 10,
                    self.props["Donors"][i] <= 5,
                    self.props["LogP"][i] <= 5,
                ]
            )
            if s >= 4:
                fulfillment.append(True)
            else:
                fulfillment.append(False)
        return fulfillment

    @property
    def _props(self):
        lp = self._get_lipinski_props()
        f = self._check_props()
        lp["Fulfill"] = f
        return lp
