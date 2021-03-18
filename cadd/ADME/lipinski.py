from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
import pandas as pd
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt


LABELS = "Molecular_weight Acceptors Donors LogP"


class LipinskiCalc:
    def __init__(self, smiles):
        self.smiles = smiles
        self.props = self._get_lipinski_props()
        self.lipinski_table = self._props

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

    @staticmethod
    def _get_stats(df):
        means = []
        plus_stds = []
        minus_stds = []
        for column in LABELS.split():
            m = mean(df[column])
            std = stdev(df[column])
            if column == "Molecular_weight":
                m = m / 100
                std = stdev(df[column]) / 100
            elif column == "Acceptors":
                m = m / 2
                std = stdev(df[column]) / 2
            means.append(m)
            plus_stds.append(m + std)
            minus_stds.append(m - std)
        means = means + [means[0]]
        plus_stds = plus_stds + [plus_stds[0]]
        minus_stds = minus_stds + [minus_stds[0]]
        return means, plus_stds, minus_stds

    def visualize(self, df):
        line1, line2, line3 = self._get_stats(df)
        label_places = np.linspace(start=0, stop=2 * np.pi, num=len(line1))
        plt.figure(figsize=(6, 6))
        ax = plt.subplot(polar=True)
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(180)
        plt.xticks(label_places, LABELS.split() + [""], fontsize=16)
        ax.fill(label_places, [5] * 5, "cornflowerblue", alpha=0.2)
        ax.plot(label_places, line1, "b", lw=3, ls="-")
        ax.plot(label_places, line2, "orange", lw=2, ls="--")
        ax.plot(label_places, line3, "orange", lw=2, ls="-.")
        labels = ("mean", "mean + std", "mean - std", "rule of five area")
        ax.legend(labels, loc=(1.1, 0.7), labelspacing=0.3, fontsize=16)
        plt.show()
