from chembl_webresource_client.new_client import new_client
import pandas as pd
import xml.etree.ElementTree as ET
from requests import request
from tqdm.notebook import tqdm
import numpy as np


class DataLoader:

    _BIOACTIVITY_API = new_client.activity

    def __init__(self, target, bioactivity_type, assay_type, target_organism):
        self._target = target
        self._bioactivity_type = bioactivity_type
        self._assay_type = assay_type
        self._target_organism = target_organism
        self._compounds = []
        self._smiles = []

    def _get_bioactivity(self):
        bioactivities = DataLoader._BIOACTIVITY_API.filter(
            target_chembl_id=self._target,
            type=f'{self._bioactivity_type}', assay_type=f'{self._assay_type}',
            target_organism=f'{self._target_organism}').only(
            "molecule_chembl_id",
            "type",
            "standard_units",
            "standard_value"
        )
        bioactivities_df = pd.DataFrame.from_records(bioactivities)
        bioactivities_df.dropna(axis=0, inplace=True)
        bioactivities_df.drop_duplicates(inplace=True)
        self._compounds = list(bioactivities_df['molecule_chembl_id'])
        return bioactivities_df

    def _get_smiles(self):
        for compound in tqdm(self._compounds):
            r = request(url=f'https://www.ebi.ac.uk/chembl/api/data/molecule/{compound}', method='get')
            tree = ET.fromstring(r.text)
            smile = tree.find('molecule_structures').find('canonical_smiles').text
            self._smiles.append(smile)
        return self._smiles

    def get_data(self):
        bioact = self._get_bioactivity()
        smiles = self._get_smiles()
        bioact.insert(1, 'smiles', smiles)
        return bioact


def calc_pIC50(data):
    try:
        values = [float(v) for v in list(data['standard_value'])]
        for index, value in enumerate(values):
            if value != 0:
                values[index] = round(-np.log10(value * 10 ** -9), 3)
            else:
                values[index] = 0
    except KeyError:
        print('standard_value column not found')
    return values


class SmilesConverter:
    def __init__(self, smiles):
        self.smiles = smiles
        self.alphabet = self._alphabet

    @property
    def _alphabet(self):
        alphabet = []
        for smile in self.smiles:
            for i, s in enumerate(smile):
                if smile[i] not in alphabet:
                    alphabet.append(s)
        return alphabet

    @staticmethod
    def _smile_to_array(smile, alphabet):
        height = len(alphabet)
        width = len(smile)
        array = np.zeros((height, width))
        for i in range(height):
            for j in range(width):
                if smile[j] == alphabet[i]:
                    array[i, j] = 1
        return array

    @staticmethod
    def _pad_data(data, max_width):
        for i, arr in enumerate(data):
            to_padd = max_width - arr.shape[1]
            data[i] = np.pad(arr, ((0, 0), (0, to_padd)), 'constant')
        return data

    def smile_to_one_hot(self):
        data = []
        max_width = 0
        for smile in self.smiles:
            arr = self._smile_to_array(smile, self.alphabet)
            data.append(arr)
            if arr.shape[1] > max_width:
                max_width = arr.shape[1]
        data = self._pad_data(data, max_width)
        data = np.array(data)
        return data


