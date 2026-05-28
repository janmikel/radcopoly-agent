from rdkit import Chem


def canonicalize_smiles(smiles: str) -> str:
    """
    Convert SMILES to RDKit canonical SMILES.
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    return Chem.MolToSmiles(mol)


if __name__ == "__main__":
    test = "C=CC1=CC=CC=C1"

    print("Input:")
    print(test)

    print("\nCanonical:")
    print(canonicalize_smiles(test))