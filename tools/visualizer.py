import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

class ChemVisualizer:
    def __init__(self):
        self.name = "ChemVisualizer"
        self.description = "Visualizes chemical molecules and reactions from SMILES strings"
        
    def detect_input_type(self, smiles_input):
        """
        Detect if the input is a reaction SMILES or molecule SMILES.
        
        Args:
            smiles_input (str): Input SMILES string
            
        Returns:
            str: 'reaction' if reaction SMILES, 'molecule' if molecule SMILES
        """
        if '>>' in smiles_input:
            return 'reaction'
        else:
            return 'molecule'

    def visualize_reaction(self, rxn_smiles, output_file=None):
        """
        Visualize a chemical reaction from a reaction SMILES string and save it to a file.
        
        Args:
            rxn_smiles (str): Reaction SMILES string (e.g., "CC(=O)O.[OH]CC>>CC(=O)OCC.O")
            output_file (str): Path to save the output image
            
        Returns:
            str: Path to the generated image file, or error message
        """
        try:
            # Create temp file if output_file is not provided
            if output_file is None:
                temp_dir = tempfile.gettempdir()
                output_file = os.path.join(temp_dir, "reaction_viz.png")
            
            # Parse the reaction SMILES
            rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
            if rxn is None:
                return f"Failed to parse reaction SMILES: {rxn_smiles}"
            
            # Check if reaction is valid
            if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
                return f"Invalid reaction - missing reactants or products: {rxn_smiles}"
            
            # Set image size based on reaction complexity
            num_components = rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates()
            width = max(800, 200 * num_components)
            height = 300
            
            # Compute 2D coordinates for the reaction
            AllChem.Compute2DCoordsForReaction(rxn)
            
            # Create a drawing object for the reaction
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            
            # Write the image to a file
            with open(output_file, 'wb') as f:
                f.write(drawer.GetDrawingText())
            
            return output_file
        
        except Exception as e:
            return f"Error visualizing reaction: {str(e)}"

    def visualize_molecule(self, smiles, output_file=None):
        """
        Visualize a molecule from a SMILES string and save it to a file.
        
        Args:
            smiles (str): SMILES string (e.g., "CC(=O)OC1=CC=CC=C1C(=O)O" for aspirin)
            output_file (str): Path to save the output image
            
        Returns:
            str: Path to the generated image file, or error message
        """
        try:
            # Create temp file if output_file is not provided
            if output_file is None:
                temp_dir = tempfile.gettempdir()
                output_file = os.path.join(temp_dir, "molecule_viz.png")
            
            # Parse the SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return f"Failed to parse SMILES: {smiles}"
            
            # Compute 2D coordinates if not already present
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)
            
            # Set image size
            width = 400
            height = 300
            
            # Create a drawing object and draw the molecule
            drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            
            # Write the image to a file
            with open(output_file, 'wb') as f:
                f.write(drawer.GetDrawingText())
            
            return output_file
        
        except Exception as e:
            return f"Error visualizing molecule: {str(e)}"

    def _run(self, smiles_input):
        """
        Visualize a SMILES string (either molecule or reaction) and save to a file.
        
        Args:
            smiles_input (str): SMILES string (either molecule or reaction)
            
        Returns:
            str: Path to the generated image file, or error message
        """
        # Detect if input is a reaction or molecule
        input_type = self.detect_input_type(smiles_input)
        
        # Visualize based on input type
        if input_type == 'reaction':
            return self.visualize_reaction(smiles_input)
        else:
            return self.visualize_molecule(smiles_input)