import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, AllChem, Crippen
from pathlib import Path
from typing import Dict, Optional, List, Union
import logging

# Configuración del sistema de registro
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CalculadoraADME:
    """Calculadora de propiedades ADME (Absorción, Distribución, Metabolismo y Excreción)."""
    
    NOMBRES_PROPIEDADES = [
        'Peso Molecular (Da)', 'Fórmula Molecular', 'LogP', 'Donantes de Enlaces-H',
        'Aceptores de Enlaces-H', 'TPSA', 'Enlaces Rotables', 'QED', 
        'Fracción Carbonos SP3', 'Violaciones Lipinski', 'Área de Superficie Polar',
        'Número de Anillos', 'Número de Anillos Aromáticos', 'Absorción GI', 'Permeable a BBB',
        'Violaciones Ghose', 'Violaciones Veber', 'Violaciones Egan', 'Violaciones Muegge'
    ]

    def __init__(self):
        self.contador_smiles_invalidos = 0
        self.contador_procesados = 0

    def calcular_propiedades_adme(self, smiles: str) -> Optional[Dict[str, Union[str, float, bool]]]:
        """
        Calcula propiedades ADME para una cadena SMILES dada.
        
        Args:
            smiles: Cadena SMILES que representa una molécula
            
        Returns:
            Diccionario de propiedades calculadas o None si falla el procesamiento
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.contador_smiles_invalidos += 1
                logger.warning(f"SMILES inválido: {smiles}")
                return {prop: 'SMILES inválido' for prop in self.NOMBRES_PROPIEDADES}

            # Calcular propiedades básicas
            peso_mol = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            don_enlaces_h = Descriptors.NumHDonors(mol)
            acep_enlaces_h = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            enlaces_rotables = Descriptors.NumRotatableBonds(mol)

            propiedades = {
                'SMILES': smiles,
                'Peso Molecular (Da)': round(peso_mol, 2),
                'Fórmula Molecular': Chem.rdMolDescriptors.CalcMolFormula(mol),
                'LogP': round(logp, 2),
                'Donantes de Enlaces-H': don_enlaces_h,
                'Aceptores de Enlaces-H': acep_enlaces_h,
                'TPSA': round(tpsa, 2),
                'Enlaces Rotables': enlaces_rotables,
                'QED': round(QED.qed(mol), 3),
                'Fracción Carbonos SP3': round(self._calcular_fsp3(mol), 3),
                'Violaciones Lipinski': self._calcular_violaciones_lipinski(mol),
                'Área de Superficie Polar': round(tpsa, 2),
                'Número de Anillos': Chem.rdMolDescriptors.CalcNumRings(mol),
                'Número de Anillos Aromáticos': Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
                
                # Propiedades adicionales
                'Absorción GI': self._predecir_absorcion_gi(mol),
                'Permeable a BBB': self._predecir_permeabilidad_bbb(mol),
                'Violaciones Ghose': self._calcular_violaciones_ghose(mol),
                'Violaciones Veber': self._calcular_violaciones_veber(mol),
                'Violaciones Egan': self._calcular_violaciones_egan(mol),
                'Violaciones Muegge': self._calcular_violaciones_muegge(mol)
            }
            
            self.contador_procesados += 1
            return propiedades

        except Exception as e:
            logger.error(f"Error al procesar SMILES {smiles}: {str(e)}")
            return None

    @staticmethod
    def _predecir_absorcion_gi(mol: Chem.Mol) -> str:
        """Predice la absorción gastrointestinal basada en reglas de Veber y TPSA."""
        tpsa = Descriptors.TPSA(mol)
        enlaces_rotables = Descriptors.NumRotatableBonds(mol)
        peso_mol = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        if (tpsa <= 140 and enlaces_rotables <= 10 and 
            peso_mol <= 500 and -0.4 <= logp <= 5.6):
            return "Alta"
        else:
            return "Baja"

    @staticmethod
    def _predecir_permeabilidad_bbb(mol: Chem.Mol) -> str:
        """Predice la permeabilidad de la barrera hematoencefálica basada en propiedades fisicoquímicas."""
        tpsa = Descriptors.TPSA(mol)
        peso_mol = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        don_enlaces_h = Descriptors.NumHDonors(mol)
        
        if (tpsa < 90 and peso_mol < 400 and 
            logp < 5 and don_enlaces_h <= 3):
            return "Sí"
        else:
            return "No"

    @staticmethod
    def _calcular_violaciones_ghose(mol: Chem.Mol) -> int:
        """Calcula violaciones del filtro de Ghose."""
        peso_mol = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        num_atomos = mol.GetNumAtoms()
        refractividad_molar = Crippen.MolMR(mol)
        
        violaciones = sum([
            not (160 <= peso_mol <= 480),
            not (-0.4 <= logp <= 5.6),
            not (20 <= num_atomos <= 70),
            not (40 <= refractividad_molar <= 130)
        ])
        return violaciones

    @staticmethod
    def _calcular_violaciones_veber(mol: Chem.Mol) -> int:
        """Calcula violaciones del filtro de Veber."""
        enlaces_rotables = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        
        violaciones = sum([
            enlaces_rotables > 10,
            tpsa > 140
        ])
        return violaciones

    @staticmethod
    def _calcular_violaciones_egan(mol: Chem.Mol) -> int:
        """Calcula violaciones del filtro de Egan."""
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)
        
        violaciones = sum([
            tpsa > 132,
            not (-1 <= logp <= 6)
        ])
        return violaciones

    @staticmethod
    def _calcular_violaciones_muegge(mol: Chem.Mol) -> int:
        """Calcula violaciones del filtro de Muegge."""
        peso_mol = Descriptors.ExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        anillos = Chem.rdMolDescriptors.CalcNumRings(mol)
        atomos_carbono = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        heteroatomos = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        enlaces_rotables = Descriptors.NumRotatableBonds(mol)
        
        violaciones = sum([
            not (200 <= peso_mol <= 600),
            not (-2 <= logp <= 5),
            tpsa < 75,
            anillos < 1,
            atomos_carbono < 7,
            heteroatomos < 2,
            enlaces_rotables > 15
        ])
        return violaciones

    @staticmethod
    def _calcular_fsp3(mol: Chem.Mol) -> float:
        """Calcula la fracción de carbonos sp3."""
        if mol is None:
            return 0.0
        num_sp3 = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetHybridization() == Chem.HybridizationType.SP3)
        num_carbono = sum(1 for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() == 6)
        return num_sp3 / num_carbono if num_carbono > 0 else 0.0

    @staticmethod
    def _calcular_violaciones_lipinski(mol: Chem.Mol) -> int:
        """Calcula violaciones de la Regla de los Cinco de Lipinski."""
        violaciones = sum([
            Descriptors.ExactMolWt(mol) > 500,
            Descriptors.MolLogP(mol) > 5,
            Descriptors.NumHDonors(mol) > 5,
            Descriptors.NumHAcceptors(mol) > 10
        ])
        return violaciones

    def procesar_archivo(self, archivo_entrada: Union[str, Path], archivo_salida: Union[str, Path], 
                    columna_smiles: str = 'SMILES') -> None:
        """
        Procesa propiedades ADME para moléculas desde un archivo de entrada.
        
        Args:
            archivo_entrada: Ruta al archivo Excel de entrada
            archivo_salida: Ruta para guardar el archivo Excel de salida
            columna_smiles: Nombre de la columna que contiene las cadenas SMILES
        """
        try:
            ruta_entrada = Path(archivo_entrada)
            if not ruta_entrada.exists():
                raise FileNotFoundError(f"Archivo de entrada no encontrado: {archivo_entrada}")

            df = pd.read_excel(ruta_entrada)
            if columna_smiles not in df.columns:
                raise ValueError(f"Columna '{columna_smiles}' no encontrada en el archivo de entrada")

            resultados = []
            total_moleculas = len(df)
            
            for idx, fila in df.iterrows():
                if idx % 100 == 0:
                    logger.info(f"Procesando molécula {idx + 1}/{total_moleculas}")
                
                smiles = fila[columna_smiles]
                if pd.isna(smiles):
                    continue
                    
                propiedades = self.calcular_propiedades_adme(str(smiles))
                if propiedades:
                    fila_combinada = {**fila.to_dict(), **propiedades}
                    resultados.append(fila_combinada)

            df_resultados = pd.DataFrame(resultados)
            df_resultados.to_excel(archivo_salida, index=False)
            
            logger.info(f"Procesamiento completado:")
            logger.info(f"- Total de moléculas procesadas: {self.contador_procesados}")
            logger.info(f"- SMILES inválidos encontrados: {self.contador_smiles_invalidos}")
            logger.info(f"- Resultados guardados en: {archivo_salida}")

        except Exception as e:
            logger.error(f"Error al procesar el archivo: {str(e)}")
            raise

def main():
    archivo_entrada = 'Moleculas.xlsx'
    archivo_salida = 'propiedades_adme.xlsx'
    
    calculadora = CalculadoraADME()
    calculadora.procesar_archivo(archivo_entrada, archivo_salida)

if __name__ == "__main__":
    main()