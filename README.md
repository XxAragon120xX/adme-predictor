# ADME Properties Calculator with RDKit
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2023.09-orange)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Este proyecto calcula propiedades farmacocinéticas **ADME** (Absorción, Distribución, Metabolismo y Excreción) a partir de moléculas representadas como **SMILES** utilizando la biblioteca **RDKit**.

## Propiedades calculadas

- **LogP** (hidrofobicidad)
- **TPSA** (Área de superficie polar topológica)
- **Peso molecular**
- **Número de donadores y aceptores de puentes de hidrógeno**
- **Número de enlaces rotables**
- **Fracción de Csp³**
- **Número de anillos y anillos aromáticos**
- **Cumplimiento de la Regla de Lipinski**
- **Fórmula molecular**
- **QED** (Quantitative Estimation of Drug-likeness)
- **Predicción de absorción gastrointestinal**
- **Predicción de permeabilidad de la barrera hematoencefálica (BBB)**
- **Filtros de medicinalidad adicionales:**
  - Ghose
  - Veber
  - Egan
  - Muegge

## Estructura del proyecto

```
├── adme_calculator.py    # Script principal con las funciones ADME
├── Moleculas.xlsx        # Archivo de entrada con SMILES 
├── propiedades_adme.xlsx # Archivo de salida con propiedades calculadas
└── README.md             # Este archivo
```

## Requisitos

- Python 3.8+
- RDKit 2023.09+
- pandas
- pathlib
- typing

## Instalación

1. Clona este repositorio:
   ```bash
   git clone https://github.com/tuusuario/adme-calculator.git
   cd adme-calculator
   ```

2. Crea un entorno virtual e instala las dependencias:
   ```bash
   python -m venv venv
   source venv/bin/activate  # En Windows: venv\Scripts\activate
   pip install rdkit pandas
   ```

## Uso

### Uso básico

Ejecuta el script principal:

```bash
python adme_calculator.py
```

Por defecto, el script procesa el archivo `Moleculas.xlsx` y guarda los resultados en `propiedades_adme.xlsx`.

### Uso como módulo

También puedes importar la clase `CalculadoraADME` en tus propios scripts:

```python
from adme_calculator import CalculadoraADME

# Inicializar la calculadora
calculadora = CalculadoraADME()

# Calcular propiedades para un SMILES individual
smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirina
propiedades = calculadora.calcular_propiedades_adme(smiles)
print(propiedades)

# Procesar un archivo completo
calculadora.procesar_archivo('entrada.xlsx', 'salida.xlsx', columna_smiles='SMILES')
```

## Filtros de medicinalidad implementados

El script implementa varios filtros de medicinalidad para evaluar compuestos:

1. **Regla de los cinco de Lipinski**: 
   - Peso molecular < 500 Da
   - LogP < 5
   - Donantes de H < 5
   - Aceptores de H < 10

2. **Filtro de Ghose**:
   - 160 ≤ Peso molecular ≤ 480
   - -0.4 ≤ LogP ≤ 5.6
   - 20 ≤ Número de átomos ≤ 70
   - 40 ≤ Refractividad molar ≤ 130

3. **Filtro de Veber**:
   - Enlaces rotables ≤ 10
   - TPSA ≤ 140 Å²

4. **Filtro de Egan (absorción pasiva)**:
   - TPSA ≤ 132 Å²
   - -1 ≤ LogP ≤ 6

5. **Filtro de Muegge**:
   - 200 ≤ Peso molecular ≤ 600
   - -2 ≤ LogP ≤ 5
   - TPSA ≥ 75 Å²
   - Anillos ≥ 1
   - Átomos de carbono ≥ 7
   - Heteroátomos ≥ 2
   - Enlaces rotables ≤ 15

## Predicciones implementadas

- **Absorción gastrointestinal**: Basada en las propiedades físico-químicas y los filtros de Veber.
- **Permeabilidad de la barrera hematoencefálica**: Basada en TPSA, peso molecular, LogP y donantes de H.

## Contribuir

Las contribuciones son bienvenidas. Por favor, siéntete libre de hacer un fork, crear una rama, y enviar pull requests.

1. Haz un fork del proyecto
2. Crea una rama para tu funcionalidad (`git checkout -b feature/amazing-feature`)
3. Haz commit de tus cambios (`git commit -m 'Add amazing feature'`)
4. Haz push a la rama (`git push origin feature/amazing-feature`)
5. Abre un Pull Request

## Licencia

Este proyecto está licenciado bajo la Licencia MIT - ver el archivo [LICENSE](LICENSE) para más detalles.
