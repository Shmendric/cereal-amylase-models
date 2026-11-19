"""Independent mass-balance example; not an experimental paper replication."""
from fractions import Fraction as Q
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def calculate():
    original_g_per_L=Q(100)
    aliquot_L=Q(1,1000)
    final_L=Q(1,10)
    density_kg_per_L=Q(11,10)
    solute_g=original_g_per_L*aliquot_L
    diluted_g_per_L=solute_g/final_L
    recovered_g_per_L=diluted_g_per_L*final_L/aliquot_L
    sample_kg=density_kg_per_L*aliquot_L
    recovered_g_per_kg=solute_g/sample_kg
    require(recovered_g_per_L==Q(100),'Volume balance failed')
    require(recovered_g_per_kg==Q(1000,11),'Mass concentration calculation failed')
    require(recovered_g_per_L==recovered_g_per_kg*density_kg_per_L,'Concentration basis conversion failed')
    example=dict(original_g_per_L=str(original_g_per_L),diluted_g_per_L=str(diluted_g_per_L),recovered_g_per_L=str(recovered_g_per_L),density_kg_per_L=str(density_kg_per_L),recovered_g_per_kg_exact=str(recovered_g_per_kg),same_numeric_density_factor_if_mislabelled_g_per_L=str(recovered_g_per_L/density_kg_per_L),relative_error_if_g_per_L=str(-Q(1,11)),conditions=['homogeneous aliquot','exact transferred aliquot volume','measured final volume','solute mass conserved'],scope='Dimensional counterexample to an unconditional additional density multiplier; not a reproduction of Fasheun et al.')
    (HERE/'units_check.json').write_text(json.dumps(example,indent=2)+'\n')
    return example

def require(condition,message):
    if not condition:raise ValueError(message)

if __name__=='__main__':print(json.dumps(calculate()))
