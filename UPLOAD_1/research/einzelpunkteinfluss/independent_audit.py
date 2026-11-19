"""Audit with a second method and intentional input corruption.

SciPy's linear programming solver checks the 34 monotone minimax problems.
The scientific input/output script itself remains standard-library only.
"""
from pathlib import Path
from fractions import Fraction as Q
import csv,hashlib,json,platform,shutil,tempfile
import numpy as np
import scipy
from scipy.optimize import linprog
import influence
import units_check

HERE=Path(__file__).resolve().parent

def lp_distance(values):
    y=np.array([float(v) for v in values]);n=len(y)
    objective=np.r_[np.zeros(n),1.]
    A=[];b=[]
    for i in range(n):
        a=np.zeros(n+1);a[i]=1.;a[-1]=-1.;A.append(a);b.append(y[i])
        a=np.zeros(n+1);a[i]=-1.;a[-1]=-1.;A.append(a);b.append(-y[i])
    for i in range(n-1):
        a=np.zeros(n+1);a[i]=1.;a[i+1]=-1.;A.append(a);b.append(0.)
    sol=linprog(objective,A_ub=A,b_ub=b,bounds=[(0,None)]*(n+1),method='highs')
    if not sol.success:raise RuntimeError(sol.message)
    return sol.fun

def rejected(mutator):
    with tempfile.TemporaryDirectory() as temp:
        dest=Path(temp)
        for name in ['nominal_data.csv','Datenherkunft.json']:
            shutil.copyfile(HERE/name,dest/name)
        shutil.copytree(HERE/'data',dest/'data')
        mutator(dest)
        old_here,old_data=influence.HERE,influence.DATA
        try:
            influence.HERE,influence.DATA=dest,dest/'nominal_data.csv'
            try:influence.load_validated_rows()
            except ValueError:return True
            return False
        finally:influence.HERE,influence.DATA=old_here,old_data

def rehash(dest):
    p=dest/'Datenherkunft.json';meta=json.loads(p.read_text())
    meta['nominal_data_sha256']=hashlib.sha256((dest/'nominal_data.csv').read_bytes()).hexdigest()
    p.write_text(json.dumps(meta))

def change_rows(dest,func):
    p=dest/'nominal_data.csv'
    with p.open(newline='') as f:rows=list(csv.DictReader(f));fields=list(rows[0])
    func(rows)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def main():
    rows=influence.load_validated_rows();comparisons=[]
    for curve in influence.EXPECTED:
        y=[Q(r['value']) for r in rows if r['curve']==curve]
        for omit in [None]+list(range(len(y))):
            values=y if omit is None else y[:omit]+y[omit+1:]
            exact=influence.floor_and_witness(values)
            lp=lp_distance(values)
            delta=abs(lp-float(exact))
            if delta>1e-10:raise RuntimeError('LP comparison failed')
            comparisons.append(dict(curve=curve,omitted_index=omit,exact_floor=str(exact),lp_floor=lp,absolute_difference=delta))
    checks={}
    checks['unchanged_input_accepted']=True
    checks['changed_value_without_new_hash_rejected']=rejected(lambda p:change_rows(p,lambda rr:rr[0].update(value='0.2')))
    def wrong_order(p):
        change_rows(p,lambda rr:rr.__setitem__(slice(0,2),[rr[1],rr[0]]));rehash(p)
    checks['wrong_time_order_even_with_changed_hash_rejected']=rejected(wrong_order)
    def wrong_value(p):
        change_rows(p,lambda rr:rr[0].update(value='0.2'));rehash(p)
    checks['wrong_value_even_with_changed_hash_rejected_by_source_mapping']=rejected(wrong_value)
    def wrong_unit(p):
        change_rows(p,lambda rr:rr[0].update(value_unit='mg/g'));rehash(p)
    checks['mixed_unit_rejected']=rejected(wrong_unit)
    def wrong_curve_unit(p):
        def mutate(rr):
            for row in rr:
                if row['curve']=='Shintani_A':row['value_unit']='mg/g'
        change_rows(p,mutate);rehash(p)
    checks['consistently_wrong_curve_unit_rejected']=rejected(wrong_curve_unit)
    def damage_source(p):
        f=p/'data/original/Shintani_A.csv';f.write_bytes(f.read_bytes()+b'\n')
    checks['changed_original_file_rejected']=rejected(damage_source)
    checks['density_basis_example_verified']=units_check.calculate()['recovered_g_per_L']=='100'
    if not all(checks.values()):raise RuntimeError(checks)
    result=dict(date='2026-11-06',lp_comparisons=len(comparisons),max_lp_difference=max(c['absolute_difference'] for c in comparisons),input_and_unit_checks=checks,checks_passed=sum(checks.values()),environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),comparisons=comparisons)
    (HERE/'control_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='comparisons'}))

if __name__=='__main__':main()
