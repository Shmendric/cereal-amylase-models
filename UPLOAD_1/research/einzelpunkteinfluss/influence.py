"""Exact single-point influence on the known monotone minimax floor.

This is a sensitivity analysis of nominal digitized data, not predictive
cross-validation, an outlier test, or a rule for deleting observations.
Only the Python standard library is needed.
"""
from pathlib import Path
from fractions import Fraction as Q
import csv, hashlib, json

HERE=Path(__file__).resolve().parent
DATA=HERE/'nominal_data.csv'
EXPECTED={'Shintani_A':8,'Shintani_F':8,'H_STD':7,'H_HIGH':7}
EXPECTED_META={
    'Shintani_A':('h','g maltose / g input flour','GRP005'),
    'Shintani_F':('h','g maltose / g input flour','GRP005'),
    'H_STD':('min','mg maltose / g dry flour','GRP010'),
    'H_HIGH':('min','mg maltose / g dry flour','GRP010'),
}

def require(condition,message):
    if not condition:raise ValueError(message)

def load_validated_rows():
    provenance=json.loads((HERE/'Datenherkunft.json').read_text(encoding='utf-8'))
    require(hashlib.sha256(DATA.read_bytes()).hexdigest()==provenance['nominal_data_sha256'],'Nominal data hash mismatch')
    for source in provenance['source_files']:
        path=HERE/source['local_path']
        require(path.is_file(),'Missing original technical file: '+source['local_path'])
        require(hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256'],'Original technical file hash mismatch')
    with DATA.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    require(set(r['curve'] for r in rows)==set(EXPECTED),'Curve identities differ')
    require(len({r['point_id'] for r in rows})==len(rows),'Duplicate point identifier')
    for curve,count in EXPECTED.items():
        rr=[r for r in rows if r['curve']==curve]
        require(len(rr)==count,'Point count differs for '+curve)
        times=[Q(r['time']) for r in rr]
        require(times[0]>=0 and all(a<b for a,b in zip(times,times[1:])),'Times are not strictly ordered')
        require(all(Q(r['value'])>=0 for r in rr),'Negative nominal value')
        require(len({(r['time_unit'],r['value_unit'],r['source_group']) for r in rr})==1,'Mixed units or origins')
        require(all((r['time_unit'],r['value_unit'],r['source_group'])==EXPECTED_META[curve] for r in rr),'Unexpected units or origin for '+curve)
    # Verify that normalized values and times are copied from the retained files.
    original={}
    for curve,field in [('Shintani_A','time_digitized_h'),('Shintani_F','time_h')]:
        with (HERE/f'data/original/{curve}.csv').open(encoding='utf-8-sig',newline='') as f:
            for row in csv.DictReader(f):
                original[(curve,row['point_id'])]=(row[field],row['maltose_g_per_g_input_flour'])
    with (HERE/'data/original/Hackenberg.csv').open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            original[(row['curve_id'],row['curve_id']+'_'+row['time_min'])]=(row['time_min'],row['maltose_mg_per_g_dry_flour'])
    actual={(r['curve'],r['point_id']):(r['time'],r['value']) for r in rows}
    require(actual==original,'Nominal mapping differs from original technical files')
    return rows

def floor_and_witness(values):
    require(len(values)>=2 and all(v>=0 for v in values),'At least two nonnegative values required')
    floor=max([Q(0)]+[(values[i]-values[j])/2 for i in range(len(values)) for j in range(i+1,len(values))])
    witness=[max(Q(0),max(values[:i+1])-floor) for i in range(len(values))]
    require(all(witness[i]<=witness[i+1] for i in range(len(values)-1)),'Witness is not monotone')
    require(max(abs(a-b) for a,b in zip(values,witness))==floor,'Witness error differs from floor')
    return floor

def dec(value):
    # Halved differences terminate; arbitrary ratios are rounded to 80 significant digits.
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec=80
        return str(Decimal(value.numerator)/Decimal(value.denominator))

def main():
    rows=load_validated_rows()
    detailed=[];summary=[]
    for curve in dict.fromkeys(r['curve'] for r in rows):
        rr=[r for r in rows if r['curve']==curve]
        values=[Q(r['value']) for r in rr]
        original=floor_and_witness(values)
        alternatives=[]
        for i,row in enumerate(rr):
            reduced=floor_and_witness(values[:i]+values[i+1:])
            require(reduced<=original,'Removing a point increased the monotone floor')
            record=dict(curve=curve,omitted_point=row['point_id'],omitted_time=row['time'],full_floor=dec(original),reduced_floor=dec(reduced),drop_in_floor=dec(original-reduced),value_unit=row['value_unit'])
            detailed.append(record);alternatives.append((reduced,row['point_id']))
        minimum=min(v for v,p in alternatives)
        summary.append(dict(curve=curve,n_points=len(rr),full_floor=dec(original),minimum_after_single_omission=dec(minimum),largest_relative_reduction=dec((original-minimum)/original) if original>0 else '0',most_influential_points=[p for v,p in alternatives if v==minimum],value_unit=rr[0]['value_unit']))
    for name,records in [('influence_all_points.csv',detailed),('influence_summary.csv',summary)]:
        with (HERE/name).open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    result=dict(date='2026-11-06',revision='Kontrollfassung 2',data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),calculation='exact Fraction arithmetic; nonterminating ratio displays rounded to 80 significant decimal digits',source='Technical data from project novelty package dated 2026-11-05',curves=summary,checks=dict(full_and_reduced_sequences=len(summary)+len(detailed),all_constructive_witnesses_verified=True,all_single_omission_floors_no_greater_than_full=True,source_hashes_and_row_mapping_verified=True,strict_time_order_and_units_verified=True))
    (HERE/'influence_summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
