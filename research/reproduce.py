"""Reproduce all selected public calculations in an isolated file copy.

No network access, manuscript build or publication. Run without -O/-OO.
"""
from pathlib import Path
import argparse,json,shutil,subprocess,sys,tempfile,time

R=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Ohne -O oder -OO ausführen.')
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--skip-optimization',action='store_true',help='Gespeicherte E/W/P2-Parameter nachrechnen, statt neu zu optimieren.')
ap.add_argument('--output',type=Path,default=R/'ergebnisse/laufprotokolle')
a=ap.parse_args();records=[];ok=True
def run(label,args,cwd):
    start=time.monotonic()
    p=subprocess.run([sys.executable,*map(str,args)],cwd=cwd,capture_output=True,text=True)
    records.append(dict(label=label,args=list(map(str,args)),seconds=round(time.monotonic()-start,3),returncode=p.returncode,stdout=p.stdout,stderr=p.stderr))
    print(label+': '+('BESTANDEN' if p.returncode==0 else 'FEHLER'),flush=True)
    if p.returncode:print(p.stderr or p.stdout,file=sys.stderr)
    return p.returncode==0
ok=run('Integrität, Dateiauswahl und Quellenzuordnung',['verify_release.py'],R)
if ok:
    with tempfile.TemporaryDirectory(prefix='getreide_public1_') as temp:
        W=Path(temp)/'paket'
        shutil.copytree(R,W,ignore=shutil.ignore_patterns('__pycache__','laufprotokolle'))
        tasks=[
          ('Nominale Vorhersage, INFOGEST und synthetische Grundmodelle',['schlussfassung17/nachrechnen.py','--package',str(W),'--output',str(W/'ergebnisse/nominal')]),
          ('1572 gespeicherte Szenarien und Rekonstruktionsfälle',['schlussfassung17/ergaenzungen.py']),
          ('Vier rationale globale E Zertifikate',['zertifizierte_fehlergrenzen/verify_certificates.py','--check-only']),
          ('Negative Zertifikatskontrollen',['zertifizierte_fehlergrenzen/check_integrity.py']),
          ('Exakte Einzelpunktanalyse',['einzelpunkteinfluss/influence.py']),
          ('LP und Einheiten Gegenrechnung',['einzelpunkteinfluss/independent_audit.py']),
          ('Hochpräzisionsgegenrechnung und 80 INFOGEST Zahlenfelder',['schlussfassung19/gegenpruefung/audit.py']),
          ('Zusatzmessmodelle',['schlussfassung21/zusatzmessungen/analyse.py']),
          ('Getrennte Substratzugänge',['schlussfassung22/haertetest/analyse.py']),
          ('Exaktes synthetisches Kalibrierbeispiel',['schlussfassung28/kalibrierbeispiel.py'])]
        if not a.skip_optimization:tasks.append(('Nominaler begrenzter E W P2 Vergleich neu optimiert',['schlussfassung28/modellvergleich.py']))
        tasks.append(('Poollabel Rang, Fehlerkovarianz und andere LP Suchstruktur',['schlussfassung29/wissenschaft_pruefen.py']))
        tasks.extend([('Numerische Kontur- und Rundungskontrolle',['schlussfassung33/numerische_kontrolle.py']),('168 Formel- und Tabellenbedingungen',['pruefung/formeln_und_tabellen.py'])])
        for label,args in tasks:
            if not run(label,args,W):ok=False;break
        a.output.mkdir(parents=True,exist_ok=True)
        for rel in ['ergebnisse/nominal','schlussfassung17/ergaenzungen_ergebnisse.json',
            'zertifizierte_fehlergrenzen/results/integrity_checks.json','einzelpunkteinfluss/influence_summary.json',
            'einzelpunkteinfluss/control_checks.json','schlussfassung19/gegenpruefung/Ergebnisse.json',
            'schlussfassung21/zusatzmessungen/ergebnisse','schlussfassung22/haertetest/ergebnisse',
            'schlussfassung28/kalibrierpruefung.json','schlussfassung28/modellpruefung.json',
            'schlussfassung29/wissenschaft_pruefung.json','schlussfassung33/Numerische_Kontrolle.json','pruefung/Formel_und_Tabellenabnahme.json']:
            src=W/rel;dest=a.output/rel
            if src.is_dir():shutil.copytree(src,dest,dirs_exist_ok=True)
            elif src.exists():dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
        # Capture true test counts from the executed modules, without treating
        # their partially overlapping assertions as independent studies.
        summary={}
        for rel,key in [('ergebnisse/nominal/rechenpruefung.json','checks'),
            ('schlussfassung17/ergaenzungen_ergebnisse.json','checks'),
            ('schlussfassung21/zusatzmessungen/ergebnisse/pruefergebnisse.json','checks'),
            ('schlussfassung22/haertetest/ergebnisse/pruefergebnisse.json','checks'),
            ('schlussfassung29/wissenschaft_pruefung.json','checks')]:
            p=W/rel
            if p.exists():summary[rel]=len(json.loads(p.read_text())[key])
else:summary={}
report=dict(package_version='1.0',all_passed=ok,steps=records,executed_module_check_counts=summary,
    reoptimized_nominal_models=not a.skip_optimization,
    scope='Alle ausgewählten öffentlichen Rechenmodule; interne technische Prüfung. Keine vollständige Neuoptimierung der 1572 historischen Fälle, kein Word-/PDF-Neusatz, keine frische Abhängigkeitsinstallation, keine unabhängige Fachprüfung oder Laborvalidierung. Die 168 neuen Prüfbedingungen enthalten wiederholte Gitterfälle.')
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'gesamtlauf.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
sys.exit(0 if ok else 1)
