"""Verify the exact publication files and the preservation of scientific content."""
from pathlib import Path
from hashlib import sha256
from zipfile import ZipFile
import argparse,csv,json,re
import numpy as np
import yaml,fitz
from docx import Document
from lxml import etree

R=Path(__file__).resolve().parent
EXCLUDED={'MANIFEST_SHA256.json','SHA256SUMS.txt'}
def included(p):
 rel=p.relative_to(R)
 return p.is_file() and rel.as_posix() not in EXCLUDED and '__pycache__' not in rel.parts and not rel.as_posix().startswith('ergebnisse/laufprotokolle/')
def sha(p):return sha256(p.read_bytes()).hexdigest()
def snapshot():
 paths=list(R.rglob('*'))
 if any(p.is_symlink() for p in paths):raise ValueError('Symlinks sind nicht Bestandteil des Releases.')
 return {p.relative_to(R).as_posix():sha(p) for p in sorted(paths) if included(p)}
def verify():
 expected=json.loads((R/'MANIFEST_SHA256.json').read_text());actual=snapshot()
 if actual!=expected:raise ValueError(json.dumps(dict(missing=sorted(expected.keys()-actual.keys()),added=sorted(actual.keys()-expected.keys()),changed=sorted(k for k in actual.keys()&expected.keys() if actual[k]!=expected[k]))))
 if (R/'SHA256SUMS.txt').read_text()!=''.join(v+'  '+k+'\n' for k,v in expected.items()):raise ValueError('Manifestdarstellungen weichen ab.')
 checks={}
 def check(name,ok):
  if name in checks:raise ValueError('Doppelter Prüfbereich')
  checks[name]=bool(ok)
  if not ok:raise AssertionError(name)
 read=lambda name:json.loads((R/name).read_text())
 check('Versanddateien und beide Manifestdarstellungen stimmen',True)
 rows=list(csv.DictReader((R/'DATEIRECHTE.csv').open(encoding='utf-8-sig')))
 check('Jede Nutzdatei genau einmal zugeordnet',len(rows)==len({x['path'] for x in rows}) and {x['path'] for x in rows}==set(expected))
 check('Dateirechte an Dateiprüfsummen gebunden',all(x['sha256']==actual[x['path']] if x['path']!='DATEIRECHTE.csv' else x['sha256']=='siehe MANIFEST_SHA256.json' for x in rows))
 check('Jede Datei mit Herkunft und Nutzungsbereich',all(x['source_id'] and x['source_reference'] and x['rights_scope'] and x['license_scope'] and x['changes'] for x in rows))
 cs=read('monografie/manuskript.json');bs=[b for c in cs for b in c['blocks']]
 figures={b['image'] for b in bs if 'image' in b}
 allowed_binaries={'Monografie_Getreideamylasen_Ausgabe1_0.docx','Monografie_Getreideamylasen_Ausgabe1_0.pdf','pruefung/Inhaltliche_Schlussabnahme_Fassung36.pdf'}|{'monografie/abbildungen/'+n for n in figures}
 check('Nur eigene bezeichnete Dokumente und neue wissenschaftliche Darstellungen',all(n in allowed_binaries for n in expected if Path(n).suffix.lower() in ['.pdf','.docx','.jpg','.jpeg','.png','.xlsx']))
 check('Keine Kritik Korrespondenz oder historischen Manuskriptkopien',not any(re.search(r'kritik|ausgang_fassung|korrespondenz|manuskriptfassung\d',n,re.I) for n in expected))
 texts=[(R/n).read_text(errors='replace') for n in expected if Path(n).suffix in ['.py','.json','.md','.csv','.xml','.cff','.txt']]
 check('Keine privaten Ablagekennungen oder Zugriffstoken',not any(re.search(r'libfile_[a-f0-9]{16,}|/workspace/scratch/[a-f0-9]{12}|/mnt[/]data/|[?&](?:sig|X-Amz-Signature)=',t,re.I) for t in texts))
 check('12 Kapitel 52 Tabellen und 92 nummerierte Gleichungen',len(cs)==12 and sum('table'in b for b in bs)+len(read('monografie/anhang_tabellen.json'))==52 and [b['equation_number'] for b in bs if 'eq'in b]==list(range(1,93)))
 refs=read('quellen/literatur.json')
 check('85 aktuelle Literaturangaben mit Quellenzugriff',len(refs)==85 and refs==[dict(number=i+1,reference=b['ref'],url=b['url'],access=b.get('access','')) for i,b in enumerate(cs[11]['blocks'])])
 check('101 Quellenregistereinträge erhalten',len(list(csv.DictReader((R/'quellen/quellenregister.csv').open(encoding='utf-8-sig'))))==101)
 cff=yaml.safe_load((R/'CITATION.cff').read_text());release=read('RELEASE.json');meta=read('monografie/manuskript_meta.json')
 check('Ausgabe und Autor in Text und Metadaten konsistent',cff['version']==cff['preferred-citation']['version']==release['package_version']==meta['version']=='1.0' and release['author']==meta['author']=='Stefan Traublinger')
 check('Keine erfundenen DOI ORCID Institution oder Veröffentlichung',all(release[k] is None for k in ['doi','orcid','institution']) and release['published'] is False and 'date-released' not in cff)
 docpath=R/'Monografie_Getreideamylasen_Ausgabe1_0.docx';doc=Document(docpath);preserved=read('pruefung/ERHALTUNG.json')
 with ZipFile(docpath) as z:
  check('Keine Kommentarbestandteile in der Lesefassung',not any('comments' in n.lower() for n in z.namelist()) and not doc.element.xpath('.//w:commentRangeStart|.//w:commentRangeEnd|.//w:commentReference'))
  check('Keine verbliebenen Änderungsmarkierungen',not doc.element.xpath('.//w:ins|.//w:del|.//w:moveFrom|.//w:moveTo'))
  media=sorted(sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith('word/media/'))
 omml=[sha256(etree.tostring(x,method='c14n')).hexdigest() for x in doc.element.xpath('.//m:oMathPara')]
 check('92 native Gleichungen exakt erhalten',omml==preserved['native_omml_sha256'] and len(omml)==92)
 check('Elf Medien einschließlich Unterschrift exakt erhalten',media==preserved['media_sha256'] and len(media)==11)
 table_hashes=[sha256(json.dumps([[c.text for c in row.cells] for row in t.rows],ensure_ascii=False).encode()).hexdigest() for t in doc.tables]
 check('Alle 52 Word Tabellen exakt erhalten',table_hashes==preserved['word_table_text_sha256'] and len(table_hashes)==52)
 source_hashes=[sha256(json.dumps(b['table'],ensure_ascii=False,sort_keys=True).encode()).hexdigest() for b in bs if 'table'in b]
 check('Alle 41 Haupttabellen in der Textquelle exakt erhalten',source_hashes==preserved['source_table_sha256'] and len(source_hashes)==41)
 check('Zehn wissenschaftliche Abbildungen vorhanden',len(figures)==10 and all((R/'monografie/abbildungen'/n).is_file() for n in figures))
 check('Word Autor und Ausgabenrevision korrekt',doc.core_properties.author==doc.core_properties.last_modified_by=='Stefan Traublinger' and doc.core_properties.revision==1)
 check('Titel Impressum und Zitierempfehlung vorhanden',all(any(t in p.text for p in doc.paragraphs) for t in ['Impressum','Ausgabe 1.0','info@varietas-verlag.de','Zitierempfehlung']))
 check('Zusammenfassung unverändert',all(any(p.text==a for p in doc.paragraphs) for a in meta['abstract'].split('\n\n')))
 check('Keine internen Redaktionskommentare in der strukturierten Quelle',not any(any(k.startswith('revision') and k.endswith('_comment') for k in b) for b in bs))
 pdf=fitz.open(R/'Monografie_Getreideamylasen_Ausgabe1_0.pdf')
 check('PDF mit 104 Seiten und korrektem Autor',len(pdf)==104 and pdf.metadata['author']=='Stefan Traublinger')
 check('104 aktuelle Fußzeilen',all('18. November 2026 | '+str(i+1)+' / 104' in ' '.join(p.get_text().split()) for i,p in enumerate(pdf)))
 data=list(csv.DictReader((R/'zertifizierte_fehlergrenzen/data/nominal_data.csv').open()))
 check('30 nominale Mehlwerte in vier Kurven',len(data)==30 and {x['curve'] for x in data}=={'Shintani_A','Shintani_F','H_STD','H_HIGH'})
 check('Mengen und Zeiteinheiten getrennt erhalten',all((x['time_unit'],x['value_unit'])==(('h','g maltose / g input flour') if x['curve'].startswith('Shintani') else ('min','mg maltose / g dry flour')) for x in data))
 check('Zwei Shintani Szenarioarchive ohne Pythonobjekte',all(all(not ar[k].dtype.hasobject for k in ar.files) for ar in [np.load(R/'fallstudie_shintani/pixel_scenarios.npz',allow_pickle=False),np.load(R/'fallstudie_shintani/F_pixel_scenarios.npz',allow_pickle=False)]))
 sources=read('quellen/datenquellen.json')
 check('Drei tatsächliche Datenquellengruppen benannt',{x['id'] for x in sources}=={'S','H','I'})
 check('Fremdanteile und Originallizenzgrenzen erhalten',sources[0]['original_license']=='CC-BY-NC-SA-4.0' and sources[2]['original_license']=='CC-BY-4.0' and 'Keine allgemeine freie Lizenz' in sources[1]['original_license'])
 check('Vier dokumentierte Volltextlücken erhalten',read('quellen/Lesetiefe_und_Grenzen.json')['remaining_full_text_gaps']==[27,28,30,85])
 check('Vier redaktionelle Änderungen ohne neue Ergebnisse',read('AENDERUNGEN.json')['scientific_results_changed'] is False and len(read('AENDERUNGEN.json')['changes'])==4)
 check('Lizenzzuordnung der einzelnen Monografieteile vorhanden',(R/'INHALTSRECHTE.csv').is_file() and 'Shintani' in (R/'INHALTSRECHTE.csv').read_text())
 return dict(status='PASS',version='1.0',files=len(expected),check_count=len(checks),checks=checks)

if __name__=='__main__':
 if not __debug__:raise RuntimeError('Ohne -O oder -OO ausführen.')
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write',action='store_true');a=ap.parse_args()
 if a.write:
  hashes=snapshot();(R/'MANIFEST_SHA256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2)+'\n')
  (R/'SHA256SUMS.txt').write_text(''.join(v+'  '+k+'\n' for k,v in hashes.items()))
  print(json.dumps(dict(status='WRITTEN',files=len(hashes))))
 else:print(json.dumps(verify(),ensure_ascii=False,indent=2))
