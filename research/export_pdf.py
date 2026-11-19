"""Export the editable Word source with local LibreOffice to a separate folder."""
from pathlib import Path
import argparse,shutil,subprocess,tempfile
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output-dir',type=Path,default=R/'ergebnisse/laufprotokolle/pdf_neusatz')
a=ap.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
source=R/'Monografie_Getreideamylasen_Ausgabe1_0.docx'
if out==R:raise ValueError('Ein separates Ausgabeverzeichnis verwenden.')
office=shutil.which('libreoffice') or shutil.which('soffice')
if not office:raise RuntimeError('LibreOffice ist nicht installiert oder nicht im Suchpfad.')
with tempfile.TemporaryDirectory(prefix='getreide_pdf_') as td:
 subprocess.run([office,'-env:UserInstallation='+Path(td).as_uri(),'--headless','--convert-to','pdf','--outdir',str(out),str(source)],check=True)
target=out/(source.stem+'.pdf')
if not target.is_file() or target.stat().st_size==0:raise RuntimeError('Keine PDF erzeugt.')
print(target)
