# GitHub-Paket veröffentlichen

Stefan Traublinger · Vorbereitungsstand 19. November 2026

## Was bereits vorbereitet ist

Das ZIP enthält den vollständigen Ordner `cereal-amylase-models`. Darin liegen die englische Projektbeschreibung mit deutscher Kurzfassung, die Zitiervorlage, der unveränderte wissenschaftliche Ausgangsbestand, die Rechenanleitung und die Lizenzabgrenzungen. Die Monografie ist mit https://zenodo.org/records/24378474 und DOI **10.5281/zenodo.24378474** verknüpft.

DNB-Kataloglink und Amazon-Links können später ergänzt werden. Ein GitHub-Repository wurde noch nicht angelegt oder veröffentlicht. Dieses Paket benötigt keine zusätzliche DOI für die Monografie.

## Veröffentlichung mit GitHub Desktop unter Windows

1. ZIP vollständig entpacken. Der Ordner `cereal-amylase-models` ist der vorbereitete Repository-Inhalt.
2. GitHub Desktop von https://desktop.github.com/ installieren und mit dem eigenen GitHub-Konto verbinden.
3. In GitHub Desktop **File → New repository** wählen. Als Namen `cereal-amylase-models` verwenden. Einen neuen lokalen Zielordner wählen. Keine zusätzliche Lizenzvorlage und keine zusätzliche README erzeugen, weil diese Dateien bereits im Paket enthalten sind.
4. **Repository → Show in Explorer** öffnen. Den **Inhalt** des entpackten Ordners `cereal-amylase-models` in diesen Repository-Ordner kopieren. Es soll dort direkt eine `README.md` neben `research`, `docs`, `metadata`, `CITATION.cff` und den Einstiegsskripten liegen. Den vorbereiteten Ordner nicht als zusätzliche Unterebene hineinkopieren. Auch `.gitattributes` und `.gitignore` mitkopieren.
5. In GitHub Desktop alle übernommenen Dateien für den Commit auswählen. Als Zusammenfassung beispielsweise **Add computational companion and Zenodo citation** eintragen und den Commit ausführen.
6. **Publish repository** wählen. Für öffentliche Zugänglichkeit die Option **Keep this code private** ausschalten und die Veröffentlichung ausführen.
7. Auf der GitHub-Seite unter **About** diese Beschreibung eintragen: `Computational companion to the monograph on cereal amylases and starch hydrolysis by Stefan Traublinger.` Als Website https://zenodo.org/records/24378474 verwenden. Passende Topics: `starch-hydrolysis`, `enzyme-kinetics`, `cereal-chemistry`, `identifiability`, `triticale`.
8. Unter **Releases → Draft a new release** den Tag `v1.0.0` anlegen. Titel: **Computational companion 1.0.0**. Den Text aus `docs/RELEASE_NOTES.md` übernehmen und den ersten Release veröffentlichen.

Die vorhandenen Lizenzen bitte beibehalten. Das Rechenpaket enthält ausdrücklich markierte Fremdanteile mit eigenen Bedingungen. Eine zusätzliche pauschale MIT-Lizenz würde diese Zuordnung nicht ersetzen.

Die `.gitattributes` schützt die unveränderten wissenschaftlichen Quelldateien vor automatischer Zeilenende-Umwandlung unter Windows. Dadurch bleiben ihre ursprünglichen Prüfsummen auch nach einem Git-Checkout erhalten.

## Nachrechnen

Im Repository-Ordner mit Python 3.12:

```bat
python -m venv .venv
.venv\Scripts\activate.bat
python -m pip install -r requirements.txt
python run_reproduction.py
```

Ergebnisse stehen unter `outputs/latest/`. Für den schnellen Integritätscheck des Forschungsbestands: `python run_reproduction.py --verify-only`. Für die Prüfsummen des unveränderten gelieferten Repository-Snapshots: `python verify_repository.py`.

Die vorhandene Laufzeit wurde geprüft. Eine frische Installation sämtlicher Abhängigkeiten und ein Lauf auf GitHub selbst wurden noch nicht geprüft. Der Rechenlauf ersetzt keine externe Fachbegutachtung oder Laborvalidierung.

## Anschließend ergänzen

Sobald der DNB-Katalogeintrag verfügbar ist, dessen beständigen Link in der README ergänzen. Amazon-Links zur deutschen und englischen Ausgabe können daneben aufgenommen werden. Für solche redaktionellen Ergänzungen ist keine neue Forschungsarbeit erforderlich; für einen neuen Release sind die Prüfsummen des Repository-Snapshots zu erneuern.

Den GitHub-Link anschließend als zugehörige Ressource im vorhandenen Zenodo-Eintrag hinterlegen. Die automatische GitHub-Zenodo-Archivierung wurde nicht eingerichtet. Die vorhandene DOI bezeichnet die Monografie; ein später separat archivierter Softwarestand wäre ein eigenes Forschungsobjekt.

## Offizielle Anleitungen

- Erstes Repository mit GitHub Desktop: https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop
- Commits und Änderungen prüfen: https://docs.github.com/en/desktop/making-changes-in-a-branch/committing-and-reviewing-changes-to-your-project-in-github-desktop
- Releases: https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases
- Zitiervorlage: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
