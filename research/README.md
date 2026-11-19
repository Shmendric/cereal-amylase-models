# Getreideamylasen und Stärkehydrolyse

Stefan Traublinger · Digitale Erstausgabe 1.0 · 18. November 2026

Die Monografie untersucht, welche Größen Maltosekurven in ausdrücklich festgelegten kinetischen Modellen bestimmen und welche zusätzliche Kalibrierung ihre chemische Deutung benötigt. Sie verbindet mathematische Nachweise mit der retrospektiven Auswertung von vier digitalisierten Weizenmehlkurven aus zwei Publikationsfamilien. Ein separater Vergleich verwendet zwanzig Paare publizierter Aktivitätsmittelwerte aus fünf Laboren. Der Triticale-Teil ist ein begründeter Messplan. Verbesserte Broteigenschaften und ein Beitrag zur Ernährungssicherung sind nicht experimentell nachgewiesen.

Die Lesefassung umfasst 104 Seiten. Sie enthält 41 Haupttabellen, elf technische Anhangstabellen, 92 native Gleichungen, zehn wissenschaftliche Abbildungen und 85 Literaturangaben. Die Word-Datei ist die bearbeitbare Satzquelle, die gleichnamige PDF die feste Lesefassung. Beide tragen die Ausgabe 1.0. Grundlage ist die intern geprüfte Manuskriptfassung 36. Originalartikel, Originalabbildungen, persönliche Rückmeldungen und historische Manuskriptkopien gehören nicht zu dieser Auswahl.

## Einstieg

Für die Lektüre Monografie_Getreideamylasen_Ausgabe1_0.pdf öffnen. Für die Rechnungen Python 3.12 verwenden und im entpackten Ordner ausführen:

    python -m pip install -r requirements.txt
    python verify_release.py
    python reproduce.py

Der Standardlauf prüft zuerst die Versanddateien und arbeitet anschließend in einer temporären Kopie. Er rechnet die nominalen E W P2 Anpassungen neu, prüft die gespeicherten technischen Szenarien, vier rationale Zertifikate, Einzelpunkteinflüsse, INFOGEST-Laborwerte, synthetische Kalibrier- und Zugangsmodelle, Fehlerkovarianz sowie die neuen Formel-, Tabellen- und Konturkontrollen. Es werden keine Quellen aus dem Netz geladen. Neue Ergebnisse stehen ausschließlich unter ergebnisse/laufprotokolle. Der Einstieg funktioniert auch über einen absoluten Skriptpfad aus einem anderen Arbeitsverzeichnis.

Für einen kurzen Lauf mit den archivierten nominalen Parametern:

    python reproduce.py --skip-optimization

Die 1572 historischen Annahmeszenarien werden aus den gespeicherten Parametern nachgerechnet, aber nicht alle neu optimiert. Technische Szenarien sind keine biologischen Wiederholungen. Überschneidende Prüfbedingungen werden nicht als unabhängige Experimente gezählt. PRUEFSTATUS.json und VEROEFFENTLICHUNGSPRUEFUNG.md dokumentieren den tatsächlich ausgeführten Paketlauf und seine Grenzen.

## Dateien und Nachweise

| Bereich | Inhalt |
|---|---|
| Monografie_Getreideamylasen_Ausgabe1_0.docx und .pdf | Saubere Word-Quelle mit 92 nativen Gleichungen und feste Lesefassung |
| monografie | Strukturierter aktueller Text, Gleichungsdefinitionen, Anhang und eigene wissenschaftliche Abbildungen |
| quellen | 85 Literaturangaben mit dokumentierter Lesetiefe, Herkunft der tatsächlichen Daten, Lizenzbelege und Nachweisliste |
| schlussfassung17 bis schlussfassung33 | Ausgewählte, funktionsbezogene Rechenmodule mit ihren versionsbezogenen Dateinamen |
| zertifizierte_fehlergrenzen | Rationale Fehlerzertifikate und ihre Begründung |
| pruefung | Erhaltungsnachweis, zusätzliche Formelkontrollen und interner Abnahmebericht des Ausgangsstands 36 |
| LICENSE.md, DATEIRECHTE.csv und INHALTSRECHTE.csv | Nutzungsbedingungen pro Datei und Zuordnung gemischter Monografieteile |
| CITATION.cff und ZITIEREN.txt | Zitierempfehlung für Monografie und Begleitpaket |

Die als data/original bezeichneten Unterordner enthalten ursprüngliche technische Projekt-CSV-Dateien und Ableseprotokolle. Sie enthalten keine Originalpublikationen. Mehrfach abgelegte Daten sind Bestandteile bestehender Prüfpfade und keine zusätzlichen Datensätze. Die nominellen dreißig Mehlwerte sind keine dreißig unabhängigen Mehlchargen. Hackenberg verwendet mg Maltose je g trockenes Mehl und Minuten, Shintani g Maltose je g eingesetztes Mehl und Stunden. Diese Bezugsgrößen bleiben getrennt.

## Neusatz und Bearbeitung

Die Word-Datei lässt sich direkt bearbeiten. monografie/manuskript.json ist zusätzlich eine strukturierte Textquelle für Nachverfolgung und Tabellenkontrollen. Das Paket enthält keinen vollständigen Word-Generator aus dieser JSON-Datei. Nach Änderungen sind strukturierte Quelle und Word-Datei gemeinsam zu pflegen, wenn die Prüfungen weiterhin auf die bearbeitete Ausgabe angewendet werden sollen.

Mit installiertem LibreOffice kann die mitgelieferte Word-Datei neu als PDF ausgegeben werden:

    python export_pdf.py

Die Ausgabe liegt unter ergebnisse/laufprotokolle/pdf_neusatz. Die feste Versand-PDF wird dabei nicht überschrieben. Zeilen- und Seitenumbrüche hängen beim Neusatz von Schriftarten und Office-Version ab. Der archivierte Satz verwendet das bestehende Format 21,59 × 27,94 cm. Ein Amazon-Drucksatz ist eine eigene Ausgabe und nicht Bestandteil dieses digitalen Pakets.

## Nutzungsbedingungen

Eigene Beiträge stehen, soweit entsprechende Rechte bestehen, unter CC BY 4.0 mit Namensnennung Stefan Traublinger. Die ausdrücklich gekennzeichneten Shintani-Anteile behalten CC BY-NC-SA 4.0. Gemischte Dateien sind komponentenbezogen zu verwenden. Die Angaben zu Hackenberg betreffen technische Zahlenablesungen und eigene Aufbereitung, keine Freigabe des Originalartikels, seiner Abbildungen oder seiner Dissertation. An bloßen wissenschaftlichen Tatsachen wird kein Ausschließlichkeitsrecht behauptet. LICENSE.md und die beiden Rechteverzeichnisse sind gemeinsam zu beachten.

## Wissenschaftliche Grenzen

Die Abnahme und alle Rechenkontrollen gehören demselben KI-gestützten Arbeitszusammenhang an. Eine menschliche Doppelablesung, externe Fachbegutachtung, experimentelle Mehlkalibrierung und unabhängige Materialvalidierung sind nicht abgeschlossen. Die Volltexte [27], [28], [30] und [85] bleiben nach dem dokumentierten Zugriffsstand offen. Historische Erstmaligkeit wird nicht als nachgewiesen beansprucht.

Die vier globalen Exponentialzertifikate gelten für die bezeichnete Modellklasse und die nominale Punktmenge. Sie sind keine experimentellen Unsicherheitsgrenzen und keine Vorhersagegarantie. Andere numerische Modellkandidaten sind nicht automatisch global optimal. Die INFOGEST-Temperaturumrechnung ist von nativer Mehlkalibrierung getrennt. Die Programme prüfen Modellumsetzung und die ausgewiesenen Datenfälle, keine allgemeine chemische Prozessgültigkeit.

Der beigefügte F36-Abnahmebericht ist ein datierter Nachweis des Ausgangsstands. Darin erwähnte historische Arbeitsarchivpfade sind nicht sämtlich Teil dieses öffentlichen Pakets. Der aktuelle Paketlauf und die Erhaltung der wissenschaftlichen Inhalte sind gesondert dokumentiert. Die vorhandene Laufzeitumgebung wurde verwendet, eine Neuinstallation der Abhängigkeiten wurde nicht getestet.

Die beiden Manifestdateien selbst, Python-Zwischenspeicher und ergebnisse/laufprotokolle sind vom Inhaltsmanifest ausgenommen. DATEIRECHTE.csv kann sich nicht selbst hashen und verweist in ihrer eigenen Zeile auf das Manifest. Veröffentlichungskennung und DOI sind noch nicht vergeben. Es wurde nichts öffentlich hochgeladen.
