# Rebuild the manuscript

The editable manuscript and PDF share one layout. The scripts convert `../manuscript.md` and its relative figure paths to an editable Word document, then export that document to tagged PDF. The scientific analysis and figure generation have their own scripts and data.

The delivered files were built using Node.js 24.22.3, Marked 17.0.5, Python 3.12.14, python-docx 1.2.0, Pillow 13.0.0, and Microsoft Word 2021 for Windows. The document uses Georgia and Arial.

Run these PowerShell commands from the manuscript directory, which contains `manuscript.md`, `figures`, and `build`. Set the paths to Node and Python first if they are outside your PATH.

```powershell
$scratch = Join-Path ([System.IO.Path]::GetTempPath()) ('tbxt-preprint-' + [guid]::NewGuid())
New-Item -ItemType Directory -Path $scratch | Out-Null
npm install --prefix $scratch marked@17.0.5
python -m pip install python-docx==1.2.0 Pillow==13.0.0
node build/build.mjs manuscript.md $scratch (Join-Path $scratch 'node_modules')
python build/build_docx.py (Join-Path $scratch 'article-tokens.json') manuscript.md TBXT-preprint.docx
& ./build/export_word.ps1 -DocumentPath TBXT-preprint.docx -PdfPath TBXT-preprint.pdf
```

The tokenization step also writes `article.html` in the temporary directory. It provides a browser-readable view with print styling. The delivered PDF is the Word export, which retains the layout and pagination of the editable manuscript.
