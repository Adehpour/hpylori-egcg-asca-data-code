"""Generate checksums for all distributed public-release files.

Run after changing code, metadata, data or the final supplementary workbook.
"""
from pathlib import Path
from hashlib import sha256

ROOT=Path(__file__).resolve().parent
entries=[]
for f in ROOT.rglob('*'):
    if not f.is_file() or f.name == 'MANIFEST_SHA256.txt':
        continue
    if any(k in f.parts for k in ('__pycache__', '.git')) or f.suffix in ('.pyc', '.ndjson'):
        continue
    entries.append((str(f.relative_to(ROOT)), sha256(f.read_bytes()).hexdigest()))
contents='# SHA-256 of all distributed release files\n'
contents+=''.join(f'{digest}  {path}\n' for path,digest in sorted(entries))
(ROOT/'MANIFEST_SHA256.txt').write_text(contents,encoding='utf-8')
print(f'Updated manifest with {len(entries)} distributed files.')
