#!/usr/bin/env python3
"""Verify packaged hashes, then independently rebuild outputs in a temp directory."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT = Path(__file__).resolve().parents[1]

def require(ok, message):
    if not ok:
        raise ValueError(message)


for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
    expected, name = line.split('  ',1)
    require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected, f'Checksum mismatch: {name}')

with tempfile.TemporaryDirectory() as directory:
    subprocess.run([sys.executable,str(ROOT/'scripts/build_features.py'),'--out',directory],check=True,stdout=subprocess.DEVNULL)
    for baseline in sorted((ROOT/'outputs').iterdir()):
        rebuilt=Path(directory)/baseline.name
        if baseline.name=='validation.json':
            a=json.loads(baseline.read_text()); b=json.loads(rebuilt.read_text())
            a.pop('runtime'); b.pop('runtime')
            require(a==b,'Validation results changed')
        else:
            require(baseline.read_bytes()==rebuilt.read_bytes(),f'Rebuild differs: {baseline.name}')

print('PASS: all package checksums; input preserved; all rebuilt TSV/schema bytes identical; validation results identical (runtime metadata excluded).')
