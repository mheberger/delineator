import tarfile
import zipfile
from pathlib import Path

def top_files(path, n=20):
    path = Path(path)
    if path.suffixes[-2:] == ['.tar', '.gz'] or path.suffix == '.gz':
        with tarfile.open(path) as t:
            items = [(m.size, m.name) for m in t.getmembers() if m.isfile()]
    else:
        with zipfile.ZipFile(path) as z:
            items = [(i.file_size, i.filename) for i in z.infolist()]

    total = sum(s for s, _ in items)
    print(f"\n{path.name}  —  {total/1e6:.1f} MB uncompressed, {len(items)} files")
    for size, name in sorted(items, reverse=True)[:n]:
        print(f"  {size/1e6:>8.1f} MB  {name}")

for f in Path("dist").iterdir():
    top_files(f)