"""
Create the SHA-256 hashes for `delienator`'s download function
using the `pooch` library

Matthew Heberger, June 2026
"""

import pooch
from pathlib import Path

registry = {}
data_dir = Path(r"C:\Users\mheberger\AppData\Local\delineator")

with open("hashes.txt", "w") as o:
    
    for f in sorted(data_dir.rglob("*")):
        if f.is_file():
            file = str(f.relative_to(data_dir))
            checksum = pooch.file_hash(str(f))  # sha256: prefix included
            print( f'    "{file}": "{checksum}",')
            o.write(f'{checksum}  {file}\n')
