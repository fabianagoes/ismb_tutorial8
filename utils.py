import os
from pathlib import Path

def get_all_subdirs(root):
    return [str(p) for p in Path(root).rglob("*") if p.is_dir() and os.path.exists(os.path.join(p,'train.csv'))]