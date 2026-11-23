"""Download the reference files of the catalogue (MetaNetX, KEGG REST, ModelSEED, BiGG, MetaCyc via gapseq, Rhea,
ECDomainMiner, TCDB, ChEBI) into DATA_DIR/<database>/... (needs internet access; run it where the sources are
reachable, e.g. the HPC access node). KEGG data are subject to KEGG's terms of use (academic use)."""
import concurrent.futures as cf
import logging
import os
import tempfile
import time
import urllib.request
from pathlib import Path

from .catalogue import DATA_CATALOGUE

log = logging.getLogger(__name__)


def plan(output_dir, sources=None):
    out = []
    for group in DATA_CATALOGUE:
        if sources and group.database not in sources:
            continue
        out += [(u, p) for u, p in group.iter_downloads(output_dir)]
    return sorted(set(out))


def fetch(url, path, retries=3, timeout=120):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for i in range(retries):
        tmp = None
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "IMP-metabolic"})
            with urllib.request.urlopen(req, timeout=timeout) as r, \
                    tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as fh:
                tmp = fh.name
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            os.replace(tmp, path)
            return {"ok": True, "url": url, "path": str(path), "bytes": path.stat().st_size}
        except Exception as e:                          # noqa: BLE001 (report every failure, retry)
            if tmp and os.path.exists(tmp):
                os.remove(tmp)
            err = repr(e)
            time.sleep(2 + 5 * i)
    return {"ok": False, "url": url, "path": str(path), "error": err}


def download(output_dir, sources=None, workers=5, skip_existing=True):
    jobs = plan(output_dir, sources)
    if skip_existing:
        jobs = [(u, p) for u, p in jobs if not (Path(p).exists() and Path(p).stat().st_size > 0)]
    res = []
    with cf.ThreadPoolExecutor(max_workers=max(1, min(workers, len(jobs) or 1))) as ex:
        for i, r in enumerate(ex.map(lambda j: fetch(*j), jobs), 1):
            res.append(r)
            log.info("%d/%d %s %s", i, len(jobs), "ok  " if r["ok"] else "FAIL", r["url"])
    return res
