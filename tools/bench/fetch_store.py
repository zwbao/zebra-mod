#!/usr/bin/env python3
"""Download the GA4GH phenopacket-store release used by the benchmark, and record what was fetched.

    python3 tools/bench/fetch_store.py [--tag 0.1.27] [--out DIR]

The release is resolved through the GitHub API (latest release unless --tag is
given); the asset `all_phenopackets.zip` is downloaded through zebra.http's
opener (same proxy handling and User-Agent as every other zebra request), its
sha256 is computed locally and compared with the digest GitHub publishes for
the asset. The archive is never extracted: tools/bench reads members straight
from the zip, so a hostile member path cannot write anywhere.

Also fetches MONDO's SSSOM mapping table at the latest MONDO release tag (pinned
by tag in the URL), used to join ORPHA and MONDO hits to the OMIM truth through
exactMatch rows only; written to DIR/../mondo/<tag>/ with its own manifest.

Writes DIR/<tag>/all_phenopackets.zip and DIR/<tag>/manifest.json:
  {repo, tag, published_at, asset_url, bytes, sha256, github_digest,
   digest_matches, retrieved_at, members (json files in the zip)}
DIR defaults to $ZEBRA_CACHE_DIR/bench/phenopacket-store (~/.cache/zebra-mod/...).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import shutil  # noqa: E402

from zebra import hpo_local  # noqa: E402
from zebra.http import USER_AGENT, _opener, cache_dir, now_iso, request  # noqa: E402

REPO = "monarch-initiative/phenopacket-store"
ASSET = "all_phenopackets.zip"
MONDO_REPO = "monarch-initiative/mondo"
# MONDO's own SSSOM mapping set (the exactMatch rows Monarch's /mappings endpoint serves),
# read at the release tag so the file is pinned; MONDO publishes no release asset for it.
MONDO_SSSOM = "https://raw.githubusercontent.com/monarch-initiative/mondo/{tag}/src/ontology/mappings/mondo.sssom.tsv"


def default_dir() -> Path:
    return cache_dir() / "bench" / "phenopacket-store"


def release_info(tag: str = None) -> dict:
    url = (f"https://api.github.com/repos/{REPO}/releases/tags/{tag}" if tag
           else f"https://api.github.com/repos/{REPO}/releases/latest")
    resp = request(url, source="GitHub releases", cache_ttl=0, accept="application/vnd.github+json")
    data = json.loads(resp.text)
    asset = next((a for a in data.get("assets") or [] if a.get("name") == ASSET), None)
    if asset is None:
        raise SystemExit(f"release {data.get('tag_name')} has no asset {ASSET}")
    return {"repo": REPO, "tag": data["tag_name"], "published_at": data.get("published_at"),
            "asset_url": asset["browser_download_url"], "bytes_published": asset.get("size"),
            "github_digest": asset.get("digest"), "release_api_url": url}


def download(url: str, dest: Path) -> str:
    """Stream `url` to `dest`; returns the sha256 hex digest."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    sha = hashlib.sha256()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with _opener().open(req, timeout=900) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            sha.update(chunk)
            fh.write(chunk)
    os.replace(tmp, dest)
    return sha.hexdigest()


def sha256_file(path: Path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()


def fetch(tag: str = None, out: Path = None, force: bool = False) -> dict:
    info = release_info(tag)
    target = (out or default_dir()) / info["tag"]
    target.mkdir(parents=True, exist_ok=True)
    zpath = target / ASSET
    if zpath.exists() and not force:
        digest = sha256_file(zpath)
        retrieved = json.loads((target / "manifest.json").read_text("utf-8")).get("retrieved_at") \
            if (target / "manifest.json").exists() else now_iso()
    else:
        digest = download(info["asset_url"], zpath)
        retrieved = now_iso()
    expected = (info.get("github_digest") or "").split(":", 1)[-1] or None
    with zipfile.ZipFile(zpath) as zf:
        members = sum(1 for n in zf.namelist() if n.endswith(".json"))
    manifest = dict(info, bytes=zpath.stat().st_size, sha256=digest, digest_matches=(expected == digest) if expected else None,
                    retrieved_at=retrieved, members=members, path=str(zpath))
    if expected and expected != digest:
        raise SystemExit(f"sha256 mismatch for {zpath}: got {digest}, GitHub says {expected}")
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", "utf-8")
    return manifest


def fetch_mappings(tag: str = None, out: Path = None, force: bool = False) -> dict:
    """MONDO's SSSOM mappings at a release tag (latest by default), with sha256 and retrieval time."""
    if tag is None:
        resp = request(f"https://api.github.com/repos/{MONDO_REPO}/releases/latest", source="GitHub releases",
                       cache_ttl=0, accept="application/vnd.github+json")
        tag = json.loads(resp.text)["tag_name"]
    target = (out or default_dir().parent / "mondo") / tag
    target.mkdir(parents=True, exist_ok=True)
    path = target / "mondo.sssom.tsv"
    url = MONDO_SSSOM.format(tag=tag)
    if path.exists() and not force and (target / "manifest.json").exists():
        old = json.loads((target / "manifest.json").read_text("utf-8"))
        if old.get("sha256") == sha256_file(path):
            return old
    digest = download(url, path)
    with open(path, encoding="utf-8") as fh:
        head = fh.read(4096)
    if "subject_id" not in head and "curie_map" not in head:
        raise SystemExit(f"{url} did not return an SSSOM table")
    manifest = {"repo": MONDO_REPO, "tag": tag, "url": url, "bytes": path.stat().st_size, "sha256": digest,
                "retrieved_at": now_iso(), "path": str(path)}
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", "utf-8")
    return manifest


HPO_REPO = "obophenotype/human-phenotype-ontology"
HPO_FILES = ("hp.json", "phenotype.hpoa", "genes_to_phenotype.txt")


def fetch_hpo(tag: str, out: Path = None, force: bool = False) -> dict:
    """The three HPO release files `zebra hpo fetch` installs, pinned to a release tag (e.g. v2026-09-01),
    each checked against the sha256 digest GitHub publishes for the asset. Point ZEBRA_HPO_DIR at the
    directory to rank with exactly the release the benchmark used."""
    resp = request(f"https://api.github.com/repos/{HPO_REPO}/releases/tags/{tag}", source="GitHub releases",
                   cache_ttl=0, accept="application/vnd.github+json")
    data = json.loads(resp.text)
    assets = {a["name"]: a for a in data.get("assets") or []}
    target = (out or default_dir().parent / "hpo") / tag
    target.mkdir(parents=True, exist_ok=True)
    files = {}
    for name in HPO_FILES:
        a = assets.get(name)
        if a is None:
            raise SystemExit(f"HPO release {tag} has no asset {name}")
        path = target / name
        expected = (a.get("digest") or "").split(":", 1)[-1] or None
        installed = hpo_local.data_dir() / name
        if path.exists() and not force:
            digest = sha256_file(path)
        elif expected and installed.exists() and sha256_file(installed) == expected:
            # the copy `zebra hpo fetch` installed IS this release asset (same sha256): no download needed
            shutil.copyfile(installed, path)
            digest = expected
        else:
            digest = download(a["browser_download_url"], path)
        if expected and digest != expected:
            raise SystemExit(f"sha256 mismatch for {path}: got {digest}, GitHub says {expected}")
        files[name] = {"bytes": path.stat().st_size, "sha256": digest, "url": a["browser_download_url"],
                       "digest_matches": expected == digest if expected else None}
    manifest = {"repo": HPO_REPO, "tag": tag, "published_at": data.get("published_at"), "files": files,
                "retrieved_at": now_iso(), "dir": str(target)}
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", "utf-8")
    return manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tag", help="release tag (default: latest)")
    ap.add_argument("--out", type=Path, help="directory (default: $ZEBRA_CACHE_DIR/bench/phenopacket-store)")
    ap.add_argument("--force", action="store_true", help="download again even if the zip is present")
    ap.add_argument("--mondo-tag", help="MONDO release tag for the SSSOM mappings (default: latest)")
    ap.add_argument("--hpo-tag", help="also fetch this HPO release (e.g. v2026-09-01, what the benchmark used)")
    args = ap.parse_args(argv)
    print(json.dumps(fetch(args.tag, args.out, args.force), indent=2))
    print(json.dumps(fetch_mappings(args.mondo_tag, args.out.parent / "mondo" if args.out else None, args.force),
                     indent=2))
    if args.hpo_tag:
        print(json.dumps(fetch_hpo(args.hpo_tag, args.out.parent / "hpo" if args.out else None, args.force), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
