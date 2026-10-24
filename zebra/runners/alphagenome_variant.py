"""One AlphaGenome predict_variant call, run by zebra in a child interpreter.

Usage: python alphagenome_variant.py '<json>' with chrom (chr-prefixed), pos (1-based),
ref, alt, ontology (list of UBERON/CL CURIEs), outputs (list of OutputType names),
timeout. Needs ALPHAGENOME_API_KEY and the official client (`pip install alphagenome`).
Prints one JSON summary on stdout; never prints the key.
"""

import json
import os
import sys

WIDTH = 16_384   # the API's smallest supported interval
LOCAL = 500      # +/- bases around the variant for the local readout


def main() -> int:
    args = json.loads(sys.argv[1])
    summary = {"status": "error", "inputs": args}
    try:
        from importlib import metadata

        import numpy as np
        from alphagenome.data import genome
        from alphagenome.models import dna_client

        version = metadata.version("alphagenome")
        summary["model"] = {"id": "alphagenome hosted API", "checkpoint": None, "version": version}
        model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=float(args.get("timeout", 120)))
        start = max(0, args["pos"] - WIDTH // 2)
        interval = genome.Interval(chromosome=args["chrom"], start=start, end=start + WIDTH)
        variant = genome.Variant(chromosome=args["chrom"], position=args["pos"],
                                 reference_bases=args["ref"], alternate_bases=args["alt"])
        output_types = []
        for name in args["outputs"]:
            ot = getattr(dna_client.OutputType, name, None)
            if ot is None:
                raise ValueError(f"unknown AlphaGenome output {name!r}")
            output_types.append(ot)
        out = model.predict_variant(interval=interval, organism=dna_client.Organism.HOMO_SAPIENS,
                                    variant=variant, requested_outputs=output_types,
                                    ontology_terms=args["ontology"])
        heads = {}
        centre = args["pos"] - 1 - start
        for name in args["outputs"]:
            r_td = getattr(out.reference, name.lower(), None)
            a_td = getattr(out.alternate, name.lower(), None)
            if r_td is None or a_td is None:
                continue
            r = np.asarray(r_td.values, dtype=float)
            a = np.asarray(a_td.values, dtype=float)
            if r.ndim == 1:
                r, a = r[:, None], a[:, None]
            lfc = np.log2((a + 1.0) / (r + 1.0))
            # the arrays may be binned (resolution > 1 bp): map the local window to rows
            res = max(1, WIDTH // max(1, r.shape[0]))
            lo, hi = max(0, (centre - LOCAL) // res), min(r.shape[0], (centre + LOCAL) // res + 1)
            local = lfc[lo:hi]
            i, j = np.unravel_index(int(np.argmax(np.abs(local))), local.shape) if local.size else (0, 0)
            names = None
            try:
                names = list(r_td.metadata["name"])
            except Exception:
                pass
            heads[name] = {
                "num_tracks": int(r.shape[1]),
                "resolution_bp": int(res),
                "mean_diff": float(np.mean(a - r)),
                "log2fc": float(np.mean(lfc)),
                "local_max_abs_log2fc": float(local[i, j]) if local.size else None,
                "local_track": (names[j] if names and j < len(names) else int(j)),
                "local_offset_bp": int((lo + i) * res - centre) if local.size else None,
            }
        if not heads:
            raise RuntimeError("predict_variant returned none of the requested outputs")
        summary.update(status="ok", head_stats=heads, interval={"start": start, "end": start + WIDTH})
    except Exception as err:
        summary["error"] = f"{type(err).__name__}: {err}"
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
