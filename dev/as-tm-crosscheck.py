"""Score a tidymarkers allele-specific panel with Biopython's nearest-neighbour Tm.

Reads the panel CSV written by ``dev/as-tm-crosscheck.R`` (columns ``primer``,
``t_ref``, ``t_alt``) and writes the same rows plus ``tm_ref``, ``tm_alt`` and
``dtm = tm_ref - tm_alt``.

Why Biopython: ``Bio.SeqUtils.MeltingTemp.Tm_NN`` applies the 3'-terminal
mismatch table (SantaLucia & Peyret 2001) when it is given the second strand as
``c_seq`` and ``strict = False``.  ``ntthal`` truncates instead, so it cannot
see a terminal mismatch at all and is not a comparator here.

The second strand goes in 3'->5' (``c_seq = target_5p3[::-1]``) and
``strict = False`` is required for a mismatch.  There is no ``table`` or
``temp_c`` keyword; the table argument is ``nn_table=<dict>`` and the default
(``DNA_NN3``, Allawi & SantaLucia 1997) is left in place.

Usage: python3 dev/as-tm-crosscheck.py <panel.csv> <out.csv>
"""

import csv
import sys

from Bio.SeqUtils import MeltingTemp as mt

# Matched to reaction_conditions(): Na 50 mM, Mg 1.5 mM, dNTPs 0.6 mM, and
# dnac1 = dnac2 = 25 nM (the PCR probe concentration `pcr_efficiency()` uses).
SETTINGS = dict(dnac1=25, dnac2=25, Na=50, K=0, Tris=0, Mg=1.5, dNTPs=0.6,
                saltcorr=5)


def main(argv):
    if len(argv) != 3:
        sys.exit("usage: as-tm-crosscheck.py <panel.csv> <out.csv>")
    rows = list(csv.DictReader(open(argv[1], newline="", encoding="utf-8")))
    if not rows:
        sys.exit("empty panel")
    for r in rows:
        # c_seq is the complement strand in 3'->5' orientation.
        tm_ref = mt.Tm_NN(r["primer"], c_seq=r["t_ref"][::-1], strict=False,
                          **SETTINGS)
        tm_alt = mt.Tm_NN(r["primer"], c_seq=r["t_alt"][::-1], strict=False,
                          **SETTINGS)
        r["tm_ref"] = round(tm_ref, 4)
        r["tm_alt"] = round(tm_alt, 4)
        r["dtm"] = round(tm_ref - tm_alt, 4)
    with open(argv[2], "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv)
