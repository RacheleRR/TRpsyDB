import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from backend.db.database import init_db, get_connection

LOCI = [
    ("chr4:3074876-3074938",    "chr4",  3074876,  3074938,  "CAG",   "HTT",    1, "HD",      "Huntington disease repeat"),
    ("chr3:63912684-63912714",  "chr3",  63912684, 63912714, "AAGGG", "RFC1",   1, "CANVAS",  "RFC1 CANVAS repeat"),
    ("chrX:147582107-147582162","chrX",  147582107,147582162,"CGG",   "FMR1",   1, "FXS",     "Fragile X CGG repeat"),
    ("chr19:45770205-45770264", "chr19", 45770205, 45770264, "CTG",   "DMPK",   1, "DM1",     "Myotonic dystrophy 1"),
    ("chr22:46191234-46191304", "chr22", 46191234, 46191304, "GAA",   "FXN",    1, "FRDA",    "Friedreich ataxia GAA repeat"),
    ("chr9:27573529-27573546",  "chr9",  27573529, 27573546, "GGGGCC","C9orf72",1, "ALS,FTD", "C9orf72 G4C2 repeat"),
]

def seed():
    init_db()
    with get_connection() as conn:
        for row in LOCI:
            conn.execute(
                """INSERT OR IGNORE INTO loci
                   (locus_id,chrom,start,end,motif,gene_symbol,evidence_tier,phenotypes,common_name)
                   VALUES (?,?,?,?,?,?,?,?,?)""", row)
        n = conn.execute("SELECT COUNT(*) FROM loci").fetchone()[0]
    print(f"[seed] {n} loci in database.")

if __name__ == "__main__":
    seed()
