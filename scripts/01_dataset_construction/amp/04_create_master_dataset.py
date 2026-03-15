from __future__ import annotations
import pandas as pd
from config import DATA_INTERMEDIATE_DIR

def main() -> None:
    dfs = []
    for db in ["CAMP", "DBAASP", "dbAMP3", "DRAMP"]:
        f = DATA_INTERMEDIATE_DIR / f"{db}_03_mic_filtered.parquet"
        if f.exists():
            dfs.append(pd.read_parquet(f))
        else:
            print(f"⚠️ missing: {f.name}")

    if not dfs:
        raise SystemExit("No DBs loaded.")

    master = pd.concat(dfs, ignore_index=True)
    out = DATA_INTERMEDIATE_DIR / "DB_MASTER_04.parquet"
    master.to_parquet(out, index=False)

    print(f"✅ master built: {len(master)} rows -> {out.name}")

if __name__ == "__main__":
    main()