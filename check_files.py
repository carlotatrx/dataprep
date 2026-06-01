import os

root = "/scratch3/PALAEO-RA/daily_data/final/"

log = open("/scratch3/PALAEO-RA/daily_data/final/filename_warnings.txt", "w", encoding="utf-8")

for dirpath, _, filenames in os.walk(root):
    for fn in filenames:
        if not fn.endswith(".tsv"):
            continue

        path = os.path.join(dirpath, fn)

        try:
            with open(path, "r", encoding="utf-8") as f:
                id_ = name = None
                for _ in range(11):  # metaheader only
                    line = f.readline()
                    if not line:
                        break
                    if line.startswith("ID\t"):
                        id_ = line.strip().split("\t", 1)[1]
                    elif line.startswith("Name\t"):
                        name = line.strip().split("\t", 1)[1]
        except Exception as e:
            print(f"[ERROR] {path}: {e}")
            continue

        if id_ is None and name is None:
            print(f"[WARN] no ID/Name in metaheader: {path}")
            continue

        missing = []

        if id_ and id_ not in fn:
            missing.append(f"ID='{id_}'")
        if name and name not in fn:
            missing.append(f"Name='{name}'")

        if missing:
            log.write(
                f"FILE: {path}\n"
                f"  Missing in filename: {', '.join(missing)}\n\n"
            )

log.close()
