import os

root = "/scratch3/PALAEO-RA/daily_data/final/"  # change this
valid_vbls = {"Tn", "Tx", "ta", "tb", "p"}

files_no_qc  = []
files_yes_qc = []
files_units_problems = []

for subdir, _, files in os.walk(root):
    for f in files:
        if not f.endswith(".tsv"):
            continue
        path = os.path.join(subdir, f)
        with open(path, "r", encoding="utf-8") as fh:
            header = []
            for line in fh:
                if line.strip().startswith("Year"):
                    break
                header.append(line.strip())

        vbl = None
        meta = ""
        units = ""
        for line in header:
            if line.startswith("Vbl"):
                vbl = line.split("\t")[1].strip()
            elif line.startswith("Meta"):
                meta = line.split("\t", 1)[1].strip() if "\t" in line else line[5:].strip()
            elif line.startswith("Units"):
                units = line.split("\t")[1].strip()

        if vbl in valid_vbls:

            if units not in {"hPa", "C"}:
                files_units_problems.append(f"{os.path.relpath(path, root)}\t{units}")
                print(f"Units problem in file: {os.path.relpath(path, root)} with units: {units}")
            else:
                if "qc" in meta.lower():
                    files_yes_qc.append(os.path.relpath(path, root))
                else:
                    files_no_qc.append(os.path.relpath(path, root))

print("Files with valid Vbl but NO QC:")
for f in files_no_qc:
    print(f)

with open('/scratch3/PALAEO-RA/daily_data/tmp/files_no_qc.txt', 'w') as f:
    f.write("\n".join(files_no_qc))

with open('/scratch3/PALAEO-RA/daily_data/tmp/files_yes_qc.txt', 'w') as f:
    f.write("\n".join(files_yes_qc))

with open('/scratch3/PALAEO-RA/daily_data/tmp/files_units_problems.txt', 'w') as f:
    f.write("\n".join(files_units_problems))
