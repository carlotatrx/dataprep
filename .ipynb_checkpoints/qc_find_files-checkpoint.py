import os

"""
Walk through all subfolders under `root` and inspect every .tsv file.

For each file:
    - Read the Meta in the SEF header.
    - Check if it has the QC identifier string ("QC software=dataresqc v1.1.0").

    1. If the header contains QC but the filename does NOT end with "_qc.tsv",
       rename the file to add "_qc" before ".tsv".

    2. If the filename ends with "_qc.tsv" but the header lacks QC information,
       print a warning and rename the file to remove "_qc".

    3. Record each processed file and its QC status ("yesQC" or "noQC") in qc_overview.csv.
"""

root = "/scratch3/PALAEO-RA/daily_data/final/"
qc_string = "QC software=dataresqc"
outfile = "qc_overview.csv"

with open(outfile, "w", encoding="utf-8") as out:
    out.write("filename,qc_flag\n")

    for subdir, _, files in os.walk(root):
        for fname in files:
            if not fname.endswith(".tsv"):
                continue

            fpath = os.path.join(subdir, fname)

            # read header until data columns start
            meta_header = ""
            with open(fpath, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line.startswith("Year"):
                        break
                    if line.startswith("Meta"):
                        meta_header = line

            qc_flag = "yesQC" if qc_string in meta_header else "noQC"
            
            # rename file if QC but missing _qc in filename
            if qc_flag == "yesQC" and not fname.endswith("_qc.tsv"):
                print(f"WARNING: {fpath} has QC but doesn't ends with _qc, renaming.")

                new_fname = fname.replace(".tsv", "_qc.tsv")
                new_path = os.path.join(subdir, new_fname)
                os.rename(fpath, new_path)
                fpath = new_path
                fname = new_fname
                
            # warn if filename ends with _qc.tsv but header does NOT contain QC string
            if fname.endswith("_qc.tsv") and qc_flag != "yesQC":
                print(f"WARNING: {fpath} ends with _qc.tsv but has no QC metadata, renaming.")
                
                # remove _qc if filename has it but header lacks QC info
                new_fname = fname.replace("_qc.tsv", ".tsv")
                new_path = os.path.join(subdir, new_fname)
                os.rename(fpath, new_path)
                fpath = new_path
                fname = new_fname

            with open(outfile, "a", encoding="utf-8") as out:
                out.write(f"{fpath},{qc_flag}\n")
