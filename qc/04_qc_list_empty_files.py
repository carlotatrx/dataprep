import os

path = "/scratch3/PALAEO-RA/daily_data/tmp/sef_tests/qc_tests/"

empty_files = []
nonempty_files = []

for file in os.listdir(path):
    if file.endswith(".txt"):
        with open(os.path.join(path, file), "r") as f:
            lines = f.readlines()
            if len(lines) <= 1:
                empty_files.append(file)
            else:
                nonempty_files.append(file)

empty_files.sort()
nonempty_files.sort()

for file in empty_files:
    print(f"Empty file: {file}")

for file in nonempty_files:
    print(f"Non-empty file: {file}")
