import re, os

REPO = "/mnt/c/Users/USER/Desktop/schizochytrium-hs6-genomics"
html = open(f"{REPO}/lysis-enzymes.html", encoding="utf-8").read()

refs = set(re.findall(r'(?:src|href)="(figures/lysis_enzymes/[^"]+)"', html))
missing = []
for r in sorted(refs):
    path = f"{REPO}/{r}"
    if not os.path.exists(path):
        missing.append(r)
    elif os.path.getsize(path) < 5000:
        missing.append(r + " (suspiciously small)")

print(f"{len(refs)} unique figure references found")
if missing:
    print("MISSING/SUSPECT:")
    for m in missing:
        print(" ", m)
else:
    print("All references resolve to real, non-trivial files.")

# also check for leftover placeholder text
placeholders = re.findall(r'\[.*?inserted here.*?\]', html)
if placeholders:
    print("LEFTOVER PLACEHOLDERS:", placeholders)
else:
    print("No leftover placeholder text.")
