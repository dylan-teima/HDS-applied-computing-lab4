
python3 - <<'PY'
from pathlib import Path
import hashlib

root = Path.cwd()
manifest = root / "SHA256SUMS.txt"
AI_USAGE = root / "AI_USAGE.md"

actual = hashlib.sha256(AI_USAGE.read_bytes()).hexdigest()
lines = manifest.read_text().splitlines()
updated = []
matches = 0

for line in lines:
    if not line.strip():
        updated.append(line)
        continue

    expected, filename = line.split("  ", 1)

    if filename == "AI_USAGE.md":
        print("Previous checksum:", expected)
        updated.append(f"{actual}  AI_USAGE.md")
        matches += 1
    else:
        updated.append(line)

if matches != 1:
    raise ValueError(f"Expected one AI_USAGE.md entry; found {matches}")

manifest.write_text("\n".join(updated) + "\n")

saved = next(
    line.split("  ", 1)[0]
    for line in manifest.read_text().splitlines()
    if line.endswith("  AI_USAGE.md")
)

print("Project:", root)
print("Current checksum: ", actual)
print("Saved checksum:   ", saved)
print("AI_USAGE.md check:", "PASS" if saved == actual else "FAIL")
PY
