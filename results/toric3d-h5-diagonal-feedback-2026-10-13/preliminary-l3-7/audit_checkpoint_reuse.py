"""Independent, read-only verification of both user-directed checkpoint subsets."""
import hashlib
import json
import tomllib
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parent
repo = root.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return tomllib.loads(path.read_text())


original_path = root / "production/checkpoints/checkpoint.toml"
parent_path = root / "preliminary-l3-9/production/checkpoints/checkpoint.toml"
final_path = here / "production/checkpoints/checkpoint.toml"
assert sha(original_path) == "772da74e8dd161de97c00c5008773fb19982ecb55e4a3343489309656323d271"
assert sha(parent_path) == "c3544f0d3f72ff2cd5fe0e79379ec502b0360cdc15e76a3f7b9cff4a18faa80c"
for subdir, expected in (
    (root, "5041dd833aa8f28f318b4f5fb495f5025d61dfeb87f1a8ab8d76e5edcf6d7e00"),
    (root / "preliminary-l3-9", "dca614dc243ec8d918cb91ac27c74425c7af0fbba6f2381c88d33d3ebf025412"),
    (here, "1597b9941132d3f24c11c68ec1933ea639f5289225fe80ff5070367e81416e72"),
):
    assert sha(subdir / "selection.json") == expected

original, parent, final = map(read, (original_path, parent_path, final_path))
manifest = read(here / "production/checkpoint-reuse.toml")
indices = manifest["parent_indices_1based"]
assert len(indices) == len(final["states"]) == 195
assert len(parent["states"]) == 260 and len(original["states"]) == 390
assert manifest["source_checkpoint_sha256"] == sha(parent_path)
assert final["point_seeds"] == [parent["point_seeds"][i - 1] for i in indices]
# The four-size precursor retained the completed X prefix of the six-size run.
assert parent["point_seeds"][:80] == original["point_seeds"][:80]
assert parent["states"][:80] == original["states"][:80]
reused_points = 0
for state, index in zip(final["states"], indices):
    before = parent["states"][index - 1]
    assert state["completed_shots"] == 1000 and state["invalid_recoveries"] == 0
    if before["completed_shots"] == 1000:
        # Includes all 23 joint arrays, timings, record and signature digests.
        assert state == before
        reused_points += 1
    else:
        assert before["completed_shots"] == 0
assert reused_points == 186 and manifest["reused_paired_records"] == 186000
for checkpoint in (original, parent, final):
    assert checkpoint["config"]["seed"] == 1234
    assert checkpoint["config"]["shots"] == 1000
    assert sum(s["invalid_recoveries"] for s in checkpoint["states"]) == 0
    for filename, expected in checkpoint["config"]["source_hashes"].items():
        assert sha(repo / filename) == expected
for field in ("decoder_settings", "source_hashes", "backend_versions", "protocols"):
    assert original["config"][field] == parent["config"][field] == final["config"][field]
for channel in ("x_only", "z_only"):
    assert final["config"]["channel_grids"][channel]["dimensions"] == [[3, 3, 3], [5, 5, 5], [7, 7, 7]]
    assert original["config"]["channel_grids"][channel]["rates"] == parent["config"]["channel_grids"][channel]["rates"] == final["config"]["channel_grids"][channel]["rates"]

report = dict(
    passed=True, complete_points=195, paired_records=195000,
    exact_reused_points=reused_points, exact_reused_records=186000,
    newly_completed_records=9000, inherited_point_seeds_verified=True,
    original_to_four_size_exact_X_points=80,
    original_preserved_records=sum(s["completed_shots"] for s in original["states"]),
    four_size_preserved_records=sum(s["completed_shots"] for s in parent["states"]),
    original_checkpoint_sha256=sha(original_path),
    parent_checkpoint_sha256=sha(parent_path), final_checkpoint_sha256=sha(final_path),
    selection_sha256=sha(here / "selection.json"), audit_script_sha256=sha(Path(__file__)),
    source_hashes_verified=True, invalid_recoveries=0,
    seed_rule=manifest["seed_rule"],
)
(here / "production/independent-reuse-audit.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
