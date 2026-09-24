"""Exact stochastic Clifford trajectories and correlated-Y frame sampling.

All qubit indices in exported TOML/fault fixtures are one-based. Internal
indices are zero-based. Joint count order is 00,10,01,11 (X bit first).
"""
import hashlib
import json
import math
from pathlib import Path
import time
import tomllib

import numpy as np
import pymatching
import stim


def load_circuits(path):
    records = tomllib.loads(Path(path).read_text())["circuits"]
    return {c["distance"]: c for c in records}


def stable_seed(*parts):
    key = json.dumps(parts, separators=(",", ":"), ensure_ascii=True).encode()
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "little")


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


class Experiment:
    def __init__(self, circuit):
        self.circuit = circuit
        self.d = circuit["distance"]
        self.n = self.d ** 2
        assert circuit["construction"] == "bp" and circuit["clock"] == "gate_layer"
        assert circuit["input_state"] == "zero" and circuit["boundary_orientation"] == "x_ns"
        self.layers = circuit["layers"]
        self.locations = [(i, q - 1) for i, layer in enumerate(self.layers)
                          for q in layer["active_qubits"]]
        self.location_steps = np.array([i for i, _ in self.locations])
        self.location_qubits = np.array([q for _, q in self.locations])
        self.a = [np.array(s) - 1 for s in circuit["a_s"]]
        self.b = [np.array(s) - 1 for s in circuit["b_p"]]
        self.lx = np.array(circuit["logical_x"]) - 1
        self.lz = np.array(circuit["logical_z"]) - 1
        self.layer_text = []
        for layer in self.layers:
            lines = []
            for operation in layer["operations"]:
                gate = "CX" if operation["gate"] == "CNOT" else operation["gate"]
                if gate not in ("H", "CX", "X", "Z"):
                    raise ValueError(f"Unsupported ideal gate {gate}")
                lines.append(gate + " " + " ".join(str(q - 1) for q in operation["qubits"]) + "\n")
            self.layer_text.append("".join(lines))
        self.ideal_text = "".join(self.layer_text)
        self.offsets = np.cumsum([0] + [len(s) for s in self.layer_text]).tolist()
        self.ideal = stim.Circuit(self.ideal_text)
        self.observables = ([self.pauli(s, "Z") for s in self.a]
                            + [self.pauli(s, "X") for s in self.b]
                            + [self.pauli(self.lz, "Z", reference=True),
                               self.pauli(self.lx, "X", reference=True)])
        self.readout = stim.Circuit()
        for observable in self.observables:
            targets = []
            for qubit, kind in enumerate(observable):
                if kind:
                    if targets:
                        targets.append(stim.target_combiner())
                    targets.append({1: stim.target_x, 2: stim.target_y, 3: stim.target_z}[kind](qubit))
            self.readout.append("MPP", targets)
        # This is exactly (U_Bp† ⊗ I_R)|Phi_L+>, prepared without noise.
        self.initial = stim.TableauSimulator(seed=0)
        self.initial.set_state_from_stabilizers(self.observables)
        self.initial.do(self.ideal.inverse())
        self.matcher_x = self.matcher(self.a, self.lz)
        self.matcher_z = self.matcher(self.b, self.lx)

    def pauli(self, support, kind, reference=False):
        pauli = stim.PauliString(self.n + 1)
        for q in support:
            pauli[int(q)] = kind
        if reference:
            pauli[self.n] = kind
        return pauli

    def matcher(self, checks, logical_support):
        matrix = np.zeros((len(checks), self.n), dtype=np.uint8)
        for row, support in enumerate(checks):
            matrix[row, support] = 1
        assert np.all(matrix.sum(axis=0) <= 2)
        observable = np.zeros((1, self.n), dtype=np.uint8)
        observable[0, logical_support] = 1
        return pymatching.Matching.from_check_matrix(matrix, weights=1.0, faults_matrix=observable)

    def replay(self, faults, seed):
        """Replay an explicit operator history, retaining duplicate events in fixtures."""
        after = {}
        for fault in faults:
            step, qubit = int(fault["step"]) - 1, int(fault["qubit"]) - 1
            if (step, qubit) not in self.locations or fault["gate"] not in ("H", "Y"):
                raise ValueError("Fault is not an eligible H/Y event")
            after.setdefault(step, []).append(f"{fault['gate']} {qubit}\n")
        text = "".join(gates + "".join(after.get(i, [])) for i, gates in enumerate(self.layer_text))
        simulator = self.initial.copy(seed=seed)
        simulator.do(stim.Circuit(text))
        return simulator

    def decode(self, measurements):
        measurements = np.asarray(measurements, dtype=np.uint8)
        a_count = len(self.a)
        b_count = len(self.b)
        px = self.matcher_x.decode_batch(np.ascontiguousarray(measurements[:, :a_count]),
                                         enable_correlations=False).reshape(-1)
        pz = self.matcher_z.decode_batch(np.ascontiguousarray(measurements[:, a_count:a_count+b_count]),
                                         enable_correlations=False).reshape(-1)
        return np.column_stack((px ^ measurements[:, -2], pz ^ measurements[:, -1]))

    def measure_and_decode(self, simulator):
        simulator.do(self.readout)
        bits = simulator.current_measurement_record()[-len(self.observables):]
        return tuple(int(v) for v in self.decode([bits])[0])

    def y_measurements(self, masks):
        """One mask bit injects both frame components at the same location."""
        masks = np.asarray(masks, dtype=bool)
        x = np.zeros((len(masks), self.n), dtype=bool)
        z = np.zeros_like(x)
        offset = 0
        for layer in self.layers:
            for operation in layer["operations"]:
                gate = operation["gate"]
                qubits = [q - 1 for q in operation["qubits"]]
                if gate == "H":
                    q = qubits[0]
                    x[:, q], z[:, q] = z[:, q].copy(), x[:, q].copy()
                elif gate == "CNOT":
                    control, target = qubits
                    x[:, target] ^= x[:, control]
                    z[:, control] ^= z[:, target]
            support = [q - 1 for q in layer["active_qubits"]]
            events = masks[:, offset:offset + len(support)]
            x[:, support] ^= events
            z[:, support] ^= events
            offset += len(support)
        measurements = [np.logical_xor.reduce(x[:, s], axis=1) for s in self.a]
        measurements += [np.logical_xor.reduce(z[:, s], axis=1) for s in self.b]
        measurements += [np.logical_xor.reduce(x[:, self.lz], axis=1),
                         np.logical_xor.reduce(z[:, self.lx], axis=1)]
        return np.column_stack(measurements).astype(np.uint8)

    def h_measurements(self, masks, measurement_seeds):
        measurements = np.empty((len(masks), len(self.observables)), dtype=np.uint8)
        for shot, mask in enumerate(masks):
            locations = np.flatnonzero(mask)
            parts = []
            previous = 0
            i = 0
            while i < len(locations):
                step = int(self.location_steps[locations[i]])
                end = i + 1
                while end < len(locations) and self.location_steps[locations[end]] == step:
                    end += 1
                parts.append(self.ideal_text[self.offsets[previous]:self.offsets[step + 1]])
                parts.append("H " + " ".join(str(self.location_qubits[k]) for k in locations[i:end]) + "\n")
                previous = step + 1
                i = end
            parts.append(self.ideal_text[self.offsets[previous]:])
            noisy = stim.Circuit("".join(parts))
            noisy += self.readout
            simulator = self.initial.copy(seed=int(measurement_seeds[shot]))
            simulator.do(noisy)
            measurements[shot] = simulator.current_measurement_record()
        return measurements

    def sample(self, channel, p, shots, seed):
        if channel not in ("y_only", "hadamard"):
            raise ValueError("channel must be y_only or hadamard")
        if not math.isfinite(p) or not 0 <= p <= 1:
            raise ValueError("p must be a finite total fault probability in [0,1]")
        if not isinstance(shots, int) or shots <= 0:
            raise ValueError("shots must be a positive integer")
        rng = np.random.default_rng(seed)
        masks = rng.random((shots, len(self.locations))) < p
        if channel == "y_only":
            measurements = self.y_measurements(masks)
        else:
            seeds = rng.integers(0, 2**64, size=shots, dtype=np.uint64)
            measurements = self.h_measurements(masks, seeds)
        failures = self.decode(measurements)
        return np.bincount(failures[:, 0] + 2 * failures[:, 1], minlength=4).astype(np.int64)


def run_point(*, circuit, channel, p, stage, shots, root_seed, batch_size,
              fingerprint, directory):
    """Atomic per-point checkpoints; fixed chunk seeds make worker order irrelevant."""
    if shots <= 0 or batch_size <= 0 or shots % batch_size:
        raise ValueError("shots must be a positive multiple of the fixed batch size")
    distance = circuit["distance"]
    p_key = int(round(p * 1_000_000))
    if abs(p - p_key / 1_000_000) > 1e-12:
        raise ValueError("p must be represented exactly on the 1e-6 grid")
    identity = dict(channel=channel, distance=distance, p=p_key / 1_000_000,
                    stage=stage, root_seed=root_seed, batch_size=batch_size,
                    fingerprint=fingerprint)
    path = Path(directory) / f"{stage}-{channel}-d{distance}-p{p_key:06d}.json"
    if path.exists():
        record = json.loads(path.read_text())
        if any(record[k] != v for k, v in identity.items()):
            raise ValueError(f"Checkpoint identity mismatch: {path}")
        if sum(record["counts"]) != record["shots"] or record["shots"] % batch_size:
            raise ValueError(f"Invalid checkpoint counts: {path}")
    else:
        record = {**identity, "shots": 0, "counts": [0, 0, 0, 0], "seconds": 0.0}
    if record["shots"] > shots:
        raise ValueError("Requested target is below the saved shot count")
    if record["shots"] == shots:
        return record
    experiment = Experiment(circuit)
    while record["shots"] < shots:
        start = time.perf_counter()
        seed = stable_seed(root_seed, stage, channel, distance, p_key, record["shots"])
        counts = experiment.sample(channel, p, batch_size, seed)
        record["counts"] = (np.array(record["counts"]) + counts).tolist()
        record["shots"] += batch_size
        record["seconds"] += time.perf_counter() - start
        atomic_json(path, record)
    return record
