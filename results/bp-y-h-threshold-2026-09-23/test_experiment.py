"""Physics and reproducibility checks for the standalone Bp experiment."""
import importlib
import json
from pathlib import Path
import tempfile
import tomllib
import unittest

import numpy as np
import stim

RUN = Path(__file__).resolve().parent


class SamplerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.module = importlib.import_module("sampler")
        except ModuleNotFoundError:
            cls.module = None

    def setUp(self):
        self.assertIsNotNone(self.module, "Exact operator sampler is not implemented")
        self.circuits = self.module.load_circuits(RUN / "circuits.toml")

    def test_noiseless_reference_all_distances(self):
        for d, circuit in self.circuits.items():
            with self.subTest(distance=d):
                experiment = self.module.Experiment(circuit)
                s = experiment.initial.copy(seed=11)
                branch = s.copy(seed=12)
                branch.postselect_z(experiment.n, desired_value=False)
                self.assertTrue(all(branch.peek_z(q) == 1 for q in range(experiment.n)))
                s.do(experiment.ideal)
                self.assertTrue(all(s.peek_observable_expectation(p) == 1
                                    for p in experiment.observables))
                for channel in ("y_only", "hadamard"):
                    counts = experiment.sample(channel, 0.0, 16, 23)
                    np.testing.assert_array_equal(counts, [16, 0, 0, 0])

    def test_fixed_histories_match_yao_states(self):
        fixtures = tomllib.loads((RUN / "oracle_fixtures.toml").read_text())["fixtures"]
        experiment = self.module.Experiment(self.circuits[3])
        for fixture in fixtures:
            with self.subTest(case=fixture["name"]):
                s = experiment.replay(fixture["faults"], seed=31)
                actual = s.state_vector(endian="little").astype(complex)
                expected = np.array(fixture["state_real"]) + 1j * np.array(fixture["state_imag"])
                overlap = np.vdot(expected, actual)
                self.assertAlmostEqual(abs(overlap), 1.0, places=6)
                if fixture["channel"] == "y_only":
                    measured = experiment.measure_and_decode(s)
                    self.assertEqual(measured, tuple(fixture["logical_failures"]))
                    masks = np.zeros((1, len(experiment.locations)), dtype=bool)
                    for fault in fixture["faults"]:
                        index = experiment.locations.index((fault["step"] - 1, fault["qubit"] - 1))
                        masks[0, index] ^= True
                    raw = experiment.y_measurements(masks)[0]
                    np.testing.assert_array_equal(raw, fixture["raw_measurements"])
                else:
                    events = [(f["step"] - 1, f["qubit"] - 1) for f in fixture["faults"]]
                    if len(set(events)) != len(events):
                        continue
                    masks = np.zeros((1, len(experiment.locations)), dtype=bool)
                    for event in events:
                        masks[0, experiment.locations.index(event)] = True
                    s.do(experiment.readout)
                    np.testing.assert_array_equal(
                        experiment.h_measurements(masks, [31])[0], s.current_measurement_record())

    def test_hadamard_squared_is_identity(self):
        experiment = self.module.Experiment(self.circuits[3])
        step, qubit = experiment.locations[4]
        fault = dict(step=step + 1, qubit=qubit + 1, gate="H")
        s = experiment.replay([fault, fault], seed=5)
        self.assertTrue(all(s.peek_observable_expectation(p) == 1
                            for p in experiment.observables))
        # Twirling each H separately would instead introduce a Y with probability 1/2.
        bell = stim.TableauSimulator(seed=6)
        bell.h(0)
        bell.cx(0, 1)
        exact = bell.copy(seed=7)
        exact.h(0)
        exact.h(0)
        self.assertEqual(exact.peek_observable_expectation(stim.PauliString("ZZ")), 1)
        twirl_expectations = []
        for a in ("X", "Z"):
            for b in ("X", "Z"):
                branch = bell.copy(seed=8)
                branch.do(stim.Circuit(f"{a} 0\n{b} 0"))
                twirl_expectations.append(branch.peek_observable_expectation(stim.PauliString("ZZ")))
        self.assertEqual(np.mean(twirl_expectations), 0)

    def test_probability_endpoints_and_invalid_inputs(self):
        experiment = self.module.Experiment(self.circuits[3])
        for channel in ("y_only", "hadamard"):
            counts = experiment.sample(channel, 1.0, 20, 19)
            self.assertEqual(sum(counts), 20)
            for invalid in (-0.1, 1.1, float("nan")):
                with self.assertRaises(ValueError):
                    experiment.sample(channel, invalid, 20, 19)
        with self.assertRaises(ValueError):
            experiment.sample("xz_mixture", 0.1, 20, 19)

    def test_seeded_replay_and_resume(self):
        experiment = self.module.Experiment(self.circuits[3])
        for channel in ("y_only", "hadamard"):
            np.testing.assert_array_equal(experiment.sample(channel, 0.13, 128, 77),
                                          experiment.sample(channel, 0.13, 128, 77))
            with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
                args = dict(circuit=self.circuits[3], channel=channel, p=0.13,
                            stage="test", root_seed=83, batch_size=32,
                            fingerprint="test-fixture")
                self.module.run_point(directory=Path(a), shots=64, **args)
                resumed = self.module.run_point(directory=Path(a), shots=128, **args)
                fresh = self.module.run_point(directory=Path(b), shots=128, **args)
                self.assertEqual(resumed["counts"], fresh["counts"])
                self.assertEqual(resumed["shots"], 128)
                with self.assertRaises(ValueError):
                    self.module.run_point(directory=Path(a), shots=128,
                                          **{**args, "fingerprint": "changed-physics"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
