#!/usr/bin/env python3
"""Run real Julia adapter/cache code across late-loaded backend world ages.

Only the expensive TopoNoise backend is replaced with a small Julia fixture.
The production worker, cache, settings construction, and adapter execute without
modification. No tensor-network contraction or Python protocol mock is used.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import unittest


HERE = Path(__file__).resolve().parent

# A direct access to Runner.TopoNoise or its newly defined keyword constructor
# from the worker's older world must break this test, as it did on Julia 1.12
# before the first production contraction. Load the fixture inside exercise()
# so top-level evaluation cannot silently advance the caller's world age.
JULIA_REGRESSION = r'''
using Test
include(ARGS[1])

function load_late_backend_fixture!()
    Core.eval(Runner, :(module TopoNoise
            struct Toric3DZTNSettings
                max_bond_dimension::Int
                cutoff::Float64
                negative_tolerance::Float64
            end
            Toric3DZTNSettings(; max_bond_dimension, cutoff, negative_tolerance) =
                Toric3DZTNSettings(max_bond_dimension, cutoff, negative_tolerance)
    end))
    Core.eval(Runner, :(module RepairedSeedTN
            build_decoder(context, p; settings) = (; context, p, settings)
            function decode_record(decoder, observed)
                @assert decoder.p == 0.006
                @assert decoder.settings.cutoff == 1e-18
                @assert observed.recovery_mask == 0
                return (valid=true, probabilities=[0.7, 0.1, 0.1, 0.1],
                        logical_mask=0, log_joint_syndrome=-2.0)
            end
    end))
end

function exercise_late_backend()
    load_late_backend_fixture!()
    config = (p=0.006, cutoff=1e-18)
    context = (fixture_context=true,)
    record = (seed_syndrome=Bool[], seed_recovery=Bool[],
              recovery_mask=0, final_syndrome=Bool[])

    @testset "late-loaded backend entry points" begin
      @testset "persistent decoder construction after backend load" begin
        cache = Dict{Int,Any}()
        decoder, hit, elapsed = cached_decoder!(cache, context, config, 32)
        @test hit === false
        @test elapsed >= 0.0
        @test decoder.settings.max_bond_dimension == 32
        @test decoder.settings.cutoff == 1e-18
        @test decoder.settings.negative_tolerance == 1e-10
        reused, hit, elapsed = cached_decoder!(cache, context, config, 32)
        @test hit === true
        @test reused === decoder
        @test elapsed == 0.0
        @test tn_result(decoder, record)["status"] == "complete"
      end

      @testset "one-record adapter after backend load" begin
        result = tn_result(context, record, 0.006, 64, 1e-18)
        @test result["status"] == "complete"
        @test result["valid"] === true
        @test result["logical_mask"] == 0
        @test result["probabilities"] == [0.7, 0.1, 0.1, 0.1]
      end
    end
end

exercise_late_backend()
'''


class JuliaLateLoadingTests(unittest.TestCase):
    def test_adapter_and_cache_work_after_runtime_backend_loading(self):
        julia = os.environ.get("TNML_TEST_JULIA") or shutil.which("julia")
        if not julia:
            self.skipTest("Julia is required for the real world-age regression")
        result = subprocess.run(
            [julia, "--startup-file=no", "--history-file=no", "--threads=1",
             "--depwarn=error", "-e", JULIA_REGRESSION,
             str(HERE / "persistent_worker.jl")],
            text=True, capture_output=True, timeout=60, check=False,
        )
        self.assertEqual(
            result.returncode, 0,
            "Late-loaded backend must work through the production Julia worker "
            "and adapter:\n" + result.stdout + result.stderr,
        )


if __name__ == "__main__":
    unittest.main()
