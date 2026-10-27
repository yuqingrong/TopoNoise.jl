# Run the TN scan on HKUST(GZ) HPC4

This package runs the existing repaired-seed, Z-error, diagonal H+6CNOT preparation experiment. It includes the initial seed repair and applies no later feedback. The TN decoder conditions on the seed and final syndrome records and the known repair. It uses the existing Julia TN implementation; it does not change the circuit or the noise model.

The confirmed login architecture is `aarch64`, and `julia/1.12.4` is available. The [Guangzhou Phase-IV guide](https://docs.hpc.hkust-gz.edu.cn/docs/hpc4/domestic/quickstart-submit-slurm-job/) identifies the CPU partition as `hpc`. System Python 3.9.9 cannot run these controllers, which require Python 3.11 or 3.12. Setup creates a separate native Python environment through CondaPkg. Do not copy the Mac `.CondaPkg` directory or Julia depot.

HPC can increase throughput by running independent shot batches concurrently. Each contraction uses one Julia thread and one BLAS thread. Allocating many CPUs to one task will not parallelize this implementation. Single-core speed on Kunpeng, ARM dependency installation, account limits, and actual memory use still need to be measured on the cluster. The resource requests below are starting values, not verified site limits.

## 1. Transfer a fresh source copy

The supplied transfer archive contains the current scientific sources, dependency files, experiment scripts, and this HPC directory. It excludes environments, credentials, Git metadata, and previous results. It preserves the working source files rather than assuming they are all in a Git commit.

On your Mac, upload the supplied archive:

```bash
scp /Users/rongyuqing/jcode/TopoNoise.jl/results/tnml-hpc4-transfer-2026-10-27/TopoNoise-hpc4-20261027.tar.gz yrong265@hpc4login.hpc.hkust-gz.edu.cn:~/
```

In the HPC terminal, create a **new** directory and extract there. `mkdir` deliberately fails if that directory already exists. You may use a persistent shared project/storage directory instead of home; keep its absolute path unchanged after setup.

```bash
mkdir "$HOME/tnml-hpc-20261027"
tar -xzf "$HOME/TopoNoise-hpc4-20261027.tar.gz" -C "$HOME/tnml-hpc-20261027"
cd "$HOME/tnml-hpc-20261027"
```

Always submit the jobs below from this repository directory. Slurm copies the batch script to a spool directory, so the scripts locate the checkout through `SLURM_SUBMIT_DIR`.

## 2. Build and check the native environment

```bash
sbatch examples/tnml_repaired_seed/hpc/setup.sbatch
squeue -u "$USER"
```

Wait for this job to finish successfully. Its logs are `slurm-tnml-setup-JOBID.out` and `.err`; its detailed report is under `hpc_environment/setup-JOBID/`. Success publishes `.hpc-python-path`. Setup uses `julia/1.12.4`, creates `.hpc-julia-depot`, resolves the native CondaPkg environment, imports TopoNoise and Python extensions, and checks the native PyMatching binding. It also saves the package versions and a Python dependency listing.

The source manifest was generated with Julia 1.12.6; this cluster setup checks the available 1.12.4 module explicitly. Freeze the study only after setup, so it records the actual installed dependency files. PyMatching or other dependencies may need a native compiler build on ARM. If setup fails, inspect the logs before submitting any simulations. Network access, native build tools, and the two-hour/eight-GiB setup allocation are not yet verified on HPC4.

After setup succeeds, load the environment in your terminal and create a fresh study manifest. These are lightweight controller commands; numerical work is submitted to Slurm.

```bash
module load julia/1.12.4
source examples/tnml_repaired_seed/hpc/env.sh
export TNML_STUDY="$HOME/tnml-hpc-study-20261027"
"$TNML_PYTHON" examples/tnml_repaired_seed/hpc/manifest.py create --root "$TNML_STUDY"
```

Keep this checkout and environment fixed for the study. Every task checks frozen scientific and controller file hashes. Changing sources or runtime paths requires a separate study, not resuming old outputs. Cluster seeds are separate from the Mac study; keep their results separate.

## 3. Submit compatibility and exact-comparison checks

```bash
sbatch --array=0-4%4 examples/tnml_repaired_seed/hpc/array.sbatch zero "$TNML_STUDY"
sbatch --array=0-2%3 examples/tnml_repaired_seed/hpc/array.sbatch smoke "$TNML_STUDY"
```

- `zero`: one zero-noise diagnostic for each of L=3,5,7,9,11. These are diagnostics, not 1,000 simulated shots. The analytical failure probability at p=0 is zero, without a Monte Carlo error bar.
- `smoke`: two new histories at each of p=0.6%, 1%, and 2% for L=3, with an independent exact ML comparison for every history.

After the jobs finish, collect their status:

```bash
sbatch examples/tnml_repaired_seed/hpc/collect.sbatch "$TNML_STUDY"
```

Read the collector output and the study's `collected/` directory. Small smoke checks establish that the cluster environment runs the model; they do not establish a threshold or prove convergence across all noise rates. Production continues to compare every L=3 shot with exact ML. Production dispatch refuses to start unless all smoke tasks pass their audit.

## 4. Run L=3 production in parallel

The frozen plan uses χ=240, cutoff=1e-26, and negative tolerance=1e-10. Each nonzero-p point has **20 independent batches of 50 shots**, totaling **1,000 shots**. The 22 positive rates, expressed here as percentages, are:

```text
0.25, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75,
3.00, 3.25, 3.50, 3.75, 4.00, 4.25
```

There are 440 production tasks and 22,000 production histories. Task IDs 0–19 cover the first probability; 20–39 cover the next, and so on. Start with one point at four concurrent tasks:

```bash
sbatch --array=0-19%4 examples/tnml_repaired_seed/hpc/array.sbatch production "$TNML_STUDY"
```

After confirming its completion, exact audits, memory use, and runtime, submit the remaining tasks if your account's submission limits permit it:

```bash
sbatch --array=20-439%4 examples/tnml_repaired_seed/hpc/array.sbatch production "$TNML_STUDY"
```

`%4` caps concurrent tasks in that array. It is adjustable to your allocation. Multiple submitted arrays each have their own cap. If the site limits array size or queued tasks, submit smaller disjoint ranges (for example `20-39%4`). The [Slurm array documentation](https://slurm.schedmd.com/job_array.html) explains the syntax. No array is submitted automatically by this package.

Default array resources are one CPU, 10 GiB requested memory, and eight hours. The production supervisor allows 8 GiB worker RSS, 600 seconds per contraction, and six hours per task. Scheduler time and memory should exceed the supervisor limits to leave room for orderly checkpointing and Python overhead. Limits are frozen in the study manifest; raising them requires a deliberate new configuration/study, not silently changing a running experiment.

## 5. Measure larger-L feasibility separately

```bash
sbatch --array=0-5%2 --mem=18G --time=03:00:00 examples/tnml_repaired_seed/hpc/array.sbatch feasibility "$TNML_STUDY"
```

This runs L=5 at p=0.6%, 1%, and 2%, each with three histories decoded at χ=240 and χ=480. Each pair uses the same RNG seed, and the collector verifies the saved histories match. These are convergence/runtime diagnostics, not production samples. They use a 16 GiB worker RSS limit, 1,800 seconds per contraction, and two hours per task. Agreement on three histories is only an initial check; disagreement or a timeout means the larger scan is not yet ready.

The requested L=3,5,7,9,11 figure remains the eventual target. This package enables the full L=3 probability range and starts larger-size diagnostics. It deliberately does not launch thousands of unvalidated L=5–11 contractions. HPC parallelism reduces the time needed for independent work but does not remove growth in the cost of a single tensor contraction.

## 6. Gather, inspect, and resume

```bash
sbatch examples/tnml_repaired_seed/hpc/collect.sbatch "$TNML_STUDY"
squeue -u "$USER"
```

For a finished array, replace `JOBID` with its actual number:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

The collector audits every **planned** batch, including batches that never started. Missing, invalid, and unaudited histories remain unknown. It reports a failure-rate point and Wilson interval only when all 1,000 planned histories qualify. Diagnostics are kept separate from production. It compares saved TN posteriors/decisions with separately saved exact L3 results; it never substitutes an exact decision for the TN decision.

Each task has its own output directory, lock, status, raw corpus, and checkpoints. To resume an interrupted task, resubmit only its stage and task ID with the same source/runtime paths and study manifest. For example:

```bash
sbatch --array=7 examples/tnml_repaired_seed/hpc/array.sbatch production "$TNML_STUDY"
```

An invalid posterior or exact-ML disagreement needs investigation; repeatedly resubmitting does not fix it. Task locks use operating-system `flock`: they release when their owning processes exit. The empty `.task.lock` file may remain and should be left in place; deleting it can allow duplicate workers to run. The setup job uses a separate `.hpc-setup.lock` directory, which may need manual removal after a forced kill, but only after checking that the setup job has ended.

To download the study, run on your Mac (use the actual storage path if you changed it):

```bash
scp -r yrong265@hpc4login.hpc.hkust-gz.edu.cn:~/tnml-hpc-study-20261027 /Users/rongyuqing/jcode/TopoNoise.jl/results/
```

## Local verification and remaining cluster checks

Python control tests and shell syntax can be checked locally without running contractions:

```bash
python3 -m unittest discover -s examples/tnml_repaired_seed/hpc -p 'test_*.py'
bash -n examples/tnml_repaired_seed/hpc/env.sh
bash -n examples/tnml_repaired_seed/hpc/setup.sbatch
bash -n examples/tnml_repaired_seed/hpc/array.sbatch
bash -n examples/tnml_repaired_seed/hpc/collect.sbatch
```

The test command needs Python >=3.11. Preparing and testing these files locally does not verify ARM package installation, scheduler acceptance, compute-node networking, convergence, or performance. The setup and diagnostic stages provide those checks before the larger allocation is used.
