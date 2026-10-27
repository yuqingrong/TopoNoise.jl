#!/usr/bin/env julia
# Run through setup.sbatch, or after sourcing env.sh in a compute allocation.
# This check imports native dependencies; it does not run a TN experiment.
import Dates, LinearAlgebra, SHA, TOML

length(ARGS) == 1 || error("Usage: environment_check.jl NEW_REPORT_DIRECTORY")
haskey(ENV, "TNML_REPO") || error("Source hpc/env.sh first")
const TNML_ROOT = realpath(ENV["TNML_REPO"])
const TNML_REPORT_DIR = abspath(only(ARGS))
mkpath(TNML_REPORT_DIR)
const TNML_REPORT_PATH = joinpath(TNML_REPORT_DIR, "environment.toml")
isfile(TNML_REPORT_PATH) && error("Refusing to overwrite an existing environment report")
const TNML_SETUP = get(ENV, "TNML_HPC_SETUP", "0") == "1"
const TNML_REPORT = Dict{String,Any}(
    "valid" => false,
    "stage" => "preflight",
    "checked_at_utc" => string(Dates.now(Dates.UTC)),
    "repository" => TNML_ROOT,
    "mode" => TNML_SETUP ? "setup" : "runtime",
    "julia_version" => string(VERSION),
    "julia_architecture" => string(Sys.ARCH),
    "kernel" => string(Sys.KERNEL),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", ""),
)

function tnml_atomic_text(path::String, text::String)
    temporary, io = mktemp(dirname(path))
    try
        write(io, text)
        close(io)
        mv(temporary, path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

function tnml_write_report()
    io = IOBuffer()
    TOML.print(io, TNML_REPORT; sorted=true)
    tnml_atomic_text(TNML_REPORT_PATH, String(take!(io)))
end

# Failed checks retain a valid=false report and their stage; check.log contains
# the exception. The Python path is published only after every check succeeds.
atexit() do
    isfile(TNML_REPORT_PATH) || tnml_write_report()
end

Sys.islinux() && Sys.ARCH == :aarch64 || error("Expected native Linux aarch64 Julia")
VERSION == v"1.12.4" || error("Expected the julia/1.12.4 module, got $VERSION")
realpath(dirname(Base.active_project())) == TNML_ROOT || error("Wrong Julia project")
Threads.nthreads() == 1 || error("Julia must have one worker thread")
LinearAlgebra.BLAS.get_num_threads() == 1 || error("BLAS must have one thread")
get(ENV, "JULIA_CPU_TARGET", "") == "generic" || error("JULIA_CPU_TARGET must be generic")
DEPOT_PATH == [joinpath(TNML_ROOT, ".hpc-julia-depot")] || error("Wrong shared Julia depot")
TNML_REPORT["julia_threads"] = Threads.nthreads()
TNML_REPORT["blas_threads"] = LinearAlgebra.BLAS.get_num_threads()
TNML_REPORT["julia_cpu_target"] = ENV["JULIA_CPU_TARGET"]
TNML_REPORT["depot_path"] = copy(DEPOT_PATH)
TNML_REPORT["build_tools"] = Dict(tool => something(Sys.which(tool), "unavailable")
    for tool in ("c++", "g++", "gcc", "cc", "cmake", "make"))

TNML_REPORT["stage"] = "native Python selection"
import CondaPkg
const TNML_BACKEND = CondaPkg.backend()
const TNML_PYTHON = if TNML_SETUP
    TNML_BACKEND in (:Pixi, :MicroMamba) || error("Expected managed native CondaPkg backend")
    something(CondaPkg.which("python"), "")
else
    TNML_BACKEND == :Null || error("Runtime must disable CondaPkg resolution")
    get(ENV, "TNML_PYTHON", "")
end
isabspath(TNML_PYTHON) && isfile(TNML_PYTHON) || error("No absolute native Python executable")
startswith(realpath(TNML_PYTHON), TNML_ROOT * "/") || error("Python is outside this repository")
ENV["JULIA_PYTHONCALL_EXE"] = TNML_PYTHON
TNML_REPORT["condapkg_backend"] = string(TNML_BACKEND)
TNML_REPORT["resolved_python"] = TNML_PYTHON

TNML_REPORT["stage"] = "PythonCall and native extension imports"
import PythonCall, TopoNoise
const TNML_PYSYS = PythonCall.pyimport("sys")
const TNML_PYPLATFORM = PythonCall.pyimport("platform")
const TNML_PYVERSION = VersionNumber(
    PythonCall.pyconvert(Int, TNML_PYSYS.version_info.major),
    PythonCall.pyconvert(Int, TNML_PYSYS.version_info.minor),
    PythonCall.pyconvert(Int, TNML_PYSYS.version_info.micro),
)
const TNML_PYARCH = PythonCall.pyconvert(String, TNML_PYPLATFORM.machine())
const TNML_ACTUAL_PYTHON = PythonCall.pyconvert(String, TNML_PYSYS.executable)
TNML_REPORT["python_version"] = string(TNML_PYVERSION)
TNML_REPORT["python_architecture"] = TNML_PYARCH
TNML_REPORT["pythoncall_executable"] = TNML_ACTUAL_PYTHON
v"3.11" <= TNML_PYVERSION < v"3.13" || error(
    "The controller requires Python >=3.11,<3.13; CondaPkg selected $TNML_PYVERSION. " *
    "Stop here and prepare a fresh HPC copy with Python 3.12 before freezing a study.")
lowercase(TNML_PYARCH) == "aarch64" || error("Python is not native aarch64")
realpath(TNML_ACTUAL_PYTHON) == realpath(TNML_PYTHON) || error("PythonCall selected a different Python")
realpath(pathof(TopoNoise)) == realpath(joinpath(TNML_ROOT, "src", "TopoNoise.jl")) || error("Loaded another TopoNoise checkout")

const TNML_PM = PythonCall.pyimport("pymatching")
const TNML_PM_NATIVE = PythonCall.pyimport("pymatching._cpp_pymatching")
const TNML_NUMPY = PythonCall.pyimport("numpy")
const TNML_SCIPY = PythonCall.pyimport("scipy")
const TNML_LDPC = PythonCall.pyimport("ldpc")
const TNML_METADATA = PythonCall.pyimport("importlib.metadata")
TNML_REPORT["python_packages"] = Dict(name => PythonCall.pyconvert(String, TNML_METADATA.version(name))
    for name in ("pymatching", "numpy", "scipy", "ldpc", "pip"))
TNML_REPORT["python_packages"]["pymatching"] == "2.4.0" || error("Expected pymatching==2.4.0")
TNML_REPORT["python_packages"]["ldpc"] == "2.4.1" || error("Expected ldpc==2.4.1")
TNML_REPORT["pymatching_native_extension"] = PythonCall.pyconvert(String, TNML_PM_NATIVE.__file__)

# Exercise the native decoder binding without invoking any circuit simulation.
TNML_REPORT["stage"] = "native PyMatching smoke check"
const TNML_MATCHER = TNML_PM.Matching()
TNML_MATCHER.add_edge(0, 1; fault_ids=0, weight=1.0)
const TNML_DECODED = TNML_MATCHER.decode(TNML_NUMPY.array([1, 1]; dtype="uint8"))
PythonCall.pyconvert(Vector{Int}, TNML_DECODED) == [1] || error("Native PyMatching smoke check failed")

TNML_REPORT["stage"] = "dependency record"
const TNML_FREEZE = read(`$TNML_PYTHON -m pip freeze`, String)
tnml_atomic_text(joinpath(TNML_REPORT_DIR, "pip-freeze.txt"), TNML_FREEZE)
TNML_REPORT["pip_freeze_sha256"] = bytes2hex(SHA.sha256(TNML_FREEZE))
TNML_REPORT["file_sha256"] = Dict(name => bytes2hex(SHA.sha256(read(joinpath(TNML_ROOT, name))))
    for name in ("Project.toml", "Manifest.toml", "CondaPkg.toml"))
TNML_REPORT["valid"] = true
TNML_REPORT["stage"] = "complete"
tnml_write_report()
if TNML_SETUP
    ispath(joinpath(TNML_ROOT, ".hpc-python-path")) && error("A validated Python path already exists")
    tnml_atomic_text(joinpath(TNML_ROOT, ".hpc-python-path"), TNML_PYTHON * "\n")
end
println("Julia: ", VERSION, " (", Sys.ARCH, ")")
println("PythonCall: ", TNML_ACTUAL_PYTHON, " — Python ", TNML_PYVERSION, " (", TNML_PYARCH, ")")
println("PyMatching: ", TNML_REPORT["python_packages"]["pymatching"])
println("Julia / BLAS threads: ", Threads.nthreads(), " / ", LinearAlgebra.BLAS.get_num_threads())
println("Environment check passed: ", TNML_REPORT_PATH)
