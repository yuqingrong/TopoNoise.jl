# Runnable Toric-Code Circuit Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a command-line Julia example that constructs the executable toric-code Yao circuit and exports native SVG, PDF, and PNG drawings.

**Architecture:** A self-contained module in `examples/generate_toric_circuit.jl` exposes a testable `main(args; output_dir, io, error_io)` function. It uses the package's public PEPS, schedule, and `yao_circuit` APIs, never allocates a state vector, and calls `Yao.vizcircuit` directly for each output format.

**Tech Stack:** Julia 1.10+, TopoNoise.jl, Yao 0.9, Julia `Test`.

## Global Constraints

- Use the repository root `Project.toml`; do not create another Julia environment.
- Default to a `2x2` lattice when no dimensions are supplied.
- Accept either zero or two command-line arguments; both dimensions must be positive integers.
- Write generated files below `examples/output/` by default and keep that directory ignored by Git.
- Do not simulate a register or materialize the global circuit matrix.
- Keep the existing Yao and CairoMakie APIs unchanged.

---

### Task 1: Test the runnable example contract

**Files:**
- Create: `test/example.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `toric_code_peps`, `sequential_circuit_graph`, `circuit_layers`, `yao_circuit`, and `Yao.vizcircuit`.
- Produces test expectations for `GenerateToricCircuitExample.main(args=String[]; output_dir, io, error_io)::Int` and `GenerateToricCircuitExample.generate(rows, cols; output_dir, io)`.

- [ ] **Step 1: Write the failing example tests**

Create `test/example.jl`. Assert that the requested example file exists before
including it, so the RED phase is a deliberate failed assertion rather than a
load error:

```julia
example_path = joinpath(
    @__DIR__, "..", "examples", "generate_toric_circuit.jl")

@testset "Runnable toric-code circuit example" begin
    @test isfile(example_path)

    if isfile(example_path)
        include(example_path)

        mktempdir() do directory
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                String[]; output_dir=directory, io=output, error_io=output) == 0

            text = String(take!(output))
            @test contains(text, "lattice: 2x2")
            @test contains(text, "qubits: 20")
            @test contains(text, "layer 1: (2,1)")
            @test contains(text, "layer 2: (2,2), (1,1)")
            @test contains(text, "layer 3: (1,2)")

            for extension in ("svg", "pdf", "png")
                path = joinpath(directory, "toric_circuit_2x2.$extension")
                @test isfile(path)
                @test filesize(path) > 100
            end
            @test startswith(
                read(joinpath(directory, "toric_circuit_2x2.svg"), String), "<?xml")
            @test startswith(
                read(joinpath(directory, "toric_circuit_2x2.pdf"), String), "%PDF")
        end

        mktempdir() do directory
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                ["1", "1"]; output_dir=directory, io=output, error_io=output) == 0
            @test contains(String(take!(output)), "qubits: 6")
            @test isfile(joinpath(directory, "toric_circuit_1x1.svg"))
        end

        for arguments in (["2"], ["two", "2"], ["0", "2"], ["2", "-1"])
            output = IOBuffer()
            @test GenerateToricCircuitExample.main(
                arguments; output_dir="unused", io=output, error_io=output) == 1
            text = String(take!(output))
            @test contains(text, "Usage:")
            @test contains(text, "positive integers")
        end
    end
end
```

Include it from `test/runtests.jl` after `sequential_circuit.jl`:

```julia
include("example.jl")
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
julia --project=. -e 'using TopoNoise, ITensors, Test; include("test/example.jl")'
```

Expected: one failed `isfile(example_path)` assertion because
`examples/generate_toric_circuit.jl` does not exist.

- [ ] **Step 3: Commit the failing tests only if requested**

Do not stage or commit by default because the checkout already contains user-owned uncommitted work. Preserve the RED output as the TDD evidence.

---

### Task 2: Implement the command-line generator

**Files:**
- Create: `examples/generate_toric_circuit.jl`

**Interfaces:**
- Consumes: the public TopoNoise APIs and `Yao.nqubits`/`Yao.vizcircuit`.
- Produces: `GenerateToricCircuitExample.generate(rows::Integer, cols::Integer; output_dir::AbstractString, io::IO=stdout)` and `GenerateToricCircuitExample.main(args=ARGS; output_dir=joinpath(@__DIR__, "output"), io::IO=stdout, error_io::IO=stderr)::Int`.

- [ ] **Step 1: Add the minimal example module**

Create the script with a module so tests can include it without triggering CLI exit:

```julia
module GenerateToricCircuitExample

using TopoNoise
import Yao

const USAGE = "Usage: julia --project=. examples/generate_toric_circuit.jl [rows cols] (rows and cols must be positive integers)"

function _dimensions(args)
    isempty(args) && return (2, 2)
    length(args) == 2 || throw(ArgumentError("expected zero or two dimensions"))
    dimensions = try
        parse.(Int, args)
    catch
        throw(ArgumentError("rows and cols must be positive integers"))
    end
    all(>(0), dimensions) ||
        throw(ArgumentError("rows and cols must be positive integers"))
    return Tuple(dimensions)
end

function _site_text(site)
    row, col = site
    return "($row,$col)"
end

function generate(
        rows::Integer, cols::Integer;
        output_dir::AbstractString=joinpath(@__DIR__, "output"),
        io::IO=stdout)
    peps = toric_code_peps(rows, cols)
    circuit = sequential_circuit_graph(peps)
    block = yao_circuit(circuit)

    println(io, "lattice: $(rows)x$(cols)")
    println(io, "qubits: $(Yao.nqubits(block))")
    for (number, sites) in enumerate(circuit_layers(circuit))
        println(io, "layer $number: ", join(_site_text.(sites), ", "))
    end

    mkpath(output_dir)
    basename = "toric_circuit_$(rows)x$(cols)"
    outputs = String[]
    for format in (:svg, :pdf, :png)
        path = abspath(joinpath(output_dir, "$basename.$format"))
        Yao.vizcircuit(block; format=format, filename=path)
        push!(outputs, path)
        println(io, "wrote: $path")
    end
    return (; peps, circuit, block, outputs)
end

function main(
        args=ARGS;
        output_dir::AbstractString=joinpath(@__DIR__, "output"),
        io::IO=stdout,
        error_io::IO=stderr)
    dimensions = try
        _dimensions(args)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, "error: ", error.msg)
        println(error_io, USAGE)
        return 1
    end
    generate(dimensions...; output_dir=output_dir, io=io)
    return 0
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(GenerateToricCircuitExample.main())
end
```

- [ ] **Step 2: Run the focused test and verify GREEN**

Run:

```bash
julia --project=. -e 'using TopoNoise, ITensors, Test; include("test/example.jl")'
```

Expected: the runnable example testset passes, including native SVG/PDF/PNG output.

- [ ] **Step 3: Exercise the real CLI**

Run:

```bash
julia --project=. examples/generate_toric_circuit.jl 1 1
```

Expected: exit 0, report 6 qubits and layer 1 `(1,1)`, and write the three `toric_circuit_1x1` files below `examples/output/`.

- [ ] **Step 4: Verify invalid CLI arguments**

Run:

```bash
julia --project=. examples/generate_toric_circuit.jl 0 2
```

Expected: exit 1 and print `Usage:` plus the positive-integer constraint.

---

### Task 3: Integrate the example into the repository

**Files:**
- Modify: `.gitignore`
- Modify: `README.md`

**Interfaces:**
- Consumes: the CLI established in Task 2.
- Produces: ignored generated output and user-facing setup/run instructions.

- [ ] **Step 1: Ignore generated example figures**

Append this root-anchored entry to `.gitignore`:

```gitignore
/examples/output/
```

- [ ] **Step 2: Document setup, execution, and outputs**

Add a `Runnable example` subsection after the native Yao plotting example:

````markdown
### Runnable example

From the repository root, instantiate the project once and run the generator:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/generate_toric_circuit.jl 2 2
```

Omit `2 2` to use the default lattice size. The script prints the diagonal
schedule and writes SVG, PDF, and PNG files to `examples/output/`. It builds
and draws the circuit without allocating the exponentially large state vector.
````

- [ ] **Step 3: Re-run the focused example tests**

Run:

```bash
julia --project=. -e 'using TopoNoise, ITensors, Test; include("test/example.jl")'
```

Expected: all example assertions pass and no files are left in the repository output directory by the tests.

---

### Task 4: Final verification and review

**Files:**
- Verify: `examples/generate_toric_circuit.jl`
- Verify: `test/example.jl`
- Verify: `.gitignore`
- Verify: `README.md`

**Interfaces:**
- Consumes: all prior tasks.
- Produces: verified runnable example ready for handoff.

- [ ] **Step 1: Run the complete package suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: every existing and new testset passes with zero failures.

- [ ] **Step 2: Check repository whitespace and generated artifacts**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; `examples/output/` is absent from status because it is ignored.

- [ ] **Step 3: Request an independent code review**

Ask a read-only reviewer to check CLI behavior, Yao output generation, documentation accuracy, test quality, and preservation of existing APIs. Fix all critical and important findings, then rerun the affected focused tests.

- [ ] **Step 4: Run final verification after review fixes**

Run:

```bash
julia --project=. test/runtests.jl
git diff --check
```

Expected: all tests pass and the diff check emits no output.
