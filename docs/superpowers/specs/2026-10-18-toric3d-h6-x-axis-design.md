# Physical x-axis H+6 growth comparison

The user requests sequential growth along the physical x axis instead of the existing diagonal growth. This is not an X-only noise scan and not an aspect-ratio scan. Preserve the completed diagonal study unchanged.

## Circuit and geometry

Use a periodic square yz toric-code seed with two logical qubits, plus T layers of cubic-lattice growth. There are N=Ly*Lz vertices per plane. Seed Q0 contains 2N yz edges. Each Qr, r>=1, contains N longitudinal x links from plane r-1 to r and 2N yz edges on plane r. Total qubits (3T+2)N. Boundary stars have five edges; interior stars have six. All square plaquettes and the two transverse logical pairs define the code.

Each transfer first copies every old yz edge into its corresponding fresh yz edge (2N CNOTs). Then, for every vertex, apply H to its fresh x link and four CNOTs from that link to the incident fresh yz edges. Thus each cell uses H+6 CNOTs and thirteen independent noisy gate outputs. All copy gates precede all vertex gates. Fix all edge and vertex iteration orders for reproducibility.

The X channel starts from ideal logical |00>. The Z channel starts from logical |++>, with independent Z faults at the 2N seed-H outputs and common seed-star MWPM repair shared by all three arms. Retain every trial and any undetectable logical seed error. Measurements, corrections, idles and other seed operations are ideal, as in the diagonal study.

## Feedback

- Final-only: no growth feedback, one final complete recovery.
- X one-step: before each transfer, measure old-plane square Z plaquettes, MWPM pair defects, take the endpoint edge(s) of each deterministic shortest path, and apply the propagated X route only to fresh yz targets after the transfer. This is the square-slice analog of the diagonal endpoint-routing policy.
- Z one-step: after transfers 1:T-1, measure retiring-plane physical stars and move each violated star one physical x step using the corresponding new longitudinal Z edge. This orientation-adapted rule must be named explicitly; the diagonal matching rule cannot simply be reused.
- Two-slice: after every transfer including the last, repair only Q(r-1) union Qr using all incident physical prefix checks. X uses BP+OSD; Z uses unit-weight MWPM. Audit every completed-prefix check, but do not supply those audit results or latent faults to local decoders.

The final decoder is shared across arms within each sector: full physical plaquette BP+OSD for X, full physical star MWPM for Z. Use the same decoder settings as the archived diagonal study.

## Comparison and limits

Use T=Ly=Lz=L, L=3,5,7,9,11 and 1,000 shared fault records per point for the three axis arms. Begin with verification and excluded timing pilots; then use the completed diagonal study's physical-p grids if runtime is reasonable. New records use independent recorded seeds. Across directions this is an independent-sample comparison; only the three arms within each direction share records.

Equal L gives equal growth-cell count and transfer gate-output count, but not identical seed geometry or logical distance: axis dZ=min(Ly,Lz), diagonal dZ=2min(Ly,Lz). Both have dX=(T+1)min(Ly,Lz). Report these differences and seed exposure (2N versus 3N); do not call this a direction-only control at fixed code distance. The two-slice X observable remains a preparation outcome fixed after the first corrected window, not a memory threshold in T.

Keep original diagonal estimates diagnostic. An axis scout does not establish a threshold merely by fitting. Use PDF figures only, with the archived paper-no-grid theme, one comparison PDF and compact data/figures/archive organization.

## Required evidence

Prove CSS commutation, k=2, logical commutation, every-prefix signed stabilizer/logical pullback for |00> and |++>, exact H/CNOT/output counts, and no targets in retired planes. Test noiseless outcomes, hand-injected single faults, immutability and shared records, seed-logical retention, one-step retirement, two-slice complete-prefix repair, final-score logical consistency, and exact checkpoint resume. Freeze source, protocol, selection and runtime hashes before production.
