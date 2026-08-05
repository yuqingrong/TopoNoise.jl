# Toric-Code Local Isometry Design

## Goal

Represent the local doubled-edge toric-code tensor in the same directional
order as the reference diagram and normalize it so that it is an isometry.
The public finite-PEPS state must remain unchanged.

## Tensor convention

The stored order is

\[
(p_E,p_N,p_W,p_S,v_E,v_N,v_W,v_S)
=(i,j,k,l,\alpha,\beta,\gamma,\delta).
\]

Every copy-compatible even-parity entry is \(1/\sqrt{2}\); all other entries
are zero. `physicalinds` and all internal direction tuples use E–N–W–S.

## Isometry

Treat \((\alpha,\beta)\) as the four-dimensional input and
\((\gamma,\delta,i,j,k,l)\) as the 64-dimensional output. Permute the dense
tensor with `(7, 8, 1, 2, 3, 4, 5, 6)` and reshape it directly into

\[
T_{(\gamma\delta ijkl),(\alpha\beta)} \in \mathbb{C}^{64\times4}.
\]

Each input column has two nonzero entries of magnitude \(1/\sqrt{2}\), while
distinct columns have disjoint support. Therefore the required identity is

\[
T^\dagger T=I_4.
\]

The transposed \(4\times64\) co-isometry convention is not used in the test or
documentation. Requiring an identity on its 64-dimensional side would be
impossible because its rank is at most four. Only the \(64\times4\) isometric
orientation above is in scope.

Contracting the four physical legs with uniform covectors must directly give

\[
W=\frac{1}{\sqrt2}
\begin{pmatrix}
1&0&0&1\\
0&1&1&0\\
0&1&1&0\\
1&0&0&1
\end{pmatrix}.
\]

## Finite-state normalization

For \(V\) vertices and \(E\) internal edges, distribute the scale

\[
s=2^{1-E/(2V)}
\]

over every site. This compensates for the new \(1/\sqrt2\) local factor and
preserves all existing finite-state amplitudes and norms.

## Verification

- Exhaustively check all 256 local tensor entries.
- Reject integer local-tensor storage.
- Check the exact \(W\) matrix and the \(64\times4\) isometry \(T^\dagger T=I_4\).
- Remap every finite-patch direction-dependent test to E–N–W–S.
- Keep the complete small-patch and \(3\times3\) amplitudes and norms unchanged.
- Update README and `document/test.typ`, then compile the Typst source.
