# ToQC / mpo_unitary (MATLAB)

The ToQC/mpo_unitary folder contains a **Matrix Product Operator (MPO)** backend for **time propagation and optimal control of unitary evolution operators** (as MPOs), instead of propagating quantum states (MPS).

It provides:
- **TEBD-style propagation of an MPO** under nearest-neighbor dynamics and local controls
- **MPO compression / truncation** to control bond growth
- **Infidelity + gradient** evaluation for optimal control, including **robust / ensemble** variants

> **MPO representation:** an MPO is a `1×N` cell array where each site tensor is typically `Dl × d × d × Dr`.

---

## Requirements

- **MATLAB R2022a+** (uses `tensorprod` with name-value arguments).
- Parallel Computing Toolbox is recommended for robust routines that use `parfor`.
- Add the repository to your MATLAB path:
  ```matlab
  addpath(genpath(pwd))
  ```

---

## Folder contents

### TEBD / options
- `tebd_default_options.m` — default TEBD + truncation options (`struct`)

### Core MPO operations
- `gate_1q.m` — apply a 1-site gate to one MPO tensor
- `gate_2q.m` — apply a 2-site nearest-neighbor gate + SVD truncation
- `mpo_overlap.m` — overlap-like trace contraction (see below)
- `mpo_normalize.m` — normalize an MPO (global rescaling distributed across sites)
- `mpo_compress.m` — sweep-based recompression/truncation
- `conjtp.m` — conjugate transpose of an MPO (swaps input/output legs)

### Time evolution
- `mpo_evol.m` — TEBD evolution of an MPO with piecewise-constant controls

### Objectives / gradients (optimal control)
- `infid.m` — infidelity + gradient (optionally includes ∂/∂T via finite difference)
- `infid_nograd.m` — infidelity only (+ tracked max tensor/bond size)

### Robust / ensemble (PRL) variants
- `infid_robust_prl.m` — robust infidelity + gradient over samples (parfor)
- `infid_nograd_robust_prl.m` — robust infidelity only (+ bond stats)

---

## Data structures and conventions

### MPO (operator)
Each MPO tensor uses the index order:
```
1: left bond     (Dl)
2: bottom / input (d)
3: top / output   (d)
4: right bond    (Dr)
```

In MATLAB:
```matlab
mpo{j}    % site j tensor, typically size [Dl, d, d, Dr]
```

> Note: edge tensors may appear as 3-D arrays because MATLAB drops trailing singleton dimensions. The implementation accounts for common edge cases.

### Conjugate transpose (`conjtp`)
`conjtp(mpo)` performs a complex conjugate and swaps the bottom/top indices:
```matlab
mpo{j} -> permute(conj(mpo{j}), [1, 3, 2, 4])
```

### Two-site (nearest-neighbor) terms
There are two common Hamiltonian parameterizations in this branch:

1) **Single fixed drift term per bond** (used by `infid` / `infid_nograd`):
```matlab
H0{j}     % acts on (j, j+1), stored as a (d^2 × d^2) matrix
```

2) **Linear combination of interaction “basis” terms per bond** (used by `mpo_evol` and the robust PRL functions):
```matlab
H0{j,l}   % l = 1..m0 interaction basis terms on bond j
c0(k,l)   % piecewise-constant weights per time bin k
```
so the bond Hamiltonian in bin `k` is:
```
Hbond(j,k) = Σ_l c0(k,l) * H0{j,l}
```

### Controls (`Hc`)
Controls are passed as a struct array `Hc(1:nc)` with fields:
- `Hc(jc).sys` : vector of site indices where control `jc` acts
- `Hc(jc).op`  : cell array of local operators (`d×d`) aligned with `.sys`

Example (control `jc` acts on sites 2 and 5):
```matlab
Hc(jc).sys = [2 5];
Hc(jc).op  = {X, Z};   % X on site 2, Z on site 5
```

### Piecewise-constant controls in time
Control amplitudes are typically stored as a flattened vector and reshaped internally to `[nbin, nc]`:
```matlab
nbin = length(c) / nc;
c    = reshape(c, [nbin, nc]);
```

---

## TEBD options (`tebd_default_options`)

Create defaults:
```matlab
tebd_options = tebd_default_options();
```

Fields (defaults shown):
- `sv_min` = `1e-10` : truncation threshold based on discarded singular-value weight
- `bond_dim` = `20`  : max bond dimension during two-site gates
- `bond_comp` = `20` : max bond dimension during compression sweeps
- `num_sweep` = `2`  : number of compression sweeps
- `num_midstep` = `1`: compress/normalize every `midstep` bins
- `num_refined_step` = `1`: substeps per bin (`nt`)
- `is_compressed` = `1`: enable periodic recompression
- `is_second_order` = `0`: enable second-order Trotter (adds half-steps)

---

## Function reference

### `gate_1q`
```matlab
A = gate_1q(A0, tno)
```
Applies a one-site gate `tno` (`d×d`) to a single MPO tensor `A0` (index order `[left, bottom, top, right]`).

---

### `gate_2q`
```matlab
[A, B] = gate_2q(A0, B0, tno, sv_min, D)
```
Applies a nearest-neighbor two-site gate `tno` and truncates the resulting bond via SVD.
- `tno` is a rank-4 tensor of size `d×d×d×d` (for qubits: `2×2×2×2`)
- `sv_min` controls truncation by discarded weight
- `D` caps the bond dimension

---

### `mpo_compress`
```matlab
mpo1 = mpo_compress(mpo, sv_min, D, nsweep)
```
Sweep-based MPO recompression/truncation (left-to-right and right-to-left sweeps).

---

### `mpo_overlap`
```matlab
O = mpo_overlap(mpo, mpo0)
```
Contracts two MPOs into a scalar using a trace-like contraction.
Internally it forms `conjtp(mpo)` and contracts bottom/top legs such that `O` corresponds to a **Hilbert–Schmidt inner product** (trace overlap), i.e. proportional to:
```
O ∝ tr( mpo† * mpo0 )
```
The optimization routines use `abs(O)^2` as the fidelity-like quantity (after normalization).

---

### `mpo_normalize`
```matlab
[mpo2, new_norm, old_norm] = mpo_normalize(mpo1)
```
Normalizes an MPO by distributing a global rescaling over all sites.

---

## Time evolution

### `mpo_evol`
```matlab
[mpo, Dlist] = mpo_evol(H0, Hc, c0, c, T, mpo0, tebd_options)
```

Evolves an initial MPO `mpo0` forward in time using TEBD (odd/even two-site gates + one-site control gates).

Inputs:
- `H0` : `(N−1)×m0` cell of two-site interaction basis terms, each `(d^2×d^2)`
- `Hc` : control struct array (see above)
- `c0` : `nbin×m0` drift weights per bin (scales `H0{j,l}`)
- `c`  : flattened control amplitudes length `nbin*nc` (reshaped to `[nbin,nc]`)
- `T`  : total evolution time
- `mpo0` : initial MPO (cell array)
- `tebd_options` : options struct

Outputs:
- `mpo` : final evolved MPO
- `Dlist` : tracked max tensor size at normalization/compression checkpoints

Notes:
- `dt = T / (nbin * num_refined_step)`
- If `is_second_order == 1`, the routine performs initial/final half-steps for second-order Trotterization.

---

## Infidelity + gradients (optimal control)

### `infid`
```matlab
[iF, iG] = infid(H0, Hc, x, mpo0, mpotg, varT, tebd_options)
```

Computes infidelity and gradient for unitary-MPO control.

Parameter vector:
- `x = [c(:); T]` where `c(:)` has length `nbin*nc`

Infidelity used in this code:
- `ovl = mpo_overlap(mpotg, mpo_final)`
- `iF = 1 - |ovl|^2`

The gradient `iG` is returned as a column vector of length `nbin*nc`, and (optionally) an extra element for time:
- If `varT == 1`, the routine estimates `∂iF/∂T` by a finite difference and appends it to `iG`.
- If `varT == 0`, the appended time-gradient term is `0`.

Inputs:
- `H0` : `1×(N−1)` cell of fixed two-site drift terms (`d^2×d^2`)
- `Hc` : control struct array
- `mpo0` : initial MPO (often identity, or a starting unitary)
- `mpotg` : target MPO (desired unitary)

---

### `infid_nograd`
```matlab
[iF, maxD] = infid_nograd(H0, Hc, x, mpo0, mpotg, tebd_options)
```
Infidelity only (no gradient) and tracks the maximum observed tensor size (`maxD`).

---

## Robust / ensemble variants (PRL)

These functions evaluate performance over an ensemble of uncertain couplings. They use `parfor` by default.

### Interaction model
For sample `s` and bond `j`, the two-site Hamiltonian is:
```
Hbond(s,j) = H00{j} + Σ_l J0(s,j,l) * Hu0{j,l}
```

Inputs:
- `J0`  : `(nsampl × (N−1) × num_int)` coupling coefficients
- `H00` : `1×(N−1)` baseline two-site terms
- `Hu0` : cell array of size `((N−1) × num_int)` interaction basis terms
- `Hc0` : controls struct array
- `x`   : parameter vector `x = [c(:); T]`
- `mpo00`, `mpotg0` : initial and target MPOs
- `ismean` : toggles returning averages vs per-sample values

### `infid_robust_prl`
```matlab
[iFavg, iGavg] = infid_robust_prl(J0, H00, Hu0, Hc0, x, mpo00, mpotg0, varT, ismean, tebd_options)
```
- If `ismean == 1`: returns mean infidelity and mean gradient across samples.
- If `ismean == 0`: returns per-sample infidelities and per-sample gradients (stacked by sample).

### `infid_nograd_robust_prl`
```matlab
[iFlist, Dlist] = infid_nograd_robust_prl(J0, H00, Hu0, Hc0, x, mpo00, mpotg0, ismean, tebd_options)
```
Infidelity only (+ tensor/bond stats) across samples.

---

## Minimal usage sketch

```matlab
tebd_options = tebd_default_options();

% Build: H0 (two-site cells), Hc (controls), initial mpo0, target mpotg
% Choose nbin, nc, controls c(:) of length nbin*nc, and total time T

x = [c(:); T];

[iF, iG] = infid(H0, Hc, x, mpo0, mpotg, 1, tebd_options);
```

---

## Tips / troubleshooting

- **Older MATLAB releases:** if `tensorprod` is missing, upgrade to R2022a+ or replace `tensorprod` with an equivalent contraction routine.
- **Robust routines fail to run in parallel:** replace `parfor` with `for` (slower) or install Parallel Computing Toolbox.
- **Keeping bond dimension under control:** tune `bond_dim`, `bond_comp`, and `sv_min` and enable `is_compressed`.

---

## License / citation

See the repository root for license information.  
If you use TiQC‑ToQC in academic work, please cite the repository and associated publications (if applicable).
