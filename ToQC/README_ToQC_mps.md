# ToQC / mps (MATLAB)

This folder contains the **Matrix Product State (MPS)** backend of **TiQC‑ToQC** for simulating and optimizing 1D quantum dynamics.  
It supports **TEBD-style time evolution**, MPS **compression/truncation**, and **infidelity + gradient** evaluation for (optimal) control, including **robust / ensemble** variants.

> **MPS representation:** an MPS is a `1×N` cell array where each site tensor is typically `Dl × d × Dr`.

---

## Requirements

- **MATLAB R2022a+** (the code uses `tensorprod` with name-value arguments).
- Parallel Computing Toolbox is recommended for the robust routines that use `parfor`.
- Add the repository to your MATLAB path:
  ```matlab
  addpath(genpath(pwd))
  ```

---

## Folder contents

### TEBD / options
- `tebd_default_options.m` — default TEBD + truncation options (`struct`)

### Core MPS operations
- `mps_gate_1q.m` — apply a 1-site gate to one MPS tensor
- `mps_gate_2q.m` — apply a 2-site nearest-neighbor gate + SVD truncation
- `mps_overlap.m` — overlap ⟨mps | mps0⟩
- `mps_normalize.m` — normalize an MPS (global rescaling distributed across sites)
- `mps_compress.m` — sweep-based recompression/truncation

### Time evolution
- `mps_evol.m` — TEBD evolution with piecewise-constant controls (and optional 2nd-order Trotter)

### Objectives / gradients (optimal control)
- `mps_infid.m` — infidelity + gradient (optionally includes ∂/∂T via finite difference)
- `mps_infid_nograd.m` — infidelity only (+ tracked max tensor/bond size)

### Robust / ensemble (PRL) variants
- `mps_infid_robust_prl.m` — robust infidelity + gradient over samples (parfor)
- `mps_infid_nograd_robust_prl.m` — robust infidelity only (+ bond stats)
- `mps_infid_T_robust_prl.m` — robust infidelity + gradient, parameter vector includes total time `T`
- `mps_infid_T_nograd_robust_prl.m` — robust infidelity only (with `T` in parameter vector)

---

## Data structures and conventions

### MPS (state)
```matlab
mps{j}    % site j tensor, typically size [Dl, d, Dr]
```
MATLAB may drop trailing singleton dimensions at the edges; the routines handle common edge cases.

### Two-site (nearest-neighbor) terms
Most evolution/infidelity routines use a cell array of two-site operators:
```matlab
H0{j}     % acts on sites (j, j+1), stored as a (d^2 × d^2) matrix
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
Controls are stored as a flattened vector and reshaped internally to `[nbin, nc]`:
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

### `mps_gate_1q`
```matlab
A = mps_gate_1q(A0, tno)
```
Applies a one-site gate `tno` (`d×d`) to the physical index of tensor `A0`.

---

### `mps_gate_2q`
```matlab
[A, B] = mps_gate_2q(A0, B0, tno, sv_min, D)
```
Applies a nearest-neighbor two-site gate `tno` and truncates the resulting bond via SVD.
- `tno` is a rank-4 tensor of size `d×d×d×d` (for qubits: `2×2×2×2`)
- `sv_min` controls truncation by discarded weight
- `D` caps the bond dimension

---

### `mps_compress`
```matlab
mps1 = mps_compress(mps, sv_min, Dc, nsweep)
```
Sweep-based MPS recompression/truncation (left-to-right and right-to-left sweeps).

---

### `mps_overlap`
```matlab
O = mps_overlap(mps, mps0)
```
Returns the scalar overlap ⟨mps | mps0⟩.

---

### `mps_normalize`
```matlab
[mps2, new_norm, old_norm] = mps_normalize(mps1)
```
Normalizes an MPS by distributing a global rescaling over all sites.

---

## Time evolution

### `mps_evol`
```matlab
[mps, Dlist] = mps_evol(H0, Hc, c0, c, T, mps0, tebd_options)
```

Evolves `mps0` forward in time using TEBD (odd/even two-site gates + one-site control gates).

Inputs:
- `H0` : `1×(N−1)` cell, each entry is `(d^2 × d^2)` for bond `(j,j+1)`
- `Hc` : control struct array (see “Controls” above)
- `c0` : length-`nbin` drift scaling per time bin (multiplies `H0{j}`)
- `c`  : flattened control amplitudes of length `nbin*nc` (reshaped to `[nbin,nc]`)
- `T`  : total evolution time
- `mps0` : initial MPS (cell array)
- `tebd_options` : options struct

Outputs:
- `mps` : final evolved MPS
- `Dlist` : tracked max tensor size at normalization/compression checkpoints

Notes:
- `dt = T / (nbin * num_refined_step)`
- If `is_second_order == 1`, the routine performs initial/final half-steps for second-order Trotterization.

---

## Infidelity + gradients (optimal control)

### `mps_infid`
```matlab
[iF, iG] = mps_infid(H0, Hc, x, mps0, mpstg, varT, tebd_options)
```

Computes infidelity and gradient for controls (and optionally for total time):
- Parameter vector: `x = [c(:); T]` where `c(:)` has length `nbin*nc`
- Infidelity: `iF = 1 - |⟨mpstg | U(x) | mps0⟩|^2`
- Gradient: `iG` is `nbin×nc` (the code also supports an optional `∂/∂T` term)

`varT`:
- `varT = 1` includes a finite-difference estimate of `∂iF/∂T`
- `varT = 0` sets the `T` contribution to zero

---

### `mps_infid_nograd`
```matlab
[iF, maxD] = mps_infid_nograd(H0, Hc, x, mps0, mpstg, tebd_options)
```
Infidelity only (no gradient) and tracks the maximum observed tensor size (`maxD`).

---

## Robust / ensemble variants (PRL)

These functions evaluate performance over an ensemble of uncertain couplings. They use `parfor` by default.

### Interaction model
The per-bond Hamiltonian for sample `s` on bond `j` is constructed as:
```
Hbond(s,j) = H00{j} + Σ_l J0(s,j,l) * Hu0{j,l}
```

Inputs:
- `J0`  : `(nsampl × (N−1) × num_int)` coupling coefficients
- `H00` : `1×(N−1)` baseline two-site terms
- `Hu0` : cell array of size `((N−1) × num_int)` interaction basis terms
- `Hc0` : controls struct array (same format as `Hc`)
- `c0`  : flattened control amplitudes (`nbin*nc`), reshaped to `[nbin,nc]`
- `T`   : total time
- `mps00`, `mpstg0` : initial and target MPS
- `ismean` : toggles returning averages vs per-sample values

### `mps_infid_robust_prl`
```matlab
[iFavg, iGavg] = mps_infid_robust_prl(J0, H00, Hu0, Hc0, c0, T, mps00, mpstg0, ismean, tebd_options)
```
- If `ismean == 1`: returns mean infidelity and mean gradient across samples.
- If `ismean == 0`: returns per-sample infidelities and per-sample gradients.

### `mps_infid_nograd_robust_prl`
```matlab
[iFlist, Dlist] = mps_infid_nograd_robust_prl(J0, H00, Hu0, Hc0, c0, T, mps00, mpstg0, ismean, tebd_options)
```
Infidelity only (+ tensor/bond stats) across samples.

### Time-in-parameter-vector wrappers
- `mps_infid_T_robust_prl(J, H0, Hu, Hc, x, mps0, mpstg, varT, ismean, tebd_options)`
- `mps_infid_T_nograd_robust_prl(J0, H00, Hu0, Hc0, x, mps00, mpstg0, ismean, tebd_options)`

These variants expect `x = [c(:); T]` similarly to the non-robust `mps_infid`.

---

## Minimal usage sketch

```matlab
tebd_options = tebd_default_options();

% Build: H0 (two-site cells), Hc (controls), initial mps0, target mpstg
% Choose nbin, nc, controls c(:) of length nbin*nc, and total time T

x = [c(:); T];

[iF, iG] = mps_infid(H0, Hc, x, mps0, mpstg, 1, tebd_options);
```

---

## Tips / troubleshooting

- **Older MATLAB releases:** if `tensorprod` is missing, upgrade to R2022a+ or replace `tensorprod` with an equivalent contraction routine.
- **Robust routines fail to run in parallel:** replace `parfor` with `for` (slower) or install Parallel Computing Toolbox.
- **Keeping bond dimension under control:** tune `bond_dim`, `bond_comp`, and `sv_min` and enable `is_compressed`.

---

## License / citation
 
If you use TiQC‑ToQC in academic work, please cite the repository and associated publications (if applicable).
