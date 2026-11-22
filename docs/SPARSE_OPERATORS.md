# Shared sparse operators

`src/mathtools/sparse_operators.f90` supplies `csr_matrix` and
`diagonal_matrix`, both extending the existing `mtx::matrix` interface.
It contains no iterative solver and has no dependency on any physical model.
The hydraulic projection calls the existing `solvers::PCG` with these objects.
The transport GMRES and the default multiplication of existing `smtx` objects
are not replaced.

## Assembly and repeated multiplication

Use the existing `smtx` for assembly, then explicitly take a CSR snapshot:

```fortran
use sparsematrix, only: smtx
use sparse_operators, only: csr_matrix, diagonal_matrix
use solvers, only: PCG

type(smtx) :: assembled
type(csr_matrix) :: operator
type(diagonal_matrix) :: inverse

! Assemble using the existing matrix methods.
call operator%from_matrix(assembled)
call inverse%set_values(1.0_rkind / diagonal)
call PCG(operator, inverse, rhs, solution, errcode1=status)
```

`set_values` stores the supplied diagonal; it does not invert it automatically.
The example supplies the inverse Jacobi diagonal explicitly. CSR multiplication
walks contiguous row entries without extracting/copying each row. Diagonal
multiplication is an elementwise vector product. Both preserve the established
allocatable-result `matrix%mul` interface and optional operation counters.
CSR also supports rectangular matrices and transpose multiplication.

## Updates and lifetime

- `from_matrix` copies the source, retaining its row-entry order. Later source
  changes do not affect the snapshot. Call it again after changing the source,
  including a change in sparsity or dimensions. Two passes allocate exactly the
  required storage, including when the source is dense but contains zeros.
- `csr%set` updates an existing structural entry. Adding a new nonzero outside
  the stored pattern is rejected: assemble it in `smtx` and reload explicitly.
- `clear` zeros values but retains the CSR pattern for coefficient updates.
  `nonzero` and `getrow` report numerical nonzeros, not structural zero slots.
- `init` and `resize` both reset dimensions and storage. Zero-sized matrices
  and empty rows are supported. Products, indices and masks are checked.
- Storage belongs to each object; there are no shared saved work arrays or
  implicit source pointers/caches.

The existing base constructor `mtx::initmatrix` is exported so derived
constructors can set dimensions without recursively dispatching to their own
`init`/`resize` overrides. No binding was added to the base type's dispatch table.

## Hydraulic integration

The projection keeps its global element IDs, component-compatibility checks,
closed-component gauges, fixed fluxes and nonnegative-outlet active set.
Only the algebraic unknowns are compacted: inactive and pinned zero potentials
are omitted, and solved potentials are scattered back before recovering fluxes.
Thus 35,137 mesh elements do not require 35,137 algebraic rows when only 5,874
elements are active. A fully pinned compatible system needs no PCG invocation.
The same configured absolute residual target is supplied to shared PCG.

This optimization does not alter hydrology, storage, transport discretization,
dispersion, stabilization, time stepping or configuration formats.

## Tests and timing

- `tests/test_sparse_operators.py`: checked, floating-point-trapping builds;
  rectangular/transpose products, empty rows, masks, updates, reloads, resize,
  operation counts, dimension/index/pattern guards, and PCG equivalence.
- `tests/test_pcg.py`: independent shared-solver tests.
- `tests/run_adenc_hydroflow_checks.py`: production hydraulic, bank, conservative,
  lateral-source and compact-global-ID tests.
- `tests/fortran/benchmark_adenc_projection.f90`: replay exported production
  projection inputs against separately linked original/current object builds.
  It does not contain another PCG implementation.
- `tests/run_adenc_solver_regression.py`: isolated paired transport runs with
  immutable private executables, NetCDF symlinks, observation/field/budget checks
  and preserved outputs. Optional reference-result reuse never reruns its model.

Actual timings and numerical differences belong to the run report, not a fixed
performance guarantee. Compare equivalent input revisions: historical Case I
used different width and routing conventions from the current June cases.
# Regression budget tolerances

`tests/run_adenc_solver_regression.py` checks each run's reported and
independently reconstructed mass budget. Cancellation residuals use
`1e-9 + 1e-12 * sum(abs(balance terms))` in audit solute-amount units.
Step scales include the current and previous inventory and step exchanges;
cumulative scales include initial/current inventory and cumulative exchanges.
Both reported and reconstructed residuals, and their agreement, must pass.
Nonfinite values are rejected; the relative audit error must remain below
`1e-10`. Concentration, inventory, exchange and time comparisons are unchanged.
Use `--check-existing` to recheck saved paired runs without executing the model.

