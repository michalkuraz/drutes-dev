# Project context

Read `docs/PROJECT_CONTEXT.md` for the architecture, GUI workflow, and recent
NetCDF river-projection work recorded at the user's request. This is a dated
snapshot; check the current source and Git status before relying on it.

Runtime caution: `src/core/main.f90` clears the working directory's `out/`
contents at startup. Do not run the model in an existing results directory
merely to inspect the project.
