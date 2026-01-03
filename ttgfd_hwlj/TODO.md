- Move `alj` and `rlj` into `calculate_lj_potential_table` subroutine (global -> local variables)
- Move `cos_pphi` and `npphi` into `calculate_lj_potential_table` subroutine. The following two "write" statements to stdout should also be moved into the subroutine.
- Make subroutine for constructing `input` by opening/closing `iep` and `iep` file handles. These handles can be local instead of global.
Why are `cdens` and `ctvec` allocated on the stack?
- `initialize_computed_params` could take `bulk` as input instead of three floats.
