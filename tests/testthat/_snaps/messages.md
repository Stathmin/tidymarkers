# abort messages are stable

    Code
      with_unafloddat("/definitely/not/here", 1)
    Condition
      Error in `with_unafloddat()`:
      ! Directory not found: /definitely/not/here

---

    Code
      kasp_mix_recipe(n_reactions = 1L, mix_uL = 1, stock_uM = 1)
    Condition
      Error in `kasp_mix_recipe()`:
      ! Requested working mix (12/12/30 uM) exceeds the 1 uM stock.
      i Raise `stock_uM` or lower the targets.

