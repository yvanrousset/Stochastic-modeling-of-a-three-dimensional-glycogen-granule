Changelog
===

Renaming parameters (to match the paper; or making them more descriptive)
---

- ``size_spec_gys`` -> `` l_gs_min ``
- ``size_spec_gbe_spacing`` -> ``l_gbe_spacing``
- ``size_spec_gbe_transferred`` -> ``l_gbe_transferred``
- ``size_spec_gbe_leftover`` -> ``l_gbe_leftover``

- ``r_GN`` -> ``radius_of_gn_core``
- ``b`` -> ``radius_of_the_glucose_sphere``


Constants
---
- ``GlycogenStructure.L = 0.24``
- parameter ``initial_number_of_chains`` -> ``GlycogenStructure.MAX_INITIAL_CHAINS = 2``
  - was always used as 2; cannot be larger than 2 in current model
- parameter ``b`` -> ``GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS = 1``
  - was always used as 1

Chains
---
- There is now a ``Chain`` class.


Chain status
---
- instead of an int, chains now have a list of ``Status``:
```
class Status(Enum):
    """Enum representing the status of a chain, showing which enzymatic reactions can take place."""
    GS_SUBSTRATE = 0
    GBE_SUBSTRATE = 1
    GP_SUBSTRATE = 2
    GDE_SUBSTRATE = 3
```
and can be used like this:
```
def find_chains_for_gde(self) -> list[int]:
        return [c.id for c in self.chains if Status.GDE_SUBSTRATE in c.substrate_status]
```