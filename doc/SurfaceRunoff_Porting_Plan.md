# Plan: Port `SurfaceRunoff` to URBANxx

## Goal

Port ELM's `SurfaceRunoff` subroutine to URBANxx for all five urban surfaces
(roof, sunlit wall, shaded wall, impervious road, pervious road).  
The output `qflx_surf` is computed by URBANxx and compared against ELM's value
via a `_check` routine.

---

## Assumptions for this port

| Assumption | Value | Impact |
|---|---|---|
| `use_vichydro` | `.false.` | VIC branch not needed |
| `origflag` | `0` | `fcov = fsat`; no ice-fraction path |
| `lun_pp%ispolygon(l)` | `.false.` | Polygonal-ground runoff=0 branch skipped |
| `qflx_floodc(c)` | `0.0` | Flood term drops out everywhere |
| `use_IM2_hillslope_hydrology` | `0` | Uphill lateral-flow block skipped |
| `snl(c)` | `0` (no snow) | Snow-layer branch (`snl < 0`) never taken for roof/imperv road |
| `qflx_top_soil` for urban surfaces | `= ForcRain` (atmospheric rain) | Already in `urban->atmosphereData.ForcRain`; no new setter needed |

---

## Simplified ELM logic after applying assumptions

### Pervious road  (`filter_hydrologyc`: c = 6, 11, 16)

```
fff  = fover(g)          ! decay factor, indexed by gridcell g
fsat = wtfact(c) * exp(-0.5 * fff * zwt(c))

! Perched water table adjustment
if frost_table(c) > zwt_perched(c):
    fsat = wtfact(c) * exp(-0.5 * fff * zwt_perched(c))

qflx_surf(c) = fsat * qflx_top_soil(c)    ! = fsat * ForcRain
```

### Roof / impervious road  (`filter_urbanc`, `snl == 0` always)

```
! snl == 0 branch only (no snow layers)
xs = max(0, h2osoi_liq(c,1)/dt + qflx_top_soil(c) - qflx_evap_grnd(c) - pondmx_urban/dt)
if xs > 0:
    h2osoi_liq(c,1)  = pondmx_urban               ! TopH2OSoiLiq updated
else:
    h2osoi_liq(c,1) += (qflx_top_soil(c) - qflx_evap_grnd(c)) * dt   ! TopH2OSoiLiq updated
qflx_surf(c) = xs
```

`qflx_evap_grnd` is **already available** inside URBANxx as `QflxEvapGrnd` on
roof and imperviousRoad surfaces — it is computed by `UrbanComputeSoilFluxes`
(called in `urbanxx_soilFluxes`) before `SurfaceRunoff` runs.

### Sunlit wall / shaded wall

```
qflx_surf(c) = 0
```

---

## Data flow summary

| ELM variable | Fortran source | URBANxx field | Direction | Already available? |
|---|---|---|---|---|
| `qflx_top_soil(c)` | `= ForcRain` for urban | `atmosphereData.ForcRain` | — | **Yes** — set by `UrbanSetAtmRain` each step |
| `h2osoi_liq(c,1)` input (roof/imperv) | `col_ws%h2osoi_liq` | `roof.TopH2OSoiLiq` / `imperviousRoad.TopH2OSoiLiq` | ELM → URBANxx | Yes — `UrbanSetTopH2OSoiLiqRoof/ImperviousRoad` |
| `qflx_evap_grnd(c)` (roof/imperv) | `col_wf%qflx_evap_grnd` | `roof.QflxEvapGrnd` / `imperviousRoad.QflxEvapGrnd` | — | **Yes** — computed by `UrbanComputeSoilFluxes` |
| `wtfact(c)` | `soilstate_vars%wtfact_col` | `perviousRoad.Wtfact` *(new)* | ELM → URBANxx | No |
| `fover(g)` | `soilhydrology_vars%fover` | `perviousRoad.Fover` *(new)* | ELM → URBANxx | No |
| `frost_table(c)` | `soilhydrology_vars%frost_table_col` | `perviousRoad.FrostTable` *(new)* | ELM → URBANxx | No |
| `zwt_perched(c)` | `soilhydrology_vars%zwt_perched_col` | `perviousRoad.ZwtPerched` *(new)* | ELM → URBANxx | No |
| `zwt(c)` | `soilhydrology_vars%zwt_col` | `perviousRoad.Zwt` | ELM → URBANxx | Yes — `UrbanSetWaterTableDepth` |
| **`qflx_surf(c)`** *(output)* | `col_wf%qflx_surf` | `roof/imperviousRoad/perviousRoad/wall.QflxSurf` *(new)* | URBANxx → ELM | No |
| **`h2osoi_liq(c,1)`** *(output, roof/imperv)* | `col_ws%h2osoi_liq` | `roof.TopH2OSoiLiq` / `imperviousRoad.TopH2OSoiLiq` | URBANxx → ELM | Field exists — reused as output |

> **`snl` is not needed**: with `snl == 0` assumed, `snl` is never tested and
> no new `Snl` view is required.

---

## Implementation steps

### Step 1 — Add new Kokkos View fields to `UrbanSurfaceTypeImpl.h`

**`RoofDataType`** — 1 new field:
- `QflxSurf` — `Array1DR8, numLandunits` — surface runoff output (mm/s)

**`ImperviousRoadDataType`** — 1 new field:
- `QflxSurf` — `Array1DR8, numLandunits` — surface runoff output (mm/s)

**`PerviousRoadDataType`** — 5 new fields:
- `Wtfact` — `Array1DR8, numLandunits` — max saturated fraction
- `Fover` — `Array1DR8, numLandunits` — decay factor (m)
- `FrostTable` — `Array1DR8, numLandunits` — frost table depth (m)
- `ZwtPerched` — `Array1DR8, numLandunits` — perched water table depth (m)
- `QflxSurf` — `Array1DR8, numLandunits` — surface runoff output (mm/s)

**`WallDataType`** — 1 new field:
- `QflxSurf` — `Array1DR8, numLandunits` — always zero; allocated for uniform getter interface

### Step 2 — New C API declarations in `include/Urban.h` (and `for_e3sm/include/Urban.h`)

**New setters (4)**:

```c
URBAN_EXTERN void UrbanSetWtfactPerviousRoad(UrbanType urban, const double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetFoverPerviousRoad(UrbanType urban, const double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetFrostTablePerviousRoad(UrbanType urban, const double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanSetZwtPerchedPerviousRoad(UrbanType urban, const double *values, int length, UrbanErrorCode *status);
```

**New compute function (1)**:

```c
URBAN_EXTERN void UrbanComputeSurfaceRunoff(UrbanType urban, double dtime, UrbanErrorCode *status);
```

**New getters (7)**:

```c
URBAN_EXTERN void UrbanGetQflxSurfRoof(UrbanType urban, double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetQflxSurfImperviousRoad(UrbanType urban, double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetQflxSurfPerviousRoad(UrbanType urban, double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetQflxSurfSunlitWall(UrbanType urban, double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetQflxSurfShadedWall(UrbanType urban, double *values, int length, UrbanErrorCode *status);
// Updated TopH2OSoiLiq after ponding (roof and impervious road)
URBAN_EXTERN void UrbanGetTopH2OSoiLiqRoof(UrbanType urban, double *values, int length, UrbanErrorCode *status);
URBAN_EXTERN void UrbanGetTopH2OSoiLiqImperviousRoad(UrbanType urban, double *values, int length, UrbanErrorCode *status);
```

### Step 3 — New C++ source files in `src/`

#### `UrbanSurfaceRunoffSetters.cpp`

Four 1D-double setters using `SetView1D` for the pervious road fields:
`Wtfact`, `Fover`, `FrostTable`, `ZwtPerched`.

#### `UrbanSurfaceRunoff.cpp`

`UrbanComputeSurfaceRunoff(UrbanType urban, Real dtime)` — single
`Kokkos::parallel_for` over `nlandunits`:

```cpp
// Roof and impervious road (snl == 0 assumed)
// pondmx_urban = PONDMX_URBAN constant from UrbanConstants.h
for each l:
    // Roof
    auto top_liq_roof   = roof.TopH2OSoiLiq(l);
    auto qflx_top_roof  = atm.ForcRain(l);        // qflx_top_soil = ForcRain
    auto evap_grnd_roof = roof.QflxEvapGrnd(l);
    Real xs_roof = max(0, top_liq_roof/dtime + qflx_top_roof - evap_grnd_roof - PONDMX_URBAN/dtime);
    if (xs_roof > 0):
        roof.TopH2OSoiLiq(l) = PONDMX_URBAN;
    else:
        roof.TopH2OSoiLiq(l) = max(0, top_liq_roof + (qflx_top_roof - evap_grnd_roof)*dtime);
    roof.QflxSurf(l) = xs_roof;

    // Impervious road  (same pattern)
    ...

    // Pervious road
    Real fff  = perv.Fover(l);
    Real fsat = perv.Wtfact(l) * exp(-0.5 * fff * perv.Zwt(l));
    if (perv.FrostTable(l) > perv.ZwtPerched(l)):
        fsat = perv.Wtfact(l) * exp(-0.5 * fff * perv.ZwtPerched(l));
    perv.QflxSurf(l) = fsat * atm.ForcRain(l);

    // Walls always zero
    sunlitWall.QflxSurf(l)  = 0.0;
    shadedWall.QflxSurf(l)  = 0.0;
```

#### `UrbanSurfaceRunoffGetters.cpp`

Five 1D-double getters for `QflxSurf` on all five surfaces, plus two getters
for updated `TopH2OSoiLiq` on roof and impervious road.

### Step 4 — Fortran bindings in `include/urban_mod.F90`

For each of the 12 new functions (4 setters + 1 compute + 7 getters), add:
1. A private `_C` `interface` block binding to the C symbol
2. A public Fortran wrapper subroutine

### Step 5 — Create `biogeophys/UrbanxxSurfaceRunoffMod.F90`

Three public subroutines:

#### `urbanxx_surfaceRunoff_init(num_urbanl)`

Allocate module-level persistent 1D `c_double` buffers:
- Inputs: `wtfact`, `fover`, `frost_table`, `zwt_perched` (size `num_urbanl`)
- Outputs: `out_qflx_surf_roof`, `out_qflx_surf_imperv`, `out_qflx_surf_perv`,
  `out_qflx_surf_sunwall`, `out_qflx_surf_shadewall` (size `num_urbanl`)  
- Outputs: `out_top_h2osoi_liq_roof`, `out_top_h2osoi_liq_imperv` (size `num_urbanl`)

#### `urbanxx_surfaceRunoff(num_urbanl, num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, dtime)`

1. Loop `filter_urbanc`, pack `wtfact`, `fover` (`fover` is gridcell-indexed: use `g = col_pp%gridcell(c)` to look up `soilhydrology_vars%fover(g)`), `frost_table`, `zwt_perched` into pervious-road index-ordered buffers.
2. Call setters: `UrbanSetWtfactPerviousRoad`, `UrbanSetFoverPerviousRoad`, `UrbanSetFrostTablePerviousRoad`, `UrbanSetZwtPerchedPerviousRoad`.
   - Note: `TopH2OSoiLiq` for roof and imperv road is already set upstream (by `urbanxx_soilFluxes`), so no re-upload needed.
   - Note: `zwt` is already set by `urbanxx_soilWater` upstream.
   - Note: `ForcRain` is already set by `urbanxx_SetAtmosphericForcing` upstream.
3. `call UrbanComputeSurfaceRunoff(urbanxx, dtime, status)`
4. Retrieve: call all 7 getters into output buffers.
5. Scatter outputs back to ELM:
   - Loop `filter_urbanc` and write `qflx_surf(c)` from the appropriate per-surface buffer.
   - For roof/imperv columns also write `h2osoi_liq(c,1)` from `out_top_h2osoi_liq_roof/imperv`.
   - `qflx_surf(c) = 0` for sunlit/shaded wall (filled by URBANxx getter, but writing back is sufficient confirmation).

#### `urbanxx_surfaceRunoff_check(num_urbanl, num_urbanc, filter_urbanc)`

Compare ELM `qflx_surf(c)` and (for roof/imperv) `h2osoi_liq(c,1)` against the
contents of the output buffers.  Log max absolute difference per surface type.

### Step 6 — Call site in `biogeophys/HydrologyNoDrainageMod.F90`

Add `use UrbanxxSurfaceRunoffMod` and insert immediately after the existing
`SurfaceRunoff` call:

```fortran
call SurfaceRunoff(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
     soilhydrology_vars, soilstate_vars, dtime)
call urbanxx_surfaceRunoff(num_urbanl, num_urbanc, filter_urbanc, &
     soilhydrology_vars, soilstate_vars, dtime)
call urbanxx_surfaceRunoff_check(num_urbanl, num_urbanc, filter_urbanc)
```

### Step 7 — Initialization in `UrbanxxMod.F90`

Call `urbanxx_surfaceRunoff_init(num_urbanl)` alongside the other `_init` calls
(e.g., next to `urbanxx_soilWater_init`).

### Step 8 — CMakeLists.txt

Register the 3 new `.cpp` files in the URBANxx library target:
- `UrbanSurfaceRunoffSetters.cpp`
- `UrbanSurfaceRunoff.cpp`
- `UrbanSurfaceRunoffGetters.cpp`

---

## Files to create or modify

| File | Action |
|---|---|
| `include/private/UrbanSurfaceTypeImpl.h` | Add `QflxSurf` to all 4 structs; add `Wtfact`, `Fover`, `FrostTable`, `ZwtPerched` to `PerviousRoadDataType` |
| `include/Urban.h` | Declare 4 setters, 1 compute, 7 getters |
| `for_e3sm/include/Urban.h` | Same declarations (kept in sync) |
| `include/urban_mod.F90` | Add Fortran bindings for all 12 new functions |
| `src/UrbanSurfaceRunoffSetters.cpp` | **New** — 4 setter implementations |
| `src/UrbanSurfaceRunoff.cpp` | **New** — `UrbanComputeSurfaceRunoff` kernel |
| `src/UrbanSurfaceRunoffGetters.cpp` | **New** — 7 getter implementations |
| `CMakeLists.txt` | Register the 3 new `.cpp` files |
| `biogeophys/UrbanxxSurfaceRunoffMod.F90` | **New** — `_init`, main compute, `_check` |
| `biogeophys/HydrologyNoDrainageMod.F90` | Add `use` clause + 2 call lines after `SurfaceRunoff` |
| `biogeophys/UrbanxxMod.F90` (init module) | Add `urbanxx_surfaceRunoff_init` call |
