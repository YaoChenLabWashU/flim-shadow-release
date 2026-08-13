# FLIM shadow fitter

**Nothing of yours is modified. Not one lab file is touched, edited, or moved.
`rmpath` is a complete uninstall — there is nothing to undo, no backup to
restore, no setting to change back.**

It works by MATLAB path shadowing: two files here have the same names as two lab
fit functions, so when this folder is first on the path, MATLAB calls these
instead. Take the folder off the path and the lab code runs again, unchanged.

## What it does

Replaces the lab's double-exponential fit with a Poisson maximum-likelihood /
variable-projection fit that constrains the amplitudes to be non-negative. The
lab's Gauss-Newton least-squares fit can return negative amplitudes, which is
where unphysical population fractions come from; this one cannot.

## Install

```matlab
addpath('/path/to/flim-shadow-release', '-begin')
```

`-begin` is required — it is what puts this folder ahead of the lab tree.

## Uninstall

```matlab
rmpath('/path/to/flim-shadow-release')
```

That is the entire uninstall.

## Verify it is active — both checks, every session

```matlab
which spc_fitexp2gaussGY -all     % this folder must be line 1
which spc_fitexp2prfGY   -all     % this folder must be line 1
```

Use `-all`. Plain `which` shows only the first hit and will look correct even
when it isn't.

Then **after any fit, look for a console line starting `>>> SHADOW`**:

```
>>> SHADOW spc_fitexp2prfGY | Poisson-MLE/VARPRO | chan=1 | irf=session | pop1=0.6651 redchisq=10.6573 | src=/path/to/flim-shadow-release/spc_fitexp2prfGY.m <<<
```

**No marker means the lab fitter ran, not this one.** The line prints on failed
fits too, so silence is never ambiguous. The lab functions print nothing like it.

## Troubleshooting

If the shadow is not line 1:

```matlab
restoredefaultpath
% re-add the lab tree (zFLIM fit, utilities, display, calcs, guis), then:
addpath('/path/to/flim-shadow-release', '-begin')   % last, and with -begin
```

**A path with several copies of the zFLIM tree on it can silently win over the
shadow even when plain `which` looks correct.** There are commonly eight or more
copies of these fit functions across a machine (the acquisition tree, the FLiP
tree, publication archives). `which` reports only the first; a second copy
sitting between the shadow and the point of use is invisible to it. This cost an
hour to diagnose once — use `-all`, and use `restoredefaultpath` rather than
adding more paths on top.

## What changes vs. the lab fit

- **Amplitudes are constrained non-negative**, so `pop1` is always in [0,1]. The
  lab fit is unconstrained and frequently returns `pop1` outside it.
- **`nFreeParams` is written as 5** (4 amplitudes + a floated delta), where the
  lab writes `sum(floats)`. `spc_drawfit` divides by `N - nFreeParams`, so the
  displayed `redchisq` sits on a slightly different scale — about **1.3%** — from
  historical values for the same file. Not a fit-quality difference.
- **Both the gauss and prf menu items now use the measured PRF.** The
  GaussW / `beta6` field is read, validated, and echoed back unchanged, unused by
  the fit.
- **`frSHG` is not written**, and a stale value is removed. It is unread
  everywhere in the lab codebase — computing one here would produce a number that
  looks comparable to historical values but isn't.
- **`beta5` carries the gauss-convention delta under both shadows.** The two lab
  fitters store opposite-signed deltas in this field; both shadows agree with each
  other instead. **Un-shadow BOTH together, never one** — removing only one lets
  the remaining lab fitter reinterpret the sign silently, which propagates into
  `avgTauTrunc`, `figOffset`, and every ROI lifetime. Re-fit any file whose
  `beta5` a shadow wrote.

## What is NOT covered

**Fit1 and both single-exponential menu items still run the lab's original
code**, with nothing on screen indicating the difference. Only the two
double-exponential fitters (`spc_fitexp2gaussGY`, `spc_fitexp2prfGY`) are
shadowed. If you click Fit1 you get the lab's fit; the `>>> SHADOW` marker is how
you tell.

## Evidence

On a fair comparison against the lab's `spc_fitexp2prfGY` across 505 files, each
with its own stored lifetimes, run headless through the real lab dispatch:

- **482/505** better by the lab's own displayed Pearson metric
- **498/505** better by Poisson deviance
- **73/73** on AKAR, by both metrics

Ask Matt for the full analysis, the per-file CSVs, and the caveats — they matter
more than the headline.

## Requirements

Base MATLAB only. No toolboxes. Verified with `requiredFilesAndProducts`: these
three files plus six lab functions already on your path are the complete closure.
