# FLIM shadow fitter

Replaces the lab's double-exponential fit with a Poisson maximum-likelihood /
variable-projection fit that constrains the amplitudes to be non-negative. The
lab's Gauss-Newton least-squares fit can return negative amplitudes, which is
where unphysical population fractions come from; this one cannot.

**Nothing of yours is modified.** Not one lab file is touched, edited, or moved.
It works by MATLAB path shadowing: two files here have the same names as two lab
fit functions, so when this folder is first on the path, MATLAB calls these
instead. `rmpath` is a complete uninstall — nothing to restore, nothing to undo.

Base MATLAB only. No toolboxes.

---

## 1. Download

Go to the **Releases** page of this repository and click **Source code (zip)**
under the latest release. Unzip it anywhere you like — your Desktop is fine.

You will end up with a folder called **`flim-shadow-release-1.0`**.

## 2. Install

In MATLAB, change into that folder and add it to the path:

```matlab
cd /path/to/flim-shadow-release-1.0     % or drag the folder into the Terminal,
                                        % or copy it from Finder's Get Info
addpath(pwd, '-begin')
```

Using `cd` first and then `pwd` saves you typing the path correctly.

**`-begin` is required, not optional.** It is what puts this folder ahead of the
lab tree. Without it MATLAB keeps calling the lab's copies and nothing changes.

## 3. Verify

Run both checks:

```matlab
which spc_fitexp2gaussGY -all
which spc_fitexp2prfGY   -all
```

A **correct** result is two lines — this folder first, the lab's copy second,
marked `% Shadowed`:

```
/Users/you/Desktop/flim-shadow-release-1.0/spc_fitexp2gaussGY.m
/Users/you/lab/zFLIM/fit/spc_fitexp2gaussGY.m                     % Shadowed
```

If you see only **one** line, the lab's zFLIM tree is not on your path yet —
you haven't started the FLIM GUI. Start it, then check again. You need to see
both lines, in that order.

Use `-all`. Plain `which` shows only the first hit and looks correct even when
it isn't.

## 4. Run

Start the FLIM GUI however you normally do, open a file, and hit
**Fit with double**. Then read the console:

- A line beginning **`>>> SHADOW`** means this fitter ran:

  ```
  >>> SHADOW spc_fitexp2prfGY | Poisson-MLE/VARPRO | chan=1 | irf=session | pop1=0.6651 redchisq=10.6573 | src=/Users/you/Desktop/flim-shadow-release-1.0/spc_fitexp2prfGY.m <<<
  ```

- **No marker**, or a line reading **`fraction of photons from SHG:`**, means the
  **lab's** fitter ran, not this one. Go back to step 3.

The marker prints on failed fits too, so silence is never ambiguous.

## 5. Uninstall

```matlab
rmpath('/path/to/flim-shadow-release-1.0')
```

That is the entire uninstall. Nothing to restore.

## 6. Troubleshooting

If the shadow is not line 1:

```matlab
restoredefaultpath
% re-add the lab tree (zFLIM fit, utilities, display, calcs, guis), then:
addpath('/path/to/flim-shadow-release-1.0', '-begin')   % last, and with -begin
```

**A path with several copies of the zFLIM tree on it can silently win over the
shadow even when plain `which` looks correct.** There are commonly eight or more
copies of these fit functions across a machine (the acquisition tree, the FLiP
tree, publication archives). `which` reports only the first; a second copy
sitting between the shadow and the point of use is invisible to it. This cost an
hour to diagnose once — use `-all`.

**If your lab code lives on a network share, `restoredefaultpath` removes those
paths too**, not just the stray ones. Re-add the lab tree afterwards (or restart
the FLIM GUI, if that is what normally sets your path up) before adding this
folder with `-begin` last.

---

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

## Evidence

On a fair comparison against the lab's `spc_fitexp2prfGY` across 505 files, each
with its own stored lifetimes, run headless through the real lab dispatch:

- **482/505** better by the lab's own displayed Pearson metric
- **498/505** better by Poisson deviance
- **73/73** on AKAR, by both metrics

