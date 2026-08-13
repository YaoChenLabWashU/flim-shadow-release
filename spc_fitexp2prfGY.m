function betahat = spc_fitexp2prfGY(chan)
%SPC_FITEXP2PRFGY  Path-shadowing replacement for the lab's double-exponential
%   measured-PRF fit. Redirects that ONE fit to the Poisson-MLE / variable-
%   projection core (flim_fit_core) while leaving every lab file untouched.
%
%   INSTALL: put this directory AHEAD of .../zFLIM/fit on the MATLAB path.
%       addpath('<repo>/labshadow', '-begin');
%   UNINSTALL: rmpath it. Nothing else changes.
%
%   Verify which one is live:  which spc_fitexp2prfGY
%
%   ===========================================================================
%   *** beta5 CARRIES THE GAUSS CONVENTION, REGARDLESS OF THIS FUNCTION'S NAME ***
%   ===========================================================================
%   This shadow writes its fitted delta STRAIGHT OUT as beta5, in the same
%   convention the gauss shadow uses. It does NOT convert to the sign convention
%   the real spc_exp2prfGY expects. Both shadows therefore agree with each other,
%   and a file fitted by one can be re-fitted by the other with no sign flip.
%
%   The real lab pair does NOT have that property. In spc_exp2prfGY the delta is
%   applied as interp1(x, prf1, x + deltapeak) -- shifting the PRF by +delta moves
%   the response EARLIER -- whereas spc_exp2gaussGY uses (x - tau_d), where +delta
%   moves the response LATER. The two lab fitters store OPPOSITE-SIGNED numbers in
%   the same field. This was verified numerically, not inferred from the source.
%
%   CONSEQUENCE OF UN-SHADOWING. If this file is removed from the path and the
%   real spc_fitexp2prfGY runs on a file this shadow last fitted, beta5 is read
%   back with the opposite meaning. That is silent -- no warning, no failed fit --
%   and it propagates:
%     * spc_fitPostCalcs:34   tmax = (range(2)-range(1)-beta0(5))*nsPerPoint
%                             -> wrong tmax -> wrong avgTauTrunc
%     * spc_fitPostCalcs:51   average = avgTauTrunc
%     * spc_adjustTauOffset   figOffset is derived from avgTauTrunc
%     * every ROI lifetime downstream of figOffset
%   So an un-shadowed rig does not merely lose the better fitter; it can produce
%   shifted lifetimes from stored parameters. Un-shadow BOTH fitters together, and
%   re-fit any file whose beta5 was written by a shadow. Do not un-shadow one.
%
%   ---------------------------------------------------------------------------
%   WHAT THIS COVERS, AND WHAT IT DOES NOT
%   Shadowing this name captures the double-exponential PRF path only:
%     - menu Fitting > double exp + prf
%     - the Fit2 button when lastFitFunction contains 'prf'
%   EXPLICIT BOUNDARY -- these still run LAB code, unchanged:
%     - the Fit1 button and the single-exp + gauss menu item -> spc_fitexpgaussGY
%     - the single-exp + prf menu item                       -> spc_fitexpprfGY
%   Only the two double-exponential fitters are shadowed. A user who clicks Fit1
%   gets the lab's Gauss-Newton single-exponential fit with no indication that a
%   different estimator produced the Fit2 numbers beside it.
%
%   ---------------------------------------------------------------------------
%   NAME. mfilename must report 'spc_fitexp2prfGY'. Four lab functions branch on
%   this string: spc_fitFuncName (findstr 'prf'), and spc_fitParamsFromGlobal /
%   spc_betaIntoGlobal / spc_fitPostCalcs (findstr 'prf', which decides whether
%   beta6 is a scatter term or a Gaussian width, and whether beta6 is scaled by
%   nsPerPoint). Renaming this file changes lab behaviour. Do not.
%
%   ---------------------------------------------------------------------------
%   FIT WINDOW -- DELIBERATE OFF-BY-ONE, same as the gauss shadow:
%       range = round([fitstart fitend]/nsPerPoint);   % NO +1
%   fit_flim_file uses lo = round(fitstart/dt)+1, so the two differ by one bin at
%   the rising edge, where delta is determined. The lab convention is used here ON
%   PURPOSE so the decay handed to the fitter is the vector the lab has always
%   fitted, and so the curve lines up with the axis spc_drawfit builds from
%   `range`. DO NOT "fix" this.
%
%   ---------------------------------------------------------------------------
%   frSHG IS NOT WRITTEN. The real function computes it (line 41) by re-evaluating
%   the model with beta6 zeroed. This shadow's model has no lab-convention scatter
%   term to zero, so any number produced here would be computed differently from
%   every historical frSHG while looking directly comparable to them. It is
%   therefore NOT written, and a stale value left by a previous lab fit is REMOVED
%   on entry so nothing downstream can read a number belonging to a different fit.
%   This is safe: frSHG is write-only across the entire codebase (9,293 files
%   searched -- it is displayed nowhere and read by nothing).
%
%   ---------------------------------------------------------------------------
%   beta6 (scatter, in this family). This fitter has no lab-convention scatter
%   parameter -- the SHG contribution is a free amplitude solved by the inner
%   convex step, not a seeded fraction. beta6 is therefore ECHOED UNCHANGED, not
%   synthesised, and only an unusable value is repaired. NOTE the scaling
%   asymmetry: spc_betaIntoGlobal:13 multiplies beta(6) by 1 in the prf branch
%   (by nsPerPoint in the gauss branch), so the echo is placed in betahat
%   UNSCALED here, where the gauss shadow divides by dt.
%
%   KNOWN INTEROP WART (not repaired here). Files fitted by the gauss family carry
%   a beta6 that is a Gaussian WIDTH (~0.10-0.15 ns). Under a prf lastFitFunction
%   that same number is read back as a SCATTER seed. That mis-seeding is what makes
%   the real spc_fitexp2prfGY diverge on such files (pop1 far outside [0,1],
%   negative amplitudes) -- measured on 505 files. Echoing beta6 unchanged neither
%   causes nor worsens it. Writing 0 instead would actively help a subsequent LAB
%   prf fit, but it would overwrite a stored calibration field with a value this
%   fitter did not measure, so it is deliberately NOT done. One-line change if you
%   later decide otherwise.
%
%   ---------------------------------------------------------------------------
%   IRF MODE. Always the measured PRF ('session'), never the analytic Gaussian.
%   This is the PRF menu item, so a usable PRF is REQUIRED: if the channel has
%   none, this fails through the lab's failure convention rather than silently
%   falling back to a Gaussian and reporting a "prf" fit that used no PRF.
%
%   nFreeParams is written as 5 (4 amplitudes + floated delta), this fitter's
%   actual count. The real function instead writes sum(floats). spc_drawfit
%   divides chisq by (N - nFreeParams), so the displayed redchisq is on a
%   different scale from the lab's historical value for the same file.
%
%   FAILURE CONVENTION matches the gauss shadow, and DIFFERS from the real
%   function: on failure this returns BEFORE spc_betaIntoGlobal / spc_fitPostCalcs
%   / spc_drawfit, so a failed fit cannot half-update the globals or draw a curve
%   from a betahat that was never fitted. The real function calls spc_drawfit
%   unconditionally (its line 52, outside the if/else).

    global spc

    % Same no-argument guard as the gauss shadow. Nothing in the tree calls the
    % prf fitters with no argument today, but the guard keeps the two shadows
    % behaving identically if that ever changes.
    if nargin < 1 || isempty(chan)
        error('spc_fitexp2prfGY:noChannel', ...
              'spc_fitexp2prfGY requires a channel argument. Called with none.');
    end

    betahat = [];
    note = {};

    % ---- dispatch key first, exactly as the real function does -------------
    spc.fits{chan}.lastFitFunction = mfilename;
    spc.fits{chan}.fitOrder        = 2;

    % ---- frSHG: drop any stale value NOW, before anything can read it ------
    if isfield(spc.fits{chan}, 'frSHG')
        spc.fits{chan} = rmfield(spc.fits{chan}, 'frSHG');
        note{end+1} = 'stale frSHG removed; this fitter does not compute one'; %#ok<AGROW>
    end

    % ---- read the lab's own parameter marshalling -------------------------
    [betaInit, range, floats] = spc_fitParamsFromGlobal(chan);   %#ok<ASGLU>

    nsPerPoint = spc.datainfo.psPerUnit/1000;
    dt         = nsPerPoint;
    pulseI     = spc.datainfo.pulseInt / spc.datainfo.psPerUnit * 1000;

    yfull = spc.lifetimes{chan};
    if iscell(yfull), yfull = yfull{1}; end
    yfull = double(yfull(:));

    % lab window convention (see header). Clamp only to keep indexing legal.
    lo = range(1); hi = range(2);
    if lo < 1
        lo = 1; note{end+1} = 'range(1)<1 clamped to 1'; %#ok<AGROW>
    end
    if hi > numel(yfull)
        hi = numel(yfull); note{end+1} = 'range(2) clamped to numel(lifetime)'; %#ok<AGROW>
    end

    % ---- beta6 guard (read / validate / echo; see header) ------------------
    % Range test mirrors the gauss shadow's so the two agree on what "unusable"
    % means. beta6 is not used by this fit either way.
    b6 = spc.fits{chan}.beta6;
    if ~(isnumeric(b6) && isscalar(b6) && isfinite(b6) && b6 >= 0 && b6 < 4)
        b6used = 0.11;
        note{end+1} = sprintf('beta6 unusable (%s) -> repaired to 0.11', mat2str(b6)); %#ok<AGROW>
    else
        b6used = b6;
    end

    % ---- the measured PRF is MANDATORY on this path ------------------------
    prf = [];
    if isfield(spc.fits{chan},'prf')
        p = double(spc.fits{chan}.prf(:));
        if ~isempty(p) && all(isfinite(p)) && sum(p) > 0 && numel(p) == numel(yfull)
            prf = p;
        elseif ~isempty(p) && numel(p) ~= numel(yfull)
            note{end+1} = sprintf('prf length %d ~= lifetime length %d', ...
                                  numel(p), numel(yfull)); %#ok<AGROW>
        end
    end

    % ---- stage range / prf / lifetime on spc.fit, as the real function does -
    % The real function does this at its lines 18-20, BEFORE fitting, because
    % spc_exp2prfGY reads spc.fit.prf and spc.fit.range out of the global rather
    % than taking them as arguments. Nothing in this shadow's own numerics needs
    % them, but they are staged identically so that anything else reading spc.fit
    % mid-fit (or a lab function called afterwards) sees what it always has.
    spc.fit.range = range;
    if isfield(spc.fits{chan},'prf')
        spc.fit.prf = spc.fits{chan}.prf;
    end

    % Orientation matters: the lab slices spc.lifetimes{chan}, which is a ROW.
    labRow = spc.lifetimes{chan};
    if iscell(labRow), labRow = labRow{1}; end
    spc.fit.lifetime  = double(labRow(lo:hi));
    spc.fit.failedFit = 0;

    if isempty(prf)
        note{end+1} = ['no usable measured PRF on this channel -- the prf fitter ' ...
                       'will not silently fall back to an analytic Gaussian']; %#ok<AGROW>
        fail(chan, [], note);
        return;
    end

    % ---- build the in-memory data struct the core expects ------------------
    D = struct('ok',true,'reason','', ...
        'yfull',yfull,'dt',dt,'pulseI',pulseI,'prf',prf, ...
        'tau_g_b',b6used/dt,'beta6_ns',b6used, ...
        'labDelta_bins',spc.fits{chan}.beta5/dt,'labP1',getfielddef(spc.fits{chan},'pop1',NaN), ...
        'lo',lo,'hi',hi, ...
        'fitstart',spc.fits{chan}.fitstart,'fitend',spc.fits{chan}.fitend, ...
        'winReverted',false,'winNote','');

    % ---- freedoms, translated from the lab's floats mask -------------------
    %   floats = ~[0 fixtau1 0 fixtau2 fix_delta fix_g]
    fitTau1  =  floats(2);
    fitTau2  =  floats(4);
    fixDelta = ~floats(5);

    % "Fix delta" means HOLD DELTA AT THE STORED beta5. beta5 is window-relative
    % bins after /dt -- the same frame the core's seed uses, because both lab
    % fitters evaluate on x = 1:length(lifetime) (spc_fitexp2prfGY.m:13) and
    % spc_fitParamsFromGlobal:10 divides by nsPerPoint. Verified at source and
    % empirically. Only supplied when the box is checked; otherwise the core's
    % rising-edge seed is left alone.
    if fixDelta
        deltaSeed = spc.fits{chan}.beta5 / dt;
    else
        deltaSeed = [];
    end
    % floats(6) (fix_g) is not mapped: there is no lab-convention scatter
    % parameter in this model to fix or float.

    % ---- run the core ------------------------------------------------------
    failed = false;
    try
        res = flim_fit_core(D, spc.fits{chan}.beta2, spc.fits{chan}.beta4, ...
                            'session', [], false, fitTau1, fitTau2, fixDelta, false, false, ...
                            deltaSeed);
    catch ME
        note{end+1} = ['fit error: ' ME.identifier ' ' ME.message]; %#ok<AGROW>
        res = []; failed = true;
    end

    if ~failed && ~strcmp(res.status,'ok')
        note{end+1} = ['skipped: ' res.reason]; %#ok<AGROW>
        failed = true;
    end

    % A fallback to 'gaussian' cannot happen here (prf was validated above), but
    % assert it rather than trust it -- a "prf" fit that used no PRF must not pass.
    if ~failed && ~strcmp(res.irfMode,'session')
        note{end+1} = sprintf('core fell back to irfMode=%s on the prf path', res.irfMode); %#ok<AGROW>
        failed = true;
    end

    if failed
        fail(chan, res, note);
        return;
    end

    % ---- record what the fit did and did not determine ---------------------
    note{end+1} = ['beta6 (scatter) not used by this fit: the SHG term is a free ' ...
                   'amplitude, not a seeded fraction. Stored value echoed unchanged.']; %#ok<AGROW>
    if res.tauFree
        note{end+1} = ['tau1 and/or tau2 were FLOATED (lab fix-tau boxes are off): ' ...
                       'pop1/pop2 are NOT calibrated biological readouts for this fit']; %#ok<AGROW>
    end
    if res.lowQuality && ~isempty(res.warning)
        note{end+1} = res.warning; %#ok<AGROW>
    end
    res.warning = joincsv(note);
    spc.fit.shadowRes  = res;
    spc.fit.shadowNote = res.warning;

    % ---- write-back: betahat in the lab's parameter order ------------------
    %   [ampl1 tau1 ampl2 tau2 tau_d tau_g]
    % spc_betaIntoGlobal scales by [1 dt 1 dt dt 1] in the prf branch, so
    % element 6 is UNSCALED here (the gauss shadow divides by dt instead).
    % Element 5 is written in the GAUSS convention -- see the header banner.
    betahat = [ res.ampl1;
                res.tau1_ns/dt;
                res.ampl2;
                res.tau2_ns/dt;
                res.myDelta_bins;
                b6used ];

    spc.fits{chan}.nFreeParams = 5;     % 4 amplitudes + floated delta
    spc.fit.backCorr = res.bkg;

    spc.fit.failedFit        = 0;
    spc.fits{chan}.failedFit = 0;

    spc_betaIntoGlobal(betahat, chan, chan == spc_mainChannelChoice);
    spc_fitPostCalcs(chan);

    % ---- draw, on the lab's own time axis ----------------------------------
    t = (range(1):range(2)) * spc.datainfo.psPerUnit/1000;
    spc_drawfit(t, res.fit_curve(:).', spc.fit.lifetime, betahat, chan);

    mark_floating_tau(chan, res.tauFree);
    shadow_banner(chan, res);
end

% ------------------------------------------------------------------
function fail(chan, res, note)
%FAIL  The gauss shadow's failure convention: mark failed, stash the diagnosis,
%   beep, and return WITHOUT touching the globals spc_betaIntoGlobal would write.
    global spc
    spc.fit.failedFit        = 1;
    spc.fits{chan}.failedFit = 1;
    spc.fit.shadowRes  = res;
    spc.fit.shadowNote = joincsv(note);
    beep;                                   % same audible cue as spc_nlinfitGY
    fprintf('>>> SHADOW %s | FIT FAILED | chan=%d | src=%s.m <<<\n', ...
            mfilename, chan, mfilename('fullpath'));
    disp(['   reason: ' joincsv(note)]);
end

% ------------------------------------------------------------------
function mark_floating_tau(chan, tauFree)
%MARK_FLOATING_TAU  Surface a floated lifetime on gui.spc.chisq(chan), the one
%   per-channel status widget the lab already puts on the lifetime figure.
%   Identical to the gauss shadow's marker so both fitters flag it the same way.
%   Self-clearing: spc_drawfit rewrites the String on every subsequent fit.
    global gui
    MARKCOL = [0.6 0 0];
    try
        if isempty(gui) || ~isfield(gui,'spc') || ~isfield(gui.spc,'chisq'), return; end
        if numel(gui.spc.chisq) < chan || ~ishandle(gui.spc.chisq(chan)), return; end
        h = gui.spc.chisq(chan);

        c = get(h,'ForegroundColor');
        if ~isfield(gui.spc,'chisqBaseColor'), gui.spc.chisqBaseColor = {}; end
        if numel(gui.spc.chisqBaseColor) < chan || isempty(gui.spc.chisqBaseColor{chan})
            if ~isequal(c, MARKCOL)          % never stash our own marker colour
                gui.spc.chisqBaseColor{chan} = c;
            end
        end

        if tauFree
            set(h, 'String', [get(h,'String') ' t*'], 'ForegroundColor', MARKCOL);
        elseif numel(gui.spc.chisqBaseColor) >= chan && ~isempty(gui.spc.chisqBaseColor{chan})
            set(h, 'ForegroundColor', gui.spc.chisqBaseColor{chan});
        end
    catch
        % a GUI quirk must never take the fit down
    end
end

% ------------------------------------------------------------------
function shadow_banner(chan, res)
%SHADOW_BANNER  One unambiguous line per successful fit, so which fitter ran is
%   never inferred from side effects again. Prints the RESOLVED FILE PATH,
%   because several copies of this function name exist on this machine and
%   `which` shows only the first. The lab functions print nothing like this.
    p = mfilename('fullpath');
    fprintf(['>>> SHADOW %s | Poisson-MLE/VARPRO | chan=%d | irf=%s | ' ...
             'pop1=%.4f redchisq=%s | src=%s.m <<<\n'], ...
            mfilename, chan, res.irfMode, res.pop1, rcstr(chan), p);
end

function s = rcstr(chan)
    global spc
    if isfield(spc,'fits') && numel(spc.fits) >= chan && isfield(spc.fits{chan},'redchisq')
        s = sprintf('%.4f', spc.fits{chan}.redchisq);
    else
        s = 'n/a';
    end
end

% ------------------------------------------------------------------
function v = getfielddef(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end

function s = joincsv(c)
%JOINCSV  strjoin is R2013a; this file may run on an R2012 rig.
    if isempty(c), s = ''; return; end
    s = c{1};
    for k = 2:numel(c), s = [s '; ' c{k}]; end   %#ok<AGROW>
end
