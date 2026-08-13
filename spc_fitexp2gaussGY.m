function betahat = spc_fitexp2gaussGY(chan)
%SPC_FITEXP2GAUSSGY  Path-shadowing replacement for the lab's double-exponential
%   Gaussian fit. Redirects that ONE fit to the Poisson-MLE / variable-projection
%   core (flim_fit_core) while leaving every lab file untouched.
%
%   INSTALL: put this directory AHEAD of .../zFLIM/fit on the MATLAB path.
%       addpath('<repo>/labshadow', '-begin');
%   UNINSTALL: rmpath it. Nothing else changes.
%
%   Verify which one is live:  which spc_fitexp2gaussGY
%
%   ---------------------------------------------------------------------------
%   WHAT THIS COVERS, AND WHAT IT DOES NOT
%   Shadowing this name captures the double-exponential Gaussian path only:
%     - menu Fitting > double exp + gauss          (spc_main.m:211)
%     - the Fit2 button when lastFitFunction contains 'gauss'
%     - the else-fallbacks at spc_main.m:296/473/496 and spc_autoDuringAcq.m:10
%     - Next/Prev with 'fit each time' when the stored name is exactly this one
%   It does NOT cover: the Fit1 button / single-exp-gauss menu item
%   (spc_fitexpgaussGY), or either PRF menu item (spc_fitexpprfGY,
%   spc_fitexp2prfGY). NOTE: if anyone runs a PRF fit once, lastFitFunction
%   becomes a prf name and BOTH Fit buttons dispatch to the PRF fitters from then
%   on -- persistently, because it is saved into spc_backup.mat. This shadow is
%   then bypassed with no visible indication except the menu checkmark.
%
%   ===========================================================================
%   *** beta5 CARRIES THE GAUSS CONVENTION, REGARDLESS OF FUNCTION NAME ***
%   ===========================================================================
%   This shadow writes its fitted delta straight out as beta5, in this family's
%   native convention. Its sibling shadow spc_fitexp2prfGY writes THE SAME
%   convention, deliberately, so the two shadows agree with each other and a file
%   fitted by one can be re-fitted by the other with no sign flip.
%
%   The real lab pair does NOT have that property. spc_exp2gaussGY uses (x-tau_d)
%   -- +delta moves the response LATER -- while spc_exp2prfGY uses
%   interp1(x, prf1, x+deltapeak) -- +delta moves it EARLIER. The two lab fitters
%   store OPPOSITE-SIGNED numbers in the same field. Verified numerically.
%
%   CONSEQUENCE OF UN-SHADOWING. Remove these shadows and let the real
%   spc_fitexp2prfGY read a beta5 that a shadow wrote, and the sign is reinterpreted
%   silently -- no warning, no failed fit -- propagating into
%   spc_fitPostCalcs:34 (tmax) -> avgTauTrunc -> average -> spc_adjustTauOffset's
%   figOffset -> every ROI lifetime. Un-shadow BOTH fitters together, and re-fit
%   any file whose beta5 was written by a shadow. Do not un-shadow one.
%
%   ---------------------------------------------------------------------------
%   NAME. mfilename must report 'spc_fitexp2gaussGY'. Four separate lab functions
%   branch on this string: spc_fitFuncName (findstr 'gauss'), and
%   spc_fitParamsFromGlobal / spc_betaIntoGlobal / spc_fitPostCalcs (findstr
%   'prf', which decides whether beta6 is a Gaussian width or a scatter
%   fraction). Renaming this file changes lab behaviour. Do not.
%
%   ---------------------------------------------------------------------------
%   FIT WINDOW -- DELIBERATE OFF-BY-ONE. This uses the lab's own convention,
%   taken straight from spc_fitParamsFromGlobal:
%       range = round([fitstart fitend]/nsPerPoint);   % NO +1
%       lifetime = spc.lifetimes{chan}(range(1):range(2));
%   fit_flim_file uses lo = round(fitstart/dt)+1 instead, so the two differ by one
%   bin. That one bin sits exactly at the rising edge, where delta is determined.
%   The lab convention is used here ON PURPOSE so that (a) the decay handed to the
%   fitter is the same vector the lab has always fitted, and (b) the curve handed
%   to spc_drawfit lines up with the time axis spc_drawfit builds from `range`.
%   DO NOT "fix" this to match fit_flim_file: it would silently shift every
%   number the rig produces relative to all previously published fits.
%
%   ---------------------------------------------------------------------------
%   beta6. In measured-PRF mode there is no Gaussian IRF width to report -- the
%   IRF is the PRF. beta6 is therefore NOT synthesised. The stored value is read,
%   validated, and written back unchanged; only an unusable value (NaN, Inf,
%   non-scalar, or outside spc_initialValue_double's repair range of [0,4) ns) is
%   replaced, with the same 0.11 ns that spc_initialValue_double would have used.
%   NaN is repaired here explicitly because the lab's own repair test uses
%   comparisons that NaN fails, so a NaN would otherwise propagate into the next
%   real Gaussian fit and destroy it. The staleness is recorded in
%   spc.fit.shadowRes.warning (the GUI has no field that can show it).
%
%   nargin==0. spc_auto.m:69 calls this with no arguments. The real function
%   throws there and spc_auto catches it into fit_error=1; this shadow raises
%   spc_fitexp2gaussGY:noChannel so that path behaves identically. Do not
%   "helpfully" default the channel -- spc_auto then proceeds on a stale
%   spc.fit.beta0 that nothing refreshes.
%
%   nFreeParams is written as 5 (4 amplitudes + floated delta), which is this
%   fitter's actual count. spc_drawfit divides chisq by (N - nFreeParams), so the
%   redchisq shown in the GUI is on a slightly different scale from the lab's
%   historical value for the same file. side_by_side_shadow.m reports both.

    global spc

    % spc_auto.m:69 calls this with NO arguments (a long-standing latent bug --
    % it tests spc.fit.lastFitFunction, singular, which nothing ever writes, so
    % it always takes the else branch). The real function throws there, and
    % spc_auto's try/catch turns that into fit_error=1.
    %
    % That failure is PRESERVED here on purpose. spc_auto goes on to read
    % spc.fit.beta0, which neither the real function nor this shadow refreshes,
    % so succeeding would only mean proceeding with a stale value instead of
    % failing visibly. This is the lab's dead path; it is not ours to repair.
    if nargin < 1 || isempty(chan)
        error('spc_fitexp2gaussGY:noChannel', ...
              ['spc_fitexp2gaussGY requires a channel argument. Called with none ' ...
               '(see spc_auto.m:69). The real lab function fails here too; that ' ...
               'failure is preserved deliberately.']);
    end

    betahat = [];
    note = {};

    % ---- dispatch key first, exactly as the real function does -------------
    spc.fits{chan}.lastFitFunction = mfilename;
    spc.fits{chan}.fitOrder        = 2;

    % ---- read the lab's own parameter marshalling -------------------------
    [betaInit, range, floats] = spc_fitParamsFromGlobal(chan);   %#ok<ASGLU>

    nsPerPoint = spc.datainfo.psPerUnit/1000;
    dt         = nsPerPoint;
    pulseI     = spc.datainfo.pulseInt / spc.datainfo.psPerUnit * 1000;

    yfull = spc.lifetimes{chan};
    if iscell(yfull), yfull = yfull{1}; end
    yfull = double(yfull(:));

    % lab window convention (see header). Clamp only to keep indexing legal --
    % the canonical lab function would throw on range(1)==0; one lab variant
    % patches it the same way.
    lo = range(1); hi = range(2);
    if lo < 1
        lo = 1; note{end+1} = 'range(1)<1 clamped to 1'; %#ok<AGROW>
    end
    if hi > numel(yfull)
        hi = numel(yfull); note{end+1} = 'range(2) clamped to numel(lifetime)'; %#ok<AGROW>
    end

    % ---- beta6 guard (read / validate / echo) ------------------------------
    b6 = spc.fits{chan}.beta6;                     % ns in the gauss family
    if ~(isnumeric(b6) && isscalar(b6) && isfinite(b6) && b6 >= 0 && b6 < 4)
        b6used = 0.11;                             % spc_initialValue_double's value
        note{end+1} = sprintf('beta6 unusable (%s) -> repaired to 0.11 ns', mat2str(b6)); %#ok<AGROW>
    else
        b6used = b6;
    end

    % ---- measured PRF, if this channel has a usable one --------------------
    prf = [];
    if isfield(spc.fits{chan},'prf')
        p = double(spc.fits{chan}.prf(:));
        if ~isempty(p) && all(isfinite(p)) && sum(p) > 0 && numel(p) == numel(yfull)
            prf = p;
        elseif ~isempty(p) && numel(p) ~= numel(yfull)
            note{end+1} = sprintf('prf length %d ~= lifetime length %d -> ignored', ...
                                  numel(p), numel(yfull)); %#ok<AGROW>
        end
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

    % "Fix delta" in the lab GUI means HOLD DELTA AT THE STORED beta5. The core's
    % own fixDelta pins to the data's rising edge, which is a different value, so
    % the stored beta5 (bins, window-relative -- the same convention the core uses
    % for x) is handed over as the seed. Only when the box is checked: with delta
    % floating, the rising-edge seed is left alone, because that is the seed every
    % published number was produced with.
    if fixDelta
        deltaSeed = spc.fits{chan}.beta5 / dt;
    else
        deltaSeed = [];
    end
    % floats(6) (fix_g) is not mapped: in measured-PRF mode there is no Gaussian
    % width to float. If the channel has no PRF the core falls back to the
    % analytic Gaussian, and the width stays pinned at beta6.

    % Orientation matters: the lab slices spc.lifetimes{chan}, which is a ROW,
    % and stores that row in spc.fit.lifetime and (via spc_drawfit) in .curve.
    % Slice from the original row so the shapes match the real function.
    labRow = spc.lifetimes{chan};
    if iscell(labRow), labRow = labRow{1}; end
    spc.fit.lifetime = double(labRow(lo:hi));
    spc.fit.failedFit = 0;

    % ---- run the core ------------------------------------------------------
    failed = false;
    try
        res = flim_fit_core(D, spc.fits{chan}.beta2, spc.fits{chan}.beta4, ...
                            'session', [], false, fitTau1, fitTau2, fixDelta, false, false, ...
                            deltaSeed);
    catch ME
        % An ambiguous external PRF cannot arise here (this path never supplies
        % one), but translate it and anything else into the lab's failure
        % convention rather than letting it escape into the GUI.
        if strcmp(ME.identifier,'fit_flim_file:ambiguousPRF')
            note{end+1} = ['ambiguous PRF refused: ' ME.message]; %#ok<AGROW>
        else
            note{end+1} = ['fit error: ' ME.identifier ' ' ME.message]; %#ok<AGROW>
        end
        res = [];
        failed = true;
    end

    if ~failed && ~strcmp(res.status,'ok')
        note{end+1} = ['skipped: ' res.reason]; %#ok<AGROW>
        failed = true;
    end

    if failed
        spc.fit.failedFit        = 1;
        spc.fits{chan}.failedFit = 1;
        spc.fit.shadowRes  = res;
        spc.fit.shadowNote = joincsv(note);
        beep;                                   % same audible cue as spc_nlinfitGY
        fprintf('>>> SHADOW %s | FIT FAILED | chan=%d | src=%s.m <<<\n', ...
                mfilename, chan, mfilename('fullpath'));
        disp(['   reason: ' joincsv(note)]);
        return;
    end

    % ---- record what the fit did and did not determine ---------------------
    if strcmp(res.irfMode,'session')
        note{end+1} = ['beta6 (GaussW) not used by this fit: IRF is the measured ' ...
                       'PRF. Stored value echoed unchanged.']; %#ok<AGROW>
    else
        note{end+1} = 'no usable PRF -> fell back to the analytic Gaussian IRF at beta6'; %#ok<AGROW>
    end
    if res.tauFree
        note{end+1} = ['tau1 and/or tau2 were FLOATED (lab fix-tau boxes are off): '  ...
                       'pop1/pop2 are NOT calibrated biological readouts for this fit']; %#ok<AGROW>
    end
    if res.lowQuality && ~isempty(res.warning)
        note{end+1} = res.warning; %#ok<AGROW>
    end
    res.warning = joincsv(note);
    spc.fit.shadowRes  = res;
    spc.fit.shadowNote = res.warning;

    % ---- write-back: betahat in BINS, lab parameter order ------------------
    %   [ampl1 tau1 ampl2 tau2 tau_d tau_g]
    betahat = [ res.ampl1;
                res.tau1_ns/dt;
                res.ampl2;
                res.tau2_ns/dt;
                res.myDelta_bins;
                b6used/dt ];

    spc.fits{chan}.nFreeParams = 5;     % 4 amplitudes + floated delta
    spc.fit.backCorr = res.bkg;         % this fitter's fitted background level

    spc.fit.failedFit        = 0;
    spc.fits{chan}.failedFit = 0;

    spc_betaIntoGlobal(betahat, chan, chan == spc_mainChannelChoice);
    spc_fitPostCalcs(chan);

    % ---- draw, on the lab's own time axis ----------------------------------
    t = (range(1):range(2)) * spc.datainfo.psPerUnit/1000;
    % pass the curve as a ROW, as the real function does (exp2gaussGY returns a
    % row because x is a row); spc_drawfit stores it verbatim in .curve
    spc_drawfit(t, res.fit_curve(:).', spc.fit.lifetime, betahat, chan);

    mark_floating_tau(chan, res.tauFree);
    shadow_banner(chan, res);
end

% ------------------------------------------------------------------
function mark_floating_tau(chan, tauFree)
%MARK_FLOATING_TAU  Surface a floated lifetime on the ONE per-channel status
%   widget the lab already puts on the lifetime figure: gui.spc.chisq(chan),
%   created in spc_drawInit:116 and normally set by spc_drawfit to the reduced
%   chi-squared. When either fix-tau box is off, pop1/pop2 stop being calibrated
%   biological readouts, and shadowRes.warning is invisible to a GUI user.
%
%   Self-clearing: spc_drawfit rewrites the String on every subsequent fit, so
%   the suffix disappears by itself. The ForegroundColor does NOT get rewritten
%   by spc_drawfit, so the widget's original colour is stashed once and restored
%   on any fit where the lifetimes are fixed.
%
%   No new widget, no dialog, no capability removed. Silent no-op when the GUI
%   is not present (headless comparison runs).
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
