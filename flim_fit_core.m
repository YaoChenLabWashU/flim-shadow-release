function res = flim_fit_core(D, tau1_ns, tau2_ns, irfMode, externalPRF, fitGaussW, fitTau1, fitTau2, fixDelta, singleExp, fitIrfW, deltaSeed)
%FLIM_FIT_CORE  Poisson-MLE delta-floating FLIM fit, fed from memory.
%
%   res = flim_fit_core(D, tau1_ns, tau2_ns, irfMode, externalPRF, fitGaussW,
%                       fitTau1, fitTau2, fixDelta, singleExp, fitIrfW)
%
%   D is the struct returned by flim_load_file, or any struct with the same
%   fields built by hand -- which is the point of the split: a caller holding the
%   decay in memory (the lab GUI's `global spc`) can fit without touching disk.
%
%   Every other argument, and the returned res, are exactly as documented in
%   fit_flim_file. This function contains the numerics; it does no file I/O.
%
%   deltaSeed (OPTIONAL, bins, window-relative). Overrides the delta the scan is
%   centred on. Omitted or empty -> the seed is the data's own rising edge, which
%   is the default path and the one every published number was produced with.
%   Supplied -> the scan centres there instead, so fixDelta=true pins delta to
%   THAT value rather than to the rising edge. Exists because the lab GUI's
%   "Fix delta" checkbox means "hold delta at the stored beta5", which is a
%   different quantity from the rising edge; the shadow passes beta5 here.
%   Leaving it out changes nothing.

    if nargin < 1 || isempty(tau1_ns), tau1_ns = 2.14; end
    if nargin < 2 || isempty(tau2_ns), tau2_ns = 0.69; end
    if nargin < 4 || isempty(irfMode),  irfMode = 'session'; end
    if nargin < 5, externalPRF = []; end
    if nargin < 6 || isempty(fitGaussW), fitGaussW = false; end
    if nargin < 7 || isempty(fitTau1),   fitTau1 = false; end
    if nargin < 8 || isempty(fitTau2),   fitTau2 = false; end
    if nargin < 9 || isempty(fixDelta), fixDelta = false; end
    if nargin < 10 || isempty(singleExp), singleExp = false; end
    if nargin < 11 || isempty(fitIrfW),   fitIrfW = false; end
    if nargin < 12, deltaSeed = []; end
    if singleExp, fitTau2 = false; end   % no second lifetime to float
    if islogical(irfMode)               % back-compat with old useMeasuredIRF flag
        if irfMode, irfMode = 'session'; else, irfMode = 'gaussian'; end
    end
    irfMode = lower(char(irfMode));

    res = struct('status','skipped','reason','', ...
        'myP1',NaN,'myDelta_bins',NaN,'chi2_myDelta',NaN,'chi2_labConv',NaN,'chi2_labDelta',NaN, ...
        'pop1',NaN,'pop2',NaN,'ampl1',NaN,'ampl2',NaN,'amplSHG',NaN,'avgTau_ns',NaN, ...
        'labP1',NaN,'labDelta_bins',NaN, ...
        'x',[],'ywin',[],'fit_curve',[],'residual',[], ...
        'dGrid',[],'nll_prof',[],'p1_prof',[], ...
        'dt',NaN,'pulseI',NaN,'tau1_ns',tau1_ns,'tau2_ns',tau2_ns, ...
        'fitstart',NaN,'fitend',NaN,'irfMode','','tauG_ns',NaN,'fitGaussW',false, ...
        'tauFree',false,'singleExp',singleExp,'irfBroad_ns',NaN,'fitIrfW',false,'beta6_ns',NaN, ...
        'winReverted',false,'winNote','','lowQuality',false,'warning','','bkg',NaN);

    % Fields the loader establishes. Copied BEFORE the ok test because the
    % original set them before its own early returns, and callers read them.
    res.beta6_ns      = D.beta6_ns;
    res.labDelta_bins = D.labDelta_bins;
    res.labP1         = D.labP1;
    if ~D.ok, res.reason = D.reason; return; end
    res.dt = D.dt; res.pulseI = D.pulseI;
    res.winReverted = D.winReverted; res.winNote = D.winNote;
    res.fitstart = D.fitstart; res.fitend = D.fitend;

    dt = D.dt; pulseI = D.pulseI; yfull = D.yfull;
    lo = D.lo; hi = D.hi; prf = D.prf; tau_g_b = D.tau_g_b;

    tau1_b  = tau1_ns/dt;
    tau2_b  = tau2_ns/dt;

    if hi-lo < 20, res.reason = 'window too small / off range'; return; end
    ywin = yfull(lo:hi); Nwin = numel(ywin); x = (1:Nwin)';
    res.x = x; res.ywin = ywin;

    [~, pk] = max(ywin);
    halfmax = (max(ywin)+min(ywin))/2;
    edgeIdx = find(ywin(1:pk) >= halfmax, 1, 'first');
    if isempty(edgeIdx), edgeIdx = pk; end
    seed_b = edgeIdx;
    % Optional override (see header). Absent/empty/non-finite -> rising edge, i.e.
    % the default path is untouched.
    if ~isempty(deltaSeed) && isnumeric(deltaSeed) && isscalar(deltaSeed) && isfinite(deltaSeed)
        seed_b = deltaSeed;
    end

    % ----- choose the IRF model -----
    % Resolve which measured PRF (if any) to use for 'session'/'external'.
    irfPRF = [];
    if strcmp(irfMode,'session')
        irfPRF = prf;                              % from spcSave.fits.prf
    elseif strcmp(irfMode,'external')
        irfPRF = load_external_prf(externalPRF);   % vector or path to .mat with 'prf'
    end
    prfOK = ~isempty(irfPRF) && sum(irfPRF) > 0 && all(isfinite(irfPRF));
    if any(strcmp(irfMode,{'session','external'})) && ~prfOK
        irfMode = 'gaussian';                      % no usable PRF -> fall back
    end
    res.irfMode = irfMode;

    % SHG is a single sharp 1-2 bin spike at the delta in ALL modes (a linear-interp
    % unit impulse) -- not collinear with the broad decays, so no SHG/decay degeneracy.
    shgFun = @(td) max(0, 1 - abs(x - td));

    % IRF Gaussian width: fixed at beta6 unless fitGaussW (and gaussian mode), in
    % which case it is floated via a 1D grid nested inside the delta scan. tau1/tau2
    % stay fixed throughout, so the decay basis only depends on (delta, tau_g).
    % The nested 'tg' parameter means different things per IRF mode: in 'gaussian'
    % mode it is the IRF Gaussian width; in measured ('session'/'external') mode it
    % is an EXTRA Gaussian broadening sigma convolved onto the measured PRF (sigma=0
    % -> the raw PRF, i.e. unchanged behavior). Either may be floated via the grid
    % nested in the delta scan; tau1/tau2 stay fixed so pop1 stays a valid readout.
    % fitIrfW (measured modes only): false -> no broadening (baseline); true -> FLOAT
    % sigma over a grid; a positive NUMERIC value -> FIXED sigma in ns (no fitting),
    % converted to bins per-file via dt so it is fixed in TIME across rigs.
    % fitGaussW (gaussian mode only): false -> fix width at the file's beta6; true ->
    % FLOAT the width; a positive NUMERIC value -> FIX the width at that value (ns),
    % overriding beta6 (mirrors the lab GUI's editable GaussW field).
    measMode = any(strcmp(irfMode,{'session','external'}));
    fixedIW = isnumeric(fitIrfW) && ~islogical(fitIrfW) && isscalar(fitIrfW) && fitIrfW > 0 && measMode;
    fixedGW = isnumeric(fitGaussW) && ~islogical(fitGaussW) && isscalar(fitGaussW) && fitGaussW > 0 && strcmp(irfMode,'gaussian');
    floatGW = islogical(fitGaussW) && fitGaussW && strcmp(irfMode,'gaussian');
    floatIW = islogical(fitIrfW) && fitIrfW && measMode;
    res.fitGaussW = floatGW;
    res.fitIrfW   = floatIW;
    if floatGW
        tgGrid = tau_g_b * linspace(0.4, 2.2, 19);  % includes 1.0x -> never worse than fixed
    elseif fixedGW
        tgGrid = fitGaussW / dt;                     % user-set Gaussian width (ns -> bins)
    elseif floatIW
        tgGrid = [0, linspace(0.2, 3.0, 15)];       % extra PRF broadening sigma (bins); 0 = raw PRF
    elseif fixedIW
        tgGrid = fitIrfW / dt;                       % FIXED broadening sigma (ns -> bins)
    elseif measMode
        tgGrid = 0;                                 % measured IRF, no extra broadening
    else
        tgGrid = tau_g_b;                           % fixed gaussian width (file's beta6)
    end

    % --- decide which lifetimes float, optimize them if requested ---
    % make_basis(t1,t2) builds basisFun(td,tg) for the current IRF mode, so the
    % same code path serves the fixed-tau fit and every tau candidate during the
    % nested floating-tau search below.
    res.tauFree = fitTau1 || fitTau2;
    tauLo = 0.05/dt; tauHi = 5.0/dt;          % bin-domain clamp (~lab's 0.01-5 ns)

    if res.tauFree
        % Float tau1 and/or tau2 by an outer search wrapped around the delta(+tg)
        % scan: J(t1,t2) = min over delta,tg of the Poisson NLL. The surface is
        % non-convex and symmetric under the tau1<->tau2 swap, so a bounded
        % log-spaced multi-start grid lands in the right basin first, then
        % fminsearch polishes. Delta uses a coarse scan inside the search and the
        % full fine scan once for the reported fit.
        if fixDelta, dHalfOpt = 0; dStepOpt = 1; else, dHalfOpt = 10; dStepOpt = 0.5; end
        [tau1_b, tau2_b] = optimize_taus(ywin, irfMode, irfPRF, x, pulseI, Nwin, ...
            shgFun, tgGrid, seed_b, tau1_b, tau2_b, fitTau1, fitTau2, ...
            tauLo, tauHi, dHalfOpt, dStepOpt, singleExp);
    end
    basisFun = make_basis(irfMode, irfPRF, x, pulseI, Nwin, shgFun, tau1_b, tau2_b, singleExp);

    % delta is pinned to the rising-edge seed when fixDelta (single-point scan),
    % otherwise scanned/floated as usual.
    if fixDelta, dHalf = 0; else, dHalf = 12; end
    [myDelta_b, myTg_b, dGrid, nll_prof, p1_prof, ok] = scan_delta( ...
        ywin, basisFun, tgGrid, seed_b, dHalf, 0.1);
    if ~fixDelta && ok && (abs(myDelta_b-(seed_b-12))<0.2 || abs(myDelta_b-(seed_b+12))<0.2)
        [myDelta_b, myTg_b, dGrid, nll_prof, p1_prof, ok] = scan_delta( ...
            ywin, basisFun, tgGrid, myDelta_b, 20, 0.1);
    end
    if ~ok, res.reason = 'delta scan failed'; return; end
    res.dGrid = dGrid; res.nll_prof = nll_prof; res.p1_prof = p1_prof;
    res.myDelta_bins = myDelta_b;
    if floatIW || fixedIW
        res.irfBroad_ns = myTg_b * dt;   % extra PRF broadening sigma (ns), fitted or fixed
    else
        res.tauG_ns = myTg_b * dt;       % gaussian-mode IRF width (ns)
    end

    Bm = basisFun(myDelta_b, myTg_b);

    % ---- flag a rank-deficient / degenerate amplitude fit ----
    % The basis is [decay1(, decay2), SHG, background]. It loses column rank when
    % the SHG spike falls outside the fit window (all-zero column), the two decay
    % columns collapse together (tau1~=tau2, or a broad IRF smearing them over a
    % short window), or the window holds almost no signal. The amplitude solve is
    % then ill-posed -- this is the "Rank deficient, rank = 2" warning from the
    % B\y step in fit_amps -- and pop1/avgTau/fit_curve are unreliable. Mark the
    % result low-quality (status stays 'ok') with a diagnosed reason so the GUI can
    % surface it rather than displaying garbage as if it were a clean fit.
    tolR = max(size(Bm)) * eps(norm(Bm,1));
    if rank(Bm, tolR) < size(Bm,2)
        res.lowQuality = true;
        reasons = {};
        if max(abs(Bm(:,end-1))) < 1e-9      % SHG column (2nd-to-last) is ~all zeros
            reasons{end+1} = 'SHG spike outside window';
        end
        if ~singleExp                         % two decay columns collinear?
            d1 = Bm(:,1); d2 = Bm(:,2);
            if norm(d1) > 0 && norm(d2) > 0 && ...
               abs(d1.'*d2)/(norm(d1)*norm(d2)) > 1 - 1e-6
                reasons{end+1} = 'decay components collapsed';
            end
        end
        if nnz(ywin) < 20                      % almost no signal in the window
            reasons{end+1} = 'near-empty window';
        end
        if isempty(reasons), detail = '';
        else
            % explicit join instead of strjoin (R2013a); produces the identical char row
            joined = reasons{1};
            for rk = 2:numel(reasons)
                joined = [joined ', ' reasons{rk}];   %#ok<AGROW>
            end
            detail = [' (' joined ')'];
        end
        res.warning = ['rank-deficient fit' detail];
    end

    am = fit_amps(ywin, Bm, zeros(Nwin,1));
    % Background amplitude: the last basis column is ones(N,1), so am(end) is the
    % constant the model sits on. Reported because the lab records an equivalent
    % quantity (spc.fit.backCorr) and a caller writing back into `spc` otherwise
    % has to leave a stale value there. Additive only -- nothing above uses it.
    res.bkg = am(end);
    mum = Bm*am;
    res.fit_curve = mum;
    res.residual  = ywin - mum;
    if singleExp
        % basis is [decay1, SHG, background]: one lifetime, no population fraction
        res.ampl1 = am(1); res.amplSHG = am(2); res.ampl2 = NaN;
        res.pop1 = NaN; res.pop2 = NaN; res.myP1 = NaN;
        res.tau1_ns = tau1_b*dt; res.tau2_ns = NaN;
        res.avgTau_ns = res.tau1_ns;       % a single exponential's <tau> is just tau
    else
        res.ampl1 = am(1); res.ampl2 = am(2); res.amplSHG = am(3);
        res.pop1  = am(1)/(am(1)+am(2));
        res.pop2  = am(2)/(am(1)+am(2));
        res.myP1  = res.pop1;
        % report the lifetimes actually used (fitted when floated, else the inputs)
        res.tau1_ns = tau1_b*dt; res.tau2_ns = tau2_b*dt;
        res.avgTau_ns = res.pop1*res.tau1_ns + res.pop2*res.tau2_ns;
    end
    ncolFit = size(Bm,2);                  % 4 (double) or 3 (single) free amplitudes
    rc = @(mu) sum((ywin-mu).^2 ./ max(ywin,1)) / (Nwin-ncolFit);
    res.chi2_myDelta = rc(mum);
    % lab-convention reduced chi2 (Pearson: divide by model, dof = N-2),
    % matching spc_drawfit.m's redchisq = sum((y-fit).^2./max(fit,1))/(N-nFreeParams)
    res.chi2_labConv = sum((ywin - mum).^2 ./ max(mum,1)) / (Nwin - 2);

    if ~isnan(res.labDelta_bins)
        Bl = basisFun(res.labDelta_bins, myTg_b);
        al = fit_amps(ywin, Bl, zeros(Nwin,1));
        res.chi2_labDelta = rc(Bl*al);
    end

    res.status = 'ok';

end

function basisFun = make_basis(irfMode, irfPRF, x, pulseI, Nwin, shgFun, tau1_b, tau2_b, singleExp)
%MAKE_BASIS  Build basisFun(td,tg) = [decay1(, decay2), SHG, background] for the
%   given lifetimes (in bins). Gaussian mode uses the analytic erfc decay; the
%   measured-IRF modes numerically convolve each pure exp (plus pre-pulse wrap)
%   with the PRF. singleExp -> only the tau1 decay column (no decay2). Rebuilt per
%   tau candidate during the floating-tau search.
    if nargin < 9 || isempty(singleExp), singleExp = false; end
    if strcmp(irfMode,'gaussian')
        if singleExp
            decayFun = @(td,tg) expgauss_comp(1, tau1_b, td, tg, x, pulseI);
        else
            decayFun = @(td,tg) [ expgauss_comp(1, tau1_b, td, tg, x, pulseI), ...
                                  expgauss_comp(1, tau2_b, td, tg, x, pulseI) ];
        end
    else   % 'session' or 'external': numerically convolve each decay with the PRF
        prf_n = irfPRF(:) / sum(irfPRF);  L = numel(prf_n);  [~, p0] = max(prf_n);
        % Single-pulse measured response on a grid long enough to reach the tails of
        % previous pulses that wrap under the window (n out to Nwin + 2 periods).
        Npad  = ceil(Nwin + 2*pulseI + L);
        nL    = (0:Npad-1)';  gridL = (1:Npad)';
        sp1_0 = conv(exp(-nL/tau1_b), prf_n);  sp1_0 = sp1_0(1:Npad);
        if singleExp, sp2_0 = []; else, sp2_0 = conv(exp(-nL/tau2_b), prf_n); sp2_0 = sp2_0(1:Npad); end
        % tg = extra Gaussian broadening sigma (bins) tests "the measured PRF is a
        % little too narrow for the fluorescence". By associativity the broadening
        % can act on the finished single-pulse response: exp (x) PRF (x) Gauss(tg) =
        % (exp (x) PRF) (x) Gauss(tg). tg=0 returns the raw response (unchanged fit).
        % Cached per tg so each grid sigma convolves once, not once per delta.
        cache = containers.Map('KeyType','double','ValueType','any');
        % Periodic steady state: current pulse (k=0, at delta) PLUS previous pulses
        % (k=1,2 at delta-k*pulseI) whose exp tails form the pre-rise floor. p0 aligns
        % delta to the IRF peak (same seed/scan as gaussian mode).
        decayFun = @(td,tg) meas_decay(cache, sp1_0, sp2_0, tg, td, Npad, gridL, x, p0, pulseI, singleExp);
    end
    basisFun = @(td,tg) [ decayFun(td,tg), shgFun(td), ones(Nwin,1) ];
end

function [s1, s2] = broaden_resp(cache, sp1_0, sp2_0, tg, Npad)
%BROADEN_RESP  Single-pulse response(s) broadened by an extra Gaussian sigma tg
%   (bins). tg<=0 -> the raw responses. The symmetric kernel is centre-aligned so
%   the broadening adds NO net time shift (delta semantics unchanged). Cached per tg.
    if tg <= 0, s1 = sp1_0; s2 = sp2_0; return; end
    key = round(tg*1e4);
    if isKey(cache,key), v = cache(key); s1 = v{1}; s2 = v{2}; return; end
    w  = max(1, ceil(4*tg));  m = (-w:w)';  gk = exp(-m.^2/(2*tg^2));  gk = gk/sum(gk);
    c1 = conv(sp1_0, gk);  s1 = c1(w+1 : w+Npad);
    if isempty(sp2_0), s2 = [];
    else, c2 = conv(sp2_0, gk);  s2 = c2(w+1 : w+Npad); end
    cache(key) = {s1, s2};
end

function D = meas_decay(cache, sp1_0, sp2_0, tg, td, Npad, gridL, x, p0, pulseI, singleExp)
%MEAS_DECAY  Periodic-wrap decay column(s) for the measured-IRF modes at delta td
%   and broadening tg. Sums the current pulse plus two previous-pulse wraps.
    [s1, s2] = broaden_resp(cache, sp1_0, sp2_0, tg, Npad);
    wrap = @(c) interp1(gridL, c, x - td + p0,            'linear', 0) ...
              + interp1(gridL, c, x - td + p0 + pulseI,   'linear', 0) ...
              + interp1(gridL, c, x - td + p0 + 2*pulseI, 'linear', 0);
    if singleExp, D = wrap(s1); else, D = [wrap(s1), wrap(s2)]; end
end

function [t1b, t2b] = optimize_taus(ywin, irfMode, irfPRF, x, pulseI, Nwin, shgFun, ...
        tgGrid, seed_b, t1_0, t2_0, freeT1, freeT2, tauLo, tauHi, dHalf, dStep, singleExp)
%OPTIMIZE_TAUS  Minimize the delta(+tg)-profiled Poisson NLL over the free
%   lifetime(s). Bounded log-spaced multi-start grid (avoids the tau-swap local
%   minimum and other non-convex basins) followed by an fminsearch polish. Fixed
%   taus pass through unchanged. Returns lifetimes in bins, long component as t1b.
    fixed  = [t1_0 t2_0];
    mask   = logical([freeT1 freeT2]);
    clampp = @(p) min(max(p, tauLo), tauHi);
    nll = @(p) tau_profile_nll(clampp(p), fixed, mask, ywin, irfMode, irfPRF, ...
                               x, pulseI, Nwin, shgFun, tgGrid, seed_b, dHalf, dStep, singleExp);

    % --- bounded multi-start grid over the free dims ---
    cand = exp(linspace(log(tauLo), log(tauHi), 8));   % log-spaced candidates
    p0 = fixed(mask);                                  % current taus as one start
    if nnz(mask) == 2
        [A,B] = ndgrid(cand, cand);
        starts = [A(:) B(:)];
        starts = starts(starts(:,1) >= starts(:,2), :); % canonical half: t1 >= t2
    else
        starts = cand(:);
    end
    starts = [p0(:)'; starts];

    bestJ = inf; bestP = p0;
    for i = 1:size(starts,1)
        J = nll(starts(i,:));
        if J < bestJ, bestJ = J; bestP = starts(i,:); end
    end

    % --- local polish from the best grid point ---
    opt  = optimset('Display','off','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',200);
    pPol = fminsearch(nll, bestP, opt);
    if nll(pPol) <= bestJ, bestP = pPol; end
    bestP = clampp(bestP);

    t = fixed; t(mask) = bestP;
    t1b = t(1); t2b = t(2);
    if freeT1 && freeT2 && t1b < t2b, tmp = t1b; t1b = t2b; t2b = tmp; end
end

function J = tau_profile_nll(p, fixed, mask, ywin, irfMode, irfPRF, x, pulseI, ...
        Nwin, shgFun, tgGrid, seed_b, dHalf, dStep, singleExp)
%TAU_PROFILE_NLL  Objective for optimize_taus: NLL of the best delta(+tg) fit at
%   the given free lifetimes. Returns a large finite value on invalid taus so the
%   optimizer simply avoids them.
    t = fixed; t(mask) = p;
    if any(t <= 0), J = 1e18; return; end
    bf = make_basis(irfMode, irfPRF, x, pulseI, Nwin, shgFun, t(1), t(2), singleExp);
    [~, ~, ~, nll_prof] = scan_delta(ywin, bf, tgGrid, seed_b, dHalf, dStep);
    J = min(nll_prof);
    if ~isfinite(J), J = 1e18; end
end

function [bestDelta, bestTg, dGrid, nll_prof, p1_prof, ok] = scan_delta( ...
        ywin, basisFun, tgGrid, center_b, half, stepb)
    % Scan delta; at each delta, nest a 1D search over the IRF width tg (tgGrid is
    % a single value when the width is fixed). Profiles store the best-over-tg value.
    dGrid = (center_b-half) : stepb : (center_b+half);
    nd = numel(dGrid);
    nll_prof = inf(nd,1); p1_prof = nan(nd,1); tg_at = nan(nd,1);
    for k = 1:nd
        bestn = inf; bp1 = NaN; btg = tgGrid(1);
        for tg = tgGrid(:)'
            Bk = basisFun(dGrid(k), tg);
            [ak, nv] = fit_amps(ywin, Bk, zeros(numel(ywin),1));
            if nv < bestn
                bestn = nv; bp1 = ak(1)/(ak(1)+ak(2)); btg = tg;
            end
        end
        nll_prof(k) = bestn; p1_prof(k) = bp1; tg_at(k) = btg;
    end
    [mn,km] = min(nll_prof); ok = isfinite(mn);
    bestDelta = dGrid(km); bestTg = tg_at(km);
end

function p = load_external_prf(externalPRF)
%LOAD_EXTERNAL_PRF  Resolve the 'external' IRF: a numeric vector, or a path to a
%   .mat file containing a variable 'prf' (or its first numeric vector). [] if none.
    p = [];
    if isempty(externalPRF), return; end
    if isnumeric(externalPRF)
        p = double(externalPRF(:)); return;
    end
    % isa(..,'string') instead of isstring (R2016b): isa returns false for an
    % unknown class name, so on R2012 this is exactly ischar(), while on modern
    % releases it still accepts a string scalar. Same branch taken on both.
    if ischar(externalPRF) || isa(externalPRF, 'string')
        try
            E = load(char(externalPRF));
        catch
            return;
        end
        if isfield(E,'prf')
            p = double(E.prf(:));
        else
            % No variable named 'prf'. Collect every numeric-vector candidate.
            % sort() first so the candidate list is in a stable order regardless
            % of how the .mat stored its variables (MAT-file variable ordering is
            % not guaranteed and varies by release and file format).
            fn = sort(fieldnames(E));
            cand = {};
            for k = 1:numel(fn)
                v = E.(fn{k});
                if isnumeric(v) && isvector(v) && numel(v) > 1
                    cand{end+1} = fn{k};   %#ok<AGROW>
                end
            end
            if numel(cand) > 1
                % Ordering alone only makes the guess repeatable, not right.
                % With several candidates there is no basis for choosing, and
                % silently picking one means silently fitting with the wrong
                % IRF -- so refuse and name them.
                nameList = cand{1};
                for k = 2:numel(cand)
                    nameList = [nameList ', ' cand{k}];   %#ok<AGROW>
                end
                error('fit_flim_file:ambiguousPRF', ...
                      ['External PRF file has no variable named ''prf'' and holds ' ...
                       '%d numeric vectors (%s). Refusing to guess which is the IRF. ' ...
                       'Save the intended one as ''prf'', or pass the vector directly.'], ...
                      numel(cand), nameList);
            elseif numel(cand) == 1
                p = double(E.(cand{1})(:));
            end
            % zero candidates -> p stays [], and the caller falls back to the
            % gaussian IRF exactly as before
        end
    end
end

function h = expgauss_comp(A, tau, tau_d, tau_g, x, pulseI)
    main = A*exp(tau_g^2/(2*tau^2) - (x-tau_d)/tau) .* ...
           erfc((tau_g^2 - tau*(x-tau_d))./(sqrt(2)*tau*tau_g));
    pre  = A*exp(tau_g^2/(2*tau^2) - (x-tau_d+pulseI)/tau) .* ...
           erfc((tau_g^2 - tau*(x-tau_d+pulseI))./(sqrt(2)*tau*tau_g));
    h = (main + pre)/2;
end

function [a, nllval] = fit_amps(y, B, contam)
    % Lower bounds per column. Basis is [decay1(, decay2), SHG, background]; the
    % SHG coefficient (always the second-to-last column -- col 3 of 4 in double-exp,
    % col 2 of 3 in single-exp) may be negative, all others must stay >=0.
    ncol = size(B,2);
    lb = zeros(ncol,1);
    if ncol >= 3, lb(ncol-1) = -Inf; end
    pos = isfinite(lb);                 % the non-negative columns
    clamp = @(v) max(v, lb);

    % A degenerate basis makes B\y rank-deficient; silence the console warning
    % here (fit_amps runs hundreds of times during the delta scan) -- the caller
    % detects the condition on the final basis and surfaces it via res.warning.
    ws = warning('off','MATLAB:rankDeficientMatrix');
    a = clamp( B \ max(y - contam, 0) );
    warning(ws);
    a(pos) = max(a(pos), 1);            % keep clamped columns strictly positive to start
    tol = 1e-7; maxit = 80;
    for it = 1:maxit
        [f_old, g, H] = nll(a, y, B, contam);
        if any(~isfinite(g)) || any(~isfinite(H(:))), break; end
        % pinv (not \) for the Newton step: H is PSD but can be near-singular when
        % the IRF spike and IRF(x)decay columns are collinear (measured-IRF mode);
        % pinv gives the minimum-norm descent step without singular-matrix warnings.
        step = -pinv(H)*g; s = 1;
        while true
            a_try = clamp(a + s*step); mu_try = B*a_try + contam;
            if all(mu_try > 0) && nll(a_try, y, B, contam) <= f_old, break; end
            s = s/2; if s < 1e-10, break; end
        end
        a_new = clamp(a + s*step);
        if norm(a_new - a) < tol*max(norm(a),1), a = a_new; break; end
        a = a_new;
    end
    nllval = nll(a, y, B, contam);
end

function [f, g, H] = nll(a, y, B, contam)
    mu = max(B*a + contam, 1e-9);
    f  = sum(mu - y .* log(mu));
    g  = B' * (1 - y ./ mu);
    % bsxfun (not implicit expansion): (N x 1) .* (N x M) requires R2016b. bsxfun
    % computes the identical elementwise product with the identical operation order,
    % so H is bit-for-bit unchanged. R2012-compatible.
    H  = B' * bsxfun(@times, y ./ mu.^2, B);
end
