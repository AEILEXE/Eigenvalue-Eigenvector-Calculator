'use strict';

/* =========================================================================
   Eigenvalue & Eigenvector Calculator
   All linear algebra runs client-side (no server, no external library).

   Method
   -------
   1. Reduce the input matrix to upper Hessenberg form (real Householder
      reflections) -- a similarity transform, so eigenvalues are preserved.
   2. Run the shifted QR algorithm (Wilkinson shift, complex arithmetic) with
      deflation on the Hessenberg matrix to obtain all eigenvalues, real or
      complex. This is the same family of algorithm production numerical
      libraries (LAPACK, etc.) use, and is far better conditioned for
      repeated/clustered eigenvalues than forming the characteristic
      polynomial and finding its roots.
   3. Cluster numerically-equal eigenvalues to obtain algebraic multiplicity.
   4. For each distinct eigenvalue, compute a basis for the null space of
      (A - lambda*I) via complex Gaussian elimination -- this yields the
      eigenvector(s) for that eigenvalue (their count is the geometric
      multiplicity, which may be less than the algebraic multiplicity for
      defective matrices, matching how the original SymPy program behaves).
   ========================================================================= */

(function () {
  /* ----------------------------- Complex numbers ----------------------------- */
  function cx(re, im) { return { re: re, im: im || 0 }; }
  function cAdd(a, b) { return cx(a.re + b.re, a.im + b.im); }
  function cSub(a, b) { return cx(a.re - b.re, a.im - b.im); }
  function cMul(a, b) { return cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
  function cDiv(a, b) {
    var d = b.re * b.re + b.im * b.im;
    if (d === 0) d = 1e-300;
    return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
  }
  function cAbs(a) { return Math.hypot(a.re, a.im); }
  function cNeg(a) { return cx(-a.re, -a.im); }
  function cConj(a) { return cx(a.re, -a.im); }
  function cScaleR(a, s) { return cx(a.re * s, a.im * s); }
  function cSqrt(z) {
    var r = cAbs(z);
    if (r === 0) return cx(0, 0);
    var re = Math.sqrt((r + z.re) / 2);
    var im = Math.sqrt(Math.max(0, r - z.re) / 2);
    if (z.im < 0) im = -im;
    return cx(re, im);
  }

  /* ----------------------------- Complex matrices ----------------------------- */
  function realToComplexMatrix(A, n) {
    return A.map(function (row) { return row.map(function (v) { return cx(v, 0); }); });
  }

  /* ----------------------------- Hessenberg reduction ----------------------------- */
  function toHessenberg(realA, n) {
    var A = realA.map(function (row) { return row.slice(); });
    for (var k = 0; k < n - 2; k++) {
      var L = n - k - 1;
      var x = [];
      for (var i = k + 1; i < n; i++) x.push(A[i][k]);
      var normx = Math.sqrt(x.reduce(function (s, v) { return s + v * v; }, 0));
      if (normx < 1e-300) continue;
      var alpha = x[0] >= 0 ? -normx : normx;
      var v = x.slice();
      v[0] -= alpha;
      var vnorm2 = v.reduce(function (s, vv) { return s + vv * vv; }, 0);
      if (vnorm2 < 1e-300) continue;
      var beta = 2 / vnorm2;

      for (var j = 0; j < n; j++) {
        var s = 0;
        for (var l = 0; l < L; l++) s += v[l] * A[k + 1 + l][j];
        s *= beta;
        for (l = 0; l < L; l++) A[k + 1 + l][j] -= v[l] * s;
      }
      for (var r = 0; r < n; r++) {
        var s2 = 0;
        for (l = 0; l < L; l++) s2 += A[r][k + 1 + l] * v[l];
        s2 *= beta;
        for (l = 0; l < L; l++) A[r][k + 1 + l] -= s2 * v[l];
      }
    }
    return A;
  }

  /* ----------------------------- Hessenberg-preserving QR step -----------------------------
     One implicit-shift QR step on an m x m upper-Hessenberg matrix, done with
     m-1 Givens rotations (O(m^2)) instead of a dense Householder QR
     decomposition + matrix multiply (O(m^3)). A QR step on a Hessenberg
     matrix produces another Hessenberg matrix, so this is mathematically
     the same shifted-QR iteration as a dense implementation -- only the
     cost of computing each step changes, not what is computed. This keeps
     large matrices (up to 100x100) practical: a dense-QR-per-step approach
     took ~16.5s end-to-end for a 100x100 matrix; this reduces that to ~3s. */
  function complexGivens(a, b) {
    if (cAbs(b) < 1e-300) return { c: 1, s: cx(0, 0) };
    if (cAbs(a) < 1e-300) return { c: 0, s: cx(1, 0) };
    var absa = cAbs(a);
    var norm = Math.sqrt(absa * absa + cAbs(b) * cAbs(b));
    var c = absa / norm;
    var phaseA = cx(a.re / absa, a.im / absa);
    var s = cScaleR(cMul(phaseA, cConj(b)), 1 / norm);
    return { c: c, s: s };
  }

  function hessenbergQRStep(M, m, shift) {
    var i, j, k, row, col;
    for (i = 0; i < m; i++) M[i][i] = cSub(M[i][i], shift);

    var rotations = [];
    for (k = 0; k < m - 1; k++) {
      var g = complexGivens(M[k][k], M[k + 1][k]);
      rotations.push(g);
      for (col = k; col < m; col++) {
        var top = M[k][col], bot = M[k + 1][col];
        M[k][col] = cAdd(cMul(cx(g.c, 0), top), cMul(g.s, bot));
        M[k + 1][col] = cAdd(cMul(cNeg(cConj(g.s)), top), cMul(cx(g.c, 0), bot));
      }
      M[k + 1][k] = cx(0, 0);
    }

    for (k = 0; k < m - 1; k++) {
      var gg = rotations[k];
      for (row = 0; row < m; row++) {
        var left = M[row][k], right = M[row][k + 1];
        M[row][k] = cAdd(cMul(left, cx(gg.c, 0)), cMul(right, cConj(gg.s)));
        M[row][k + 1] = cAdd(cMul(left, cNeg(gg.s)), cMul(right, cx(gg.c, 0)));
      }
    }

    for (i = 0; i < m; i++) M[i][i] = cAdd(M[i][i], shift);
  }

  /* ----------------------------- Shifted QR eigenvalue solver ----------------------------- */
  function eigenvaluesQR(realA, n) {
    if (n === 1) return [cx(realA[0][0], 0)];
    var hess = toHessenberg(realA, n);
    var A = realToComplexMatrix(hess, n);
    var eigenvalues = [];
    var m = n;
    var maxTotalIter = 100 * n + 500;
    var iterCount = 0;

    while (m > 1 && iterCount < maxTotalIter) {
      var stagnant = 0;
      var prevSub = Infinity;
      while (m > 1) {
        iterCount++;
        if (iterCount > maxTotalIter) break;
        var sub = cAbs(A[m - 1][m - 2]);
        var scale = cAbs(A[m - 2][m - 2]) + cAbs(A[m - 1][m - 1]);
        if (scale < 1e-300) scale = 1;
        if (sub < 1e-13 * scale) break;

        var a = A[m - 2][m - 2], b = A[m - 2][m - 1], c = A[m - 1][m - 2], d = A[m - 1][m - 1];
        var tr = cAdd(a, d);
        var det = cSub(cMul(a, d), cMul(b, c));
        var disc = cSqrt(cSub(cMul(tr, tr), cScaleR(det, 4)));
        var e1 = cScaleR(cAdd(tr, disc), 0.5);
        var e2 = cScaleR(cSub(tr, disc), 0.5);
        var shift = cAbs(cSub(e1, d)) < cAbs(cSub(e2, d)) ? e1 : e2;

        if (Math.abs(sub - prevSub) < 1e-15 * scale) {
          stagnant++;
          if (stagnant > 10) {
            shift = cAdd(shift, cx(scale * 0.1, scale * 0.13));
            stagnant = 0;
          }
        }
        prevSub = sub;

        var subM = [];
        for (var i = 0; i < m; i++) subM.push(A[i].slice(0, m));
        hessenbergQRStep(subM, m, shift);
        for (i = 0; i < m; i++) for (var j = 0; j < m; j++) A[i][j] = subM[i][j];
      }
      eigenvalues.push(A[m - 1][m - 1]);
      m--;
    }
    if (m === 1) eigenvalues.push(A[0][0]);
    return eigenvalues;
  }

  /* ----------------------------- Grouping (algebraic multiplicity) ----------------------------- */
  function internalClean(x, scale) {
    var tol = 1e-6 * Math.max(1, scale || 1);
    if (Math.abs(x) < tol) return 0;
    var r = Math.round(x);
    if (Math.abs(x - r) < tol) return r;
    return Math.round(x * 1e6) / 1e6;
  }

  function groupEigenvalues(roots) {
    var scale = 1;
    roots.forEach(function (r) { scale = Math.max(scale, Math.abs(r.re), Math.abs(r.im)); });
    var clusterTol = Math.max(1e-4, 1e-5 * scale);
    var used = new Array(roots.length).fill(false);
    var groups = [];
    for (var i = 0; i < roots.length; i++) {
      if (used[i]) continue;
      var members = [roots[i]];
      used[i] = true;
      for (var j = i + 1; j < roots.length; j++) {
        if (used[j]) continue;
        var dist = Math.hypot(roots[i].re - roots[j].re, roots[i].im - roots[j].im);
        if (dist < clusterTol) { members.push(roots[j]); used[j] = true; }
      }
      var avgRe = members.reduce(function (s, m) { return s + m.re; }, 0) / members.length;
      var avgIm = members.reduce(function (s, m) { return s + m.im; }, 0) / members.length;
      groups.push({
        value: cx(internalClean(avgRe, scale), internalClean(avgIm, scale)),
        rawValue: cx(avgRe, avgIm),
        multiplicity: members.length
      });
    }
    groups.sort(function (a, b) {
      if (a.value.re !== b.value.re) return a.value.re - b.value.re;
      return a.value.im - b.value.im;
    });
    return groups;
  }

  /* ----------------------------- Eigenvector null space ----------------------------- */
  function nullspace(A, n, lambda) {
    var B = [];
    for (var i = 0; i < n; i++) {
      B.push([]);
      for (var j = 0; j < n; j++) {
        var v = cx(A[i][j], 0);
        if (i === j) v = cSub(v, lambda);
        B[i].push(v);
      }
    }
    var maxEntry = 0;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) maxEntry = Math.max(maxEntry, cAbs(B[i][j]));
    var eps = 1e-6 * Math.max(1, maxEntry);

    var pivotRow = 0;
    var pivotCols = [];
    var pivotRowForCol = {};
    for (var col = 0; col < n && pivotRow < n; col++) {
      var maxVal = eps;
      var maxR = -1;
      for (var r = pivotRow; r < n; r++) {
        var m = cAbs(B[r][col]);
        if (m > maxVal) { maxVal = m; maxR = r; }
      }
      if (maxR === -1) continue;
      var tmp = B[pivotRow]; B[pivotRow] = B[maxR]; B[maxR] = tmp;

      var pivotVal = B[pivotRow][col];
      for (var c2 = 0; c2 < n; c2++) B[pivotRow][c2] = cDiv(B[pivotRow][c2], pivotVal);

      for (var r2 = 0; r2 < n; r2++) {
        if (r2 === pivotRow) continue;
        var factor = B[r2][col];
        if (cAbs(factor) < 1e-300) continue;
        for (var c3 = 0; c3 < n; c3++) {
          B[r2][c3] = cSub(B[r2][c3], cMul(factor, B[pivotRow][c3]));
        }
      }
      pivotCols.push(col);
      pivotRowForCol[col] = pivotRow;
      pivotRow++;
    }

    var freeCols = [];
    for (col = 0; col < n; col++) if (pivotCols.indexOf(col) === -1) freeCols.push(col);

    var basis = [];
    freeCols.forEach(function (f) {
      var v = new Array(n).fill(null).map(function () { return cx(0, 0); });
      v[f] = cx(1, 0);
      pivotCols.forEach(function (pc) {
        var row = pivotRowForCol[pc];
        v[pc] = cNeg(B[row][f]);
      });
      basis.push(v);
    });
    return basis;
  }

  /* ----------------------------- Public entry point ----------------------------- */
  function computeEigen(matrix, n) {
    if (n === 1) {
      var val = matrix[0][0];
      return [{ value: cx(internalClean(val, Math.abs(val) || 1), 0), multiplicity: 1, vectors: [[cx(1, 0)]] }];
    }
    var roots = eigenvaluesQR(matrix, n);
    var groups = groupEigenvalues(roots);
    groups.forEach(function (g) {
      g.vectors = nullspace(matrix, n, g.rawValue);
    });
    return groups;
  }

  /* ==========================================================================
     Number / complex formatting for display
     ========================================================================== */
  function formatReal(x) {
    if (Math.abs(x) < 1e-9) x = 0;
    var rounded = Math.round(x * 1e6) / 1e6;
    if (Object.is(rounded, -0)) rounded = 0;
    var str = rounded.toFixed(6);
    str = str.replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
    if (str === '' || str === '-0') str = '0';
    return str;
  }

  function formatComplex(c) {
    var re = Math.abs(c.re) < 1e-9 ? 0 : c.re;
    var im = Math.abs(c.im) < 1e-9 ? 0 : c.im;
    if (im === 0) return formatReal(re);
    var imAbs = Math.abs(im);
    var imStr = Math.abs(imAbs - 1) < 1e-9 ? 'i' : formatReal(imAbs) + 'i';
    if (re === 0) return (im < 0 ? '-' : '') + imStr;
    return formatReal(re) + (im < 0 ? ' - ' : ' + ') + imStr;
  }

  /* ==========================================================================
     UI wiring
     ========================================================================== */
  var MIN_N = 1;
  var MAX_N = 100;
  var LARGE_MATRIX_WARN_THRESHOLD = 50;

  var sizeInput = document.getElementById('matrix-size');
  var createMatrixBtn = document.getElementById('create-matrix-btn');
  var matrixContainer = document.getElementById('matrix-container');
  var matrixSection = document.getElementById('matrix-section');
  var calculateBtn = document.getElementById('calculate-btn');
  var clearBtn = document.getElementById('clear-btn');
  var exampleBtn = document.getElementById('example-btn');
  var copyBtn = document.getElementById('copy-btn');
  var errorBox = document.getElementById('error-box');
  var largeMatrixNotice = document.getElementById('large-matrix-notice');
  var resultsSection = document.getElementById('results-section');
  var inputMatrixDisplay = document.getElementById('input-matrix-display');
  var eigenResultsContainer = document.getElementById('eigen-results');
  var calculatingIndicator = document.getElementById('calculating-indicator');

  var graphWrap = document.getElementById('graph-wrap');
  var graphHint = document.getElementById('graph-hint');
  var graphUnavailable = document.getElementById('graph-unavailable');
  var graphContainer = document.getElementById('graph-container');
  var graphInfo = document.getElementById('graph-info');
  var graphVertexInfo = document.getElementById('graph-vertex-info');
  var graphDensityNote = document.getElementById('graph-density-note');
  var graphWeightedBadge = document.getElementById('graph-weighted-badge');
  var copyGraphBtn = document.getElementById('copy-graph-btn');
  var downloadGraphBtn = document.getElementById('download-graph-btn');

  var currentN = 0;
  var lastResults = null; // holds { matrix, n, groups } for copy functionality
  var lastGraph = null; // holds { n, edges, type, weighted } for the graph feature

  function clamp(v, lo, hi) { return Math.min(hi, Math.max(lo, v)); }

  function showError(message) {
    errorBox.textContent = message;
    errorBox.hidden = !message;
  }

  function clearResults() {
    resultsSection.hidden = true;
    eigenResultsContainer.innerHTML = '';
    inputMatrixDisplay.innerHTML = '';
    lastResults = null;
    graphWrap.hidden = true;
    graphUnavailable.hidden = true;
    graphContainer.innerHTML = '';
    graphContainer.style.aspectRatio = '';
    graphInfo.innerHTML = '';
    graphVertexInfo.textContent = '';
    graphDensityNote.hidden = true;
    graphWeightedBadge.hidden = true;
    lastGraph = null;
  }

  /* ---------------- Matrix grid generation ---------------- */
  function buildMatrixGrid(n, values) {
    matrixContainer.innerHTML = '';
    if (!n || n < MIN_N) return;

    var wrapper = document.createElement('div');
    wrapper.className = 'matrix-brackets';

    var bracketLeft = document.createElement('div');
    bracketLeft.className = 'bracket bracket-left';
    bracketLeft.setAttribute('aria-hidden', 'true');

    var grid = document.createElement('div');
    grid.className = 'matrix-grid';
    grid.style.gridTemplateColumns = 'repeat(' + n + ', minmax(56px, 1fr))';
    grid.setAttribute('role', 'group');
    grid.setAttribute('aria-label', n + ' by ' + n + ' matrix entry grid');

    var bracketRight = document.createElement('div');
    bracketRight.className = 'bracket bracket-right';
    bracketRight.setAttribute('aria-hidden', 'true');

    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        var input = document.createElement('input');
        input.type = 'text';
        input.inputMode = 'decimal';
        input.autocomplete = 'off';
        input.spellcheck = false;
        input.className = 'matrix-cell';
        input.id = 'cell-' + i + '-' + j;
        input.dataset.row = i;
        input.dataset.col = j;
        input.setAttribute('aria-label', 'Row ' + (i + 1) + ', column ' + (j + 1));
        if (values && values[i] && values[i][j] !== undefined) {
          input.value = values[i][j];
        }
        input.addEventListener('keydown', onCellKeydown);
        input.addEventListener('input', onCellInput);
        grid.appendChild(input);
      }
    }

    wrapper.appendChild(bracketLeft);
    wrapper.appendChild(grid);
    wrapper.appendChild(bracketRight);
    matrixContainer.appendChild(wrapper);
  }

  function onCellInput() {
    showError('');
  }

  function getCell(row, col) {
    return document.getElementById('cell-' + row + '-' + col);
  }

  function onCellKeydown(e) {
    var row = parseInt(e.target.dataset.row, 10);
    var col = parseInt(e.target.dataset.col, 10);
    var n = currentN;
    var target = null;

    if (e.key === 'Enter') {
      e.preventDefault();
      if (col < n - 1) target = getCell(row, col + 1);
      else if (row < n - 1) target = getCell(row + 1, 0);
      else { calculateBtn.click(); return; }
    } else if (e.key === 'ArrowRight' && caretAtEnd(e.target)) {
      target = getCell(row, col + 1);
    } else if (e.key === 'ArrowLeft' && caretAtStart(e.target)) {
      target = getCell(row, col - 1);
    } else if (e.key === 'ArrowDown') {
      target = getCell(row + 1, col);
    } else if (e.key === 'ArrowUp') {
      target = getCell(row - 1, col);
    }

    if (target) {
      e.preventDefault();
      target.focus();
      target.select();
    }
  }

  function caretAtEnd(input) {
    return input.selectionStart === input.value.length;
  }
  function caretAtStart(input) {
    return input.selectionStart === 0;
  }

  /* ---------------- Matrix size handling ---------------- */
  function updateLargeMatrixNotice(n) {
    if (n >= LARGE_MATRIX_WARN_THRESHOLD) {
      largeMatrixNotice.textContent = 'Large matrix detected. Calculations may take longer depending on your device.';
      largeMatrixNotice.hidden = false;
    } else {
      largeMatrixNotice.textContent = '';
      largeMatrixNotice.hidden = true;
    }
  }

  function createMatrix() {
    var raw = sizeInput.value.trim();
    var n = parseInt(raw, 10);

    if (raw === '' || !/^\d+$/.test(raw) || isNaN(n) || n < MIN_N || n > MAX_N) {
      showError('Matrix size must be a whole number between ' + MIN_N + ' and ' + MAX_N + '.');
      currentN = 0;
      matrixContainer.innerHTML = '';
      matrixSection.hidden = true;
      updateLargeMatrixNotice(0);
      clearResults();
      return;
    }

    showError('');
    currentN = n;
    matrixSection.hidden = false;
    updateLargeMatrixNotice(n);
    buildMatrixGrid(n, null);
    clearResults();
    var first = getCell(0, 0);
    if (first) first.focus();
  }

  createMatrixBtn.addEventListener('click', createMatrix);

  sizeInput.addEventListener('keydown', function (e) {
    if (e.key === 'Enter') {
      e.preventDefault();
      createMatrix();
    }
  });

  /* ---------------- Validation & reading the matrix ---------------- */
  var NUMBER_PATTERN = /^-?(\d+\.?\d*|\.\d+)$/;

  function readMatrix() {
    var n = currentN;
    if (!n) {
      showError('Please enter a matrix size to begin.');
      return null;
    }
    var matrix = [];
    for (var i = 0; i < n; i++) {
      var row = [];
      for (var j = 0; j < n; j++) {
        var cell = getCell(i, j);
        var raw = cell.value.trim();
        if (raw === '') {
          showError('Please enter a value in every matrix cell.');
          cell.focus();
          return null;
        }
        if (!NUMBER_PATTERN.test(raw)) {
          showError('Please enter valid numerical values (e.g. 3, -2.5, 0.75).');
          cell.focus();
          return null;
        }
        var num = parseFloat(raw);
        if (!isFinite(num)) {
          showError('Please enter valid numerical values (e.g. 3, -2.5, 0.75).');
          cell.focus();
          return null;
        }
        row.push(num);
      }
      matrix.push(row);
    }
    return matrix;
  }

  /* ---------------- Rendering results ---------------- */
  function renderStaticMatrix(container, matrix, n) {
    container.innerHTML = '';
    var wrapper = document.createElement('div');
    wrapper.className = 'matrix-brackets';

    var bracketLeft = document.createElement('div');
    bracketLeft.className = 'bracket bracket-left';
    bracketLeft.setAttribute('aria-hidden', 'true');

    var grid = document.createElement('div');
    grid.className = 'matrix-grid matrix-grid--display';
    grid.style.gridTemplateColumns = 'repeat(' + n + ', minmax(48px, 1fr))';

    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        var cellEl = document.createElement('div');
        cellEl.className = 'matrix-value';
        cellEl.textContent = formatReal(matrix[i][j]);
        grid.appendChild(cellEl);
      }
    }

    var bracketRight = document.createElement('div');
    bracketRight.className = 'bracket bracket-right';
    bracketRight.setAttribute('aria-hidden', 'true');

    wrapper.appendChild(bracketLeft);
    wrapper.appendChild(grid);
    wrapper.appendChild(bracketRight);
    container.appendChild(wrapper);

    var caption = document.createElement('p');
    caption.className = 'matrix-dimension-label';
    caption.textContent = n + ' × ' + n + ' matrix';
    container.appendChild(caption);
  }

  function renderVectorBlock(vector) {
    var n = vector.length;
    var wrapper = document.createElement('div');
    wrapper.className = 'matrix-brackets matrix-brackets--vector';

    var bracketLeft = document.createElement('div');
    bracketLeft.className = 'bracket bracket-left';
    bracketLeft.setAttribute('aria-hidden', 'true');

    var grid = document.createElement('div');
    grid.className = 'vector-grid';

    for (var i = 0; i < n; i++) {
      var cellEl = document.createElement('div');
      cellEl.className = 'matrix-value';
      cellEl.textContent = formatComplex(vector[i]);
      grid.appendChild(cellEl);
    }

    var bracketRight = document.createElement('div');
    bracketRight.className = 'bracket bracket-right';
    bracketRight.setAttribute('aria-hidden', 'true');

    wrapper.appendChild(bracketLeft);
    wrapper.appendChild(grid);
    wrapper.appendChild(bracketRight);
    return wrapper;
  }

  function renderResults(matrix, n, groups) {
    renderStaticMatrix(inputMatrixDisplay, matrix, n);

    eigenResultsContainer.innerHTML = '';
    groups.forEach(function (g, idx) {
      var card = document.createElement('article');
      card.className = 'eigen-card';
      card.setAttribute('aria-label', 'Eigenvalue ' + (idx + 1));

      var header = document.createElement('div');
      header.className = 'eigen-card-header';

      var lambda = document.createElement('div');
      lambda.className = 'eigen-lambda';
      lambda.innerHTML = '<span class="eigen-lambda-symbol">λ</span><span class="eigen-lambda-eq">=</span><span class="eigen-lambda-value">' + escapeHtml(formatComplex(g.value)) + '</span>';

      var isDefective = g.vectors.length < g.multiplicity;

      var algMult = document.createElement('div');
      algMult.className = 'eigen-multiplicity' + (g.multiplicity > 1 ? ' eigen-multiplicity--repeated' : '');
      algMult.innerHTML = 'Algebraic Multiplicity: <strong>' + g.multiplicity + '</strong>' + (g.multiplicity > 1 ? ' <span class="status-badge status-badge--repeated">Repeated</span>' : '');

      var geoMult = document.createElement('div');
      geoMult.className = 'eigen-multiplicity' + (isDefective ? ' eigen-multiplicity--defective' : '');
      geoMult.innerHTML = 'Geometric Multiplicity: <strong>' + g.vectors.length + '</strong>' + (isDefective ? ' <span class="status-badge status-badge--defective">Defective</span>' : '');

      header.appendChild(lambda);
      header.appendChild(algMult);
      header.appendChild(geoMult);
      card.appendChild(header);

      var vectorsSection = document.createElement('div');
      vectorsSection.className = 'eigen-vectors';

      var vecLabel = document.createElement('p');
      vecLabel.className = 'eigen-vectors-label';
      vecLabel.textContent = g.vectors.length > 1 ? 'Eigenvectors (basis of eigenspace):' : 'Eigenvector:';
      vectorsSection.appendChild(vecLabel);

      var vecRow = document.createElement('div');
      vecRow.className = 'eigen-vectors-row';
      g.vectors.forEach(function (v) {
        vecRow.appendChild(renderVectorBlock(v));
      });
      vectorsSection.appendChild(vecRow);

      if (isDefective) {
        var note = document.createElement('p');
        note.className = 'eigen-note';
        note.textContent = 'This eigenvalue is defective: its algebraic multiplicity (' + g.multiplicity + ') is greater than its geometric multiplicity (' + g.vectors.length + '), so it has fewer independent eigenvectors than its multiplicity.';
        vectorsSection.appendChild(note);
      }

      card.appendChild(vectorsSection);
      eigenResultsContainer.appendChild(card);
    });

    resultsSection.hidden = false;
    lastResults = { matrix: matrix, n: n, groups: groups };
    if (typeof resultsSection.scrollIntoView === 'function') {
      resultsSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  }

  function escapeHtml(str) {
    var div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
  }

  /* ==========================================================================
     Laplacian Matrix -> Graph Visualization
     This feature is intentionally self-contained: it reads the same matrix
     the eigenvalue engine reads, but never touches or depends on eigenvalue
     results, so it cannot change the mathematical output above it.
     ========================================================================== */
  var LAPLACIAN_TOL = 1e-6;
  var SVG_NS = 'http://www.w3.org/2000/svg';

  /* A matrix L is a valid (weighted or unweighted) simple-graph Laplacian
     (L = D - A) when it is symmetric, every off-diagonal entry is zero or
     negative (a negative entry's magnitude is the edge weight), every
     diagonal entry is non-negative, and every row sums to zero. An
     unweighted Laplacian is simply the special case where every edge
     weight equals 1 -- so this single analysis handles both, and the
     original K4-style 0/-1 matrices validate exactly as before. Returns
     either { ok: true, edges: [[i, j, weight], ...], weighted: bool } or
     { ok: false, reason: <human-readable explanation> } -- never silently
     "fixes" or reinterprets an invalid matrix. */
  function laplacianAnalysis(matrix, n) {
    if (!n || n < 1) return { ok: false, reason: 'the matrix is empty.' };
    var i, j;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (!isFinite(matrix[i][j])) return { ok: false, reason: 'the matrix contains a non-finite value.' };
      }
    }
    for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
        if (Math.abs(matrix[i][j] - matrix[j][i]) > LAPLACIAN_TOL) {
          return {
            ok: false,
            reason: 'the matrix is not symmetric (L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
              ' but L[' + (j + 1) + '][' + (i + 1) + '] = ' + formatReal(matrix[j][i]) + '), so it cannot represent an undirected graph.'
          };
        }
      }
    }
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        if (i === j) continue;
        if (matrix[i][j] > LAPLACIAN_TOL) {
          return {
            ok: false,
            reason: 'entry L[' + (i + 1) + '][' + (j + 1) + '] = ' + formatReal(matrix[i][j]) +
              ' is positive; off-diagonal entries of a Laplacian must be zero or negative (a negative value represents an edge weight).'
          };
        }
      }
      if (matrix[i][i] < -LAPLACIAN_TOL) {
        return {
          ok: false,
          reason: 'diagonal entry L[' + (i + 1) + '][' + (i + 1) + '] is negative; Laplacian diagonal entries must be non-negative.'
        };
      }
    }
    for (i = 0; i < n; i++) {
      var rowSum = 0;
      for (j = 0; j < n; j++) rowSum += matrix[i][j];
      if (Math.abs(rowSum) > LAPLACIAN_TOL * Math.max(1, n)) {
        return {
          ok: false,
          reason: 'row ' + (i + 1) + ' sums to ' + formatReal(rowSum) + ', not 0; every row of a Laplacian matrix must sum to zero.'
        };
      }
    }

    var edges = [];
    var weighted = false;
    for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
        var val = matrix[i][j];
        if (val < -LAPLACIAN_TOL) {
          var w = -val;
          edges.push([i + 1, j + 1, w]);
          if (Math.abs(w - 1) > LAPLACIAN_TOL) weighted = true;
        }
      }
    }
    return { ok: true, edges: edges, weighted: weighted };
  }

  function computeDegrees(n, edges) {
    var deg = new Array(n).fill(0);
    edges.forEach(function (e) {
      deg[e[0] - 1]++;
      deg[e[1] - 1]++;
    });
    return deg;
  }

  /* Weighted degree: sum of incident edge weights (distinct from plain
     degree, which is the count of incident edges). */
  function computeWeightedDegrees(n, edges) {
    var wd = new Array(n).fill(0);
    edges.forEach(function (e) {
      wd[e[0] - 1] += e[2];
      wd[e[1] - 1] += e[2];
    });
    return wd;
  }

  function countComponents(n, edges) {
    var parent = [];
    for (var i = 0; i <= n; i++) parent[i] = i;
    function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
    edges.forEach(function (e) {
      var ra = find(e[0]), rb = find(e[1]);
      if (ra !== rb) parent[ra] = rb;
    });
    var roots = {};
    for (var v = 1; v <= n; v++) roots[find(v)] = true;
    return Object.keys(roots).length;
  }

  var SUBSCRIPT_DIGITS = { '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉' };
  function toSubscript(num) {
    return String(num).split('').map(function (ch) { return SUBSCRIPT_DIGITS[ch] || ch; }).join('');
  }

  /* Only classifies a graph shape when the check is unambiguous; otherwise
     returns null so the UI falls back to a plain vertex/edge count. Uses
     standard mathematical notation (Kn, Pn, Cn, K1,n-1) with Unicode
     subscripts for the vertex count. */
  function classifyGraph(n, edges, degrees) {
    var m = edges.length;
    if (n === 1) return 'Single Vertex';
    if (m === 0) return 'Empty Graph (no edges)';

    var components = countComponents(n, edges);
    if (components > 1) return 'Disconnected Graph (' + components + ' components)';

    var maxEdges = (n * (n - 1)) / 2;
    if (m === maxEdges) return 'Complete Graph K' + toSubscript(n);

    var deg1Count = 0, deg2Count = 0, otherDeg = false, maxDeg = 0;
    degrees.forEach(function (d) {
      if (d === 1) deg1Count++;
      else if (d === 2) deg2Count++;
      else otherDeg = true;
      if (d > maxDeg) maxDeg = d;
    });

    if (!otherDeg && m === n - 1 && deg1Count === 2 && deg2Count === n - 2) {
      return 'Path Graph P' + toSubscript(n);
    }
    if (!otherDeg && m === n && deg2Count === n && n >= 3) {
      return 'Cycle Graph C' + toSubscript(n);
    }
    if (m === n - 1 && maxDeg === n - 1 && deg1Count === n - 1) {
      return 'Star Graph K' + toSubscript(1) + ',' + toSubscript(n - 1);
    }
    return null;
  }

  /* ---- Recognizable-shape layout ----
     Rather than always arranging vertices on one plain circle, the layout
     recognizes common graph shapes (path, cycle, star, complete, and each
     connected component of a disconnected graph) and draws each with the
     arrangement a textbook would use. This never changes the edge set --
     it only decides where each vertex is drawn. */

  function componentsList(n, edges) {
    var parent = [];
    for (var i = 0; i <= n; i++) parent[i] = i;
    function find(x) { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }
    edges.forEach(function (e) {
      var ra = find(e[0]), rb = find(e[1]);
      if (ra !== rb) parent[ra] = rb;
    });
    var groups = {};
    for (var v = 1; v <= n; v++) {
      var r = find(v);
      (groups[r] = groups[r] || []).push(v);
    }
    var list = Object.keys(groups).map(function (k) { return groups[k]; });
    list.forEach(function (c) { c.sort(function (a, b) { return a - b; }); });
    list.sort(function (a, b) { return a[0] - b[0]; });
    return list;
  }

  /* Walks a degree-<=2 chain starting at `start`, always stepping to the
     lowest-numbered unvisited neighbor -- used to recover the actual
     path/cycle order from the edge list, which need not match vertex ID
     order in the matrix. */
  function traceChain(start, adjacency, vertexCount) {
    var order = [start];
    var visited = {};
    visited[start] = true;
    var current = start;
    while (order.length < vertexCount) {
      var neighbors = adjacency[current].slice().sort(function (a, b) { return a - b; });
      var next = neighbors.filter(function (v) { return !visited[v]; })[0];
      if (next === undefined) break;
      order.push(next);
      visited[next] = true;
      current = next;
    }
    return order;
  }

  /* Identifies a single connected component's shape from its own vertices
     and edges only (never from the whole graph), so each component of a
     disconnected graph is classified and drawn independently. */
  function detectComponentShape(vertexIds, subEdges) {
    var n = vertexIds.length;
    var m = subEdges.length;
    if (n === 1) return { type: 'single' };

    var adjacency = {};
    vertexIds.forEach(function (v) { adjacency[v] = []; });
    subEdges.forEach(function (e) {
      adjacency[e[0]].push(e[1]);
      adjacency[e[1]].push(e[0]);
    });

    var maxEdges = (n * (n - 1)) / 2;
    if (m === maxEdges) return { type: 'complete', ids: vertexIds };

    var deg1 = [], deg2Count = 0, maxDeg = 0, hub = null;
    vertexIds.forEach(function (v) {
      var d = adjacency[v].length;
      if (d === 1) deg1.push(v);
      if (d === 2) deg2Count++;
      if (d > maxDeg) { maxDeg = d; hub = v; }
    });

    if (m === n - 1 && deg1.length === 2 && deg2Count === n - 2) {
      return { type: 'path', order: traceChain(deg1[0], adjacency, n) };
    }
    if (m === n && deg2Count === n && n >= 3) {
      var startId = vertexIds[0];
      return { type: 'cycle', order: traceChain(startId, adjacency, n) };
    }
    if (m === n - 1 && maxDeg === n - 1 && deg1.length === n - 1) {
      return { type: 'star', hub: hub, leaves: vertexIds.filter(function (v) { return v !== hub; }) };
    }
    return { type: 'generic', ids: vertexIds };
  }

  /* Places one component's vertices into a local (cx, cy, radius) budget
     according to its detected shape, writing results into `pointsOut`
     keyed by vertex ID. */
  function layoutComponent(shape, vertexIds, cx, cy, radius, pointsOut) {
    var n = vertexIds.length;
    if (n === 1 || shape.type === 'single') {
      pointsOut[vertexIds[0]] = { x: cx, y: cy };
      return;
    }

    if (shape.type === 'path') {
      var order = shape.order;
      var span = radius;
      order.forEach(function (v, i) {
        var t = order.length === 1 ? 0.5 : i / (order.length - 1);
        pointsOut[v] = { x: cx - span + t * 2 * span, y: cy };
      });
      return;
    }

    if (shape.type === 'cycle' && n === 4) {
      var half = radius * 0.78;
      var corners = [
        { x: cx - half, y: cy - half },
        { x: cx + half, y: cy - half },
        { x: cx + half, y: cy + half },
        { x: cx - half, y: cy + half }
      ];
      shape.order.forEach(function (v, i) { pointsOut[v] = corners[i]; });
      return;
    }

    if (shape.type === 'cycle') {
      shape.order.forEach(function (v, i) {
        var angle = -Math.PI / 2 + (2 * Math.PI * i) / shape.order.length;
        pointsOut[v] = { x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle) };
      });
      return;
    }

    if (shape.type === 'star') {
      pointsOut[shape.hub] = { x: cx, y: cy };
      shape.leaves.forEach(function (v, i) {
        var angle = -Math.PI / 2 + (2 * Math.PI * i) / shape.leaves.length;
        pointsOut[v] = { x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle) };
      });
      return;
    }

    if (shape.type === 'complete' && n === 4) {
      var sin120 = Math.sin(2 * Math.PI / 3);
      var positions = [
        { x: cx, y: cy - radius },
        { x: cx - radius * sin120, y: cy + radius * 0.5 },
        { x: cx + radius * sin120, y: cy + radius * 0.5 },
        { x: cx, y: cy }
      ];
      vertexIds.forEach(function (v, i) { pointsOut[v] = positions[i]; });
      return;
    }

    /* Generic fallback (larger complete graphs, or any shape that doesn't
       match a recognized pattern): a symmetric circle in vertex-ID order. */
    vertexIds.forEach(function (v, i) {
      var angle = -Math.PI / 2 + (2 * Math.PI * i) / n;
      pointsOut[v] = { x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle) };
    });
  }

  /* Top-level layout: a single connected graph is drawn centered in the
     whole canvas using its detected shape. A disconnected graph instead
     gets one grid cell per component, so components are clearly separated
     rather than interleaved on one shared circle. */
  function computeLayoutPoints(n, edges, size, nodeRadius) {
    var center = size / 2;
    var components = componentsList(n, edges);
    var pts = {};

    if (components.length <= 1) {
      var subEdgesAll = edges;
      var shape = detectComponentShape(components[0] || [1], subEdgesAll);
      var radius = n === 1 ? 0 : center - nodeRadius - 34;
      layoutComponent(shape, components[0] || [1], center, center, radius, pts);
      return pts;
    }

    var k = components.length;
    var cols = Math.ceil(Math.sqrt(k));
    var rows = Math.ceil(k / cols);
    var cellW = size / cols;
    var cellH = size / rows;

    components.forEach(function (comp, idx) {
      var col = idx % cols;
      var row = Math.floor(idx / cols);
      var cellCx = cellW * col + cellW / 2;
      var cellCy = cellH * row + cellH / 2;
      var cellRadius = Math.max(12, Math.min(cellW, cellH) / 2 - nodeRadius - 20);
      var subEdges = edges.filter(function (e) { return comp.indexOf(e[0]) !== -1; });
      var shape = detectComponentShape(comp, subEdges);
      layoutComponent(shape, comp, cellCx, cellCy, cellRadius, pts);
    });

    return pts;
  }

  /* ---- Weighted layout: edge length communicates edge weight ----
     Unweighted graphs keep the shape-aware layout above untouched. When a
     graph is weighted, edge length instead needs to visually track weight
     (smaller weight -> shorter edge), which the fixed shape templates
     above cannot express since they space every edge equally. This uses a
     lightweight, damped, bounded force-directed relaxation: each edge acts
     like a spring whose rest length is the weight normalized (min-max)
     into a fixed pixel range, and vertices repel each other so they don't
     collide. Normalization keeps a handful of extreme weights from
     blowing up the drawing, and the simulation always starts from the
     same shape-aware layout used for unweighted graphs (rather than a
     random scatter), so a weighted graph still reads as recognizably the
     same topology, just stretched/compressed by weight. */
  function desiredEdgeLength(weight, minW, maxW, minLen, maxLen) {
    if (maxW - minW <= 1e-9) return (minLen + maxLen) / 2;
    var t = (weight - minW) / (maxW - minW);
    return minLen + t * (maxLen - minLen);
  }

  function forceDirectedLayout(vertexIds, subEdges, centerX, centerY, radius, nodeRadius, minW, maxW, seedPositions, pointsOut) {
    var k = vertexIds.length;
    if (k === 0) return;
    if (k === 1) {
      var only = seedPositions && seedPositions[vertexIds[0]];
      pointsOut[vertexIds[0]] = only ? { x: only.x, y: only.y } : { x: centerX, y: centerY };
      return;
    }

    var minLen = Math.max(nodeRadius * 3.4, radius * 0.42);
    var maxLen = Math.max(minLen + 1, radius * 1.5);

    var pos = {};
    var vel = {};
    vertexIds.forEach(function (v) {
      var seed = seedPositions && seedPositions[v];
      pos[v] = seed ? { x: seed.x, y: seed.y } : { x: centerX, y: centerY };
      vel[v] = { x: 0, y: 0 };
    });

    var edgeList = subEdges.map(function (e) {
      return { a: e[0], b: e[1], len: desiredEdgeLength(e[2], minW, maxW, minLen, maxLen) };
    });

    var iterations = k > 60 ? 180 : 380;
    var minSeparation = nodeRadius * 2.6;
    var repulsionStrength = radius * radius * 0.9;
    var springStrength = 0.06;
    var centerPull = 0.002;
    var damping = 0.85;
    var maxBound = radius + nodeRadius * 0.5;

    for (var iter = 0; iter < iterations; iter++) {
      var cooling = 1 - iter / iterations;
      var force = {};
      vertexIds.forEach(function (v) { force[v] = { x: 0, y: 0 }; });

      for (var i = 0; i < k; i++) {
        for (var j = i + 1; j < k; j++) {
          var vi = vertexIds[i], vj = vertexIds[j];
          var dx = pos[vi].x - pos[vj].x;
          var dy = pos[vi].y - pos[vj].y;
          var distSq = dx * dx + dy * dy;
          var dist = Math.sqrt(distSq) || 0.01;
          var rep = repulsionStrength / distSq;
          var fx = (dx / dist) * rep;
          var fy = (dy / dist) * rep;
          force[vi].x += fx; force[vi].y += fy;
          force[vj].x -= fx; force[vj].y -= fy;
        }
      }

      edgeList.forEach(function (e) {
        var pa = pos[e.a], pb = pos[e.b];
        var dx = pb.x - pa.x, dy = pb.y - pa.y;
        var dist = Math.sqrt(dx * dx + dy * dy) || 0.01;
        var diff = dist - e.len;
        var fx = (dx / dist) * diff * springStrength;
        var fy = (dy / dist) * diff * springStrength;
        force[e.a].x += fx; force[e.a].y += fy;
        force[e.b].x -= fx; force[e.b].y -= fy;
      });

      vertexIds.forEach(function (v) {
        force[v].x += (centerX - pos[v].x) * centerPull;
        force[v].y += (centerY - pos[v].y) * centerPull;
      });

      vertexIds.forEach(function (v) {
        vel[v].x = (vel[v].x + force[v].x * cooling) * damping;
        vel[v].y = (vel[v].y + force[v].y * cooling) * damping;
        pos[v].x += vel[v].x;
        pos[v].y += vel[v].y;

        var ddx = pos[v].x - centerX, ddy = pos[v].y - centerY;
        var d = Math.sqrt(ddx * ddx + ddy * ddy);
        if (d > maxBound) {
          var scale = maxBound / d;
          pos[v].x = centerX + ddx * scale;
          pos[v].y = centerY + ddy * scale;
          vel[v].x *= 0.3; vel[v].y *= 0.3;
        }
      });
    }

    for (var pass = 0; pass < 6; pass++) {
      for (var i2 = 0; i2 < k; i2++) {
        for (var j2 = i2 + 1; j2 < k; j2++) {
          var va = vertexIds[i2], vb = vertexIds[j2];
          var ddx2 = pos[vb].x - pos[va].x;
          var ddy2 = pos[vb].y - pos[va].y;
          var d2 = Math.sqrt(ddx2 * ddx2 + ddy2 * ddy2) || 0.01;
          if (d2 < minSeparation) {
            var push = (minSeparation - d2) / 2;
            var nx = ddx2 / d2, ny = ddy2 / d2;
            pos[va].x -= nx * push; pos[va].y -= ny * push;
            pos[vb].x += nx * push; pos[vb].y += ny * push;
          }
        }
      }
    }

    vertexIds.forEach(function (v) { pointsOut[v] = pos[v]; });
  }

  /* Path and cycle (and star) shapes have a clean 1-D / radial
     parametrization, so their weighted layout is computed directly rather
     than through the general force simulation -- this stays legible at
     any vertex count (no iterative physics to destabilize a large ring),
     where letting every vertex repel every other vertex tends to tear a
     large cycle apart into a tangle instead of a recognizable ring. */
  function edgeWeightLookup(subEdges) {
    var map = {};
    subEdges.forEach(function (e) {
      map[e[0] + '-' + e[1]] = e[2];
      map[e[1] + '-' + e[0]] = e[2];
    });
    return map;
  }

  function layoutWeightedPath(order, cx, cy, radius, weightMap, minW, maxW, pointsOut) {
    var n = order.length;
    var segLens = [];
    for (var i = 0; i < n - 1; i++) {
      segLens.push(desiredEdgeLength(weightMap[order[i] + '-' + order[i + 1]], minW, maxW, 1, 3));
    }
    var total = segLens.reduce(function (a, b) { return a + b; }, 0);
    var span = radius * 2;
    var scale = total > 0 ? span / total : 0;
    var x = cx - radius, y = cy;
    pointsOut[order[0]] = { x: x, y: y };
    for (i = 0; i < n - 1; i++) {
      x += segLens[i] * scale;
      pointsOut[order[i + 1]] = { x: x, y: y };
    }
  }

  function layoutWeightedCycle(order, cx, cy, radius, weightMap, minW, maxW, pointsOut) {
    var n = order.length;
    var shares = [];
    for (var i = 0; i < n; i++) {
      var a = order[i], b = order[(i + 1) % n];
      shares.push(desiredEdgeLength(weightMap[a + '-' + b], minW, maxW, 1, 3));
    }
    var total = shares.reduce(function (a, b) { return a + b; }, 0);
    var angle = -Math.PI / 2;
    for (i = 0; i < n; i++) {
      pointsOut[order[i]] = { x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle) };
      angle += (total > 0 ? shares[i] / total : 1 / n) * 2 * Math.PI;
    }
  }

  function layoutWeightedStar(hub, leaves, cx, cy, radius, weightMap, minW, maxW, pointsOut) {
    pointsOut[hub] = { x: cx, y: cy };
    var minLen = Math.max(radius * 0.35, 18);
    var maxLen = radius * 0.92;
    leaves.forEach(function (v, i) {
      var w = weightMap[hub + '-' + v];
      var len = desiredEdgeLength(w, minW, maxW, minLen, maxLen);
      var angle = -Math.PI / 2 + (2 * Math.PI * i) / leaves.length;
      pointsOut[v] = { x: cx + len * Math.cos(angle), y: cy + len * Math.sin(angle) };
    });
  }

  /* Above this size, a fully connected or unrecognized ("generic") shape
     is too dense for a force simulation to settle into anything legible
     -- it reads as a hairball whether weighted or not -- so it falls back
     to the plain equal-length layout instead of spending time distorting
     something that will look tangled either way. Weight labels are
     unaffected and still show. */
  var FORCE_DIRECTED_MAX_K = 14;

  /* Weighted counterpart to computeLayoutPoints: same component/grid
     structure (so disconnected weighted graphs still separate into their
     own cells). Path, cycle, and star components get an exact analytic
     placement so edge length reflects weight while the shape stays
     clean at any size; small complete/generic components are seeded from
     their shape-aware layout and relaxed with the bounded force
     simulation; larger complete/generic components keep the plain
     equal-length layout, since weight-driven stretching would only add
     visual noise once a graph is already dense. If every edge has the
     same weight there is nothing to visualize, so the whole graph defers
     to the plain shape-aware layout, per the requirement that a uniform
     weighted graph reads as cleanly as an unweighted one. */
  function computeWeightedLayoutPoints(n, edges, size, nodeRadius) {
    var weights = edges.map(function (e) { return e[2]; });
    var minW = weights.length ? Math.min.apply(null, weights) : 0;
    var maxW = weights.length ? Math.max.apply(null, weights) : 0;
    if (maxW - minW <= 1e-9) return computeLayoutPoints(n, edges, size, nodeRadius);

    var center = size / 2;
    var components = componentsList(n, edges);
    var pts = {};

    function layoutOneComponent(compVertices, subEdges, cx, cy, radius) {
      if (compVertices.length <= 1) {
        pts[compVertices[0]] = { x: cx, y: cy };
        return;
      }
      var shape = detectComponentShape(compVertices, subEdges);
      var weightMap = edgeWeightLookup(subEdges);

      if (shape.type === 'path') {
        layoutWeightedPath(shape.order, cx, cy, radius, weightMap, minW, maxW, pts);
        return;
      }
      if (shape.type === 'cycle') {
        layoutWeightedCycle(shape.order, cx, cy, radius, weightMap, minW, maxW, pts);
        return;
      }
      if (shape.type === 'star') {
        layoutWeightedStar(shape.hub, shape.leaves, cx, cy, radius, weightMap, minW, maxW, pts);
        return;
      }

      if (compVertices.length <= FORCE_DIRECTED_MAX_K) {
        var seed = {};
        layoutComponent(shape, compVertices, cx, cy, radius, seed);
        forceDirectedLayout(compVertices, subEdges, cx, cy, radius, nodeRadius, minW, maxW, seed, pts);
      } else {
        layoutComponent(shape, compVertices, cx, cy, radius, pts);
      }
    }

    if (components.length <= 1) {
      var radius = n === 1 ? 0 : center - nodeRadius - 34;
      layoutOneComponent(components[0] || [1], edges, center, center, radius);
      return pts;
    }

    var k = components.length;
    var cols = Math.ceil(Math.sqrt(k));
    var rows = Math.ceil(k / cols);
    var cellW = size / cols;
    var cellH = size / rows;

    components.forEach(function (comp, idx) {
      var col = idx % cols;
      var row = Math.floor(idx / cols);
      var cellCx = cellW * col + cellW / 2;
      var cellCy = cellH * row + cellH / 2;
      var cellRadius = Math.max(12, Math.min(cellW, cellH) / 2 - nodeRadius - 20);
      var subEdges = edges.filter(function (e) { return comp.indexOf(e[0]) !== -1; });
      layoutOneComponent(comp, subEdges, cellCx, cellCy, cellRadius);
    });

    return pts;
  }

  function svgEl(tag, attrs) {
    var el = document.createElementNS(SVG_NS, tag);
    if (attrs) {
      for (var key in attrs) {
        if (Object.prototype.hasOwnProperty.call(attrs, key)) el.setAttribute(key, attrs[key]);
      }
    }
    return el;
  }

  /* Node size / label visibility scale down as vertex count grows so that
     graphs up to 100 vertices stay legible instead of a solid blob. */
  function graphVisualParams(n) {
    if (n <= 12) return { r: 16, showLabels: true, strokeWidth: 1.6, fontSize: 13 };
    if (n <= 25) return { r: 12, showLabels: true, strokeWidth: 1.3, fontSize: 10 };
    if (n <= 45) return { r: 8, showLabels: true, strokeWidth: 1.0, fontSize: 7 };
    if (n <= 70) return { r: 6, showLabels: false, strokeWidth: 0.9, fontSize: 0 };
    return { r: 4.5, showLabels: false, strokeWidth: 0.7, fontSize: 0 };
  }

  /* Dense graphs (many overlapping straight edges) otherwise render as a
     solid filled disc. Fading edges by count keeps the ring of vertices
     and the overall density legible instead of a single opaque blob. */
  function edgeOpacityForCount(edgeCount) {
    if (edgeCount > 1500) return 0.08;
    if (edgeCount > 600) return 0.18;
    if (edgeCount > 200) return 0.35;
    if (edgeCount > 60) return 0.55;
    return 0.9;
  }

  function buildGraphSVG(n, edges, degrees, weighted, weightedDegrees) {
    var size = 520;
    var center = size / 2;
    var params = graphVisualParams(n);
    var pts = weighted
      ? computeWeightedLayoutPoints(n, edges, size, params.r)
      : computeLayoutPoints(n, edges, size, params.r);
    var baseOpacity = edgeOpacityForCount(edges.length);
    var showWeightLabels = weighted && params.showLabels;

    /* ---- Adaptive viewport ----
       The layout above always positions vertices within the same size x
       size logical space as before (unchanged math). What changes here is
       only how much of that space is actually shown: a small/simple graph
       (e.g. a 4-vertex path) otherwise sits as a tiny sliver in the middle
       of a mostly-empty square. Instead, fit a tight viewBox (and a
       matching container aspect ratio) around the graph's own content --
       vertices plus, when shown, weight-label boxes -- so the drawing
       fills the available space. The aspect ratio is clamped so a very
       short/degenerate graph (e.g. a single edge) doesn't become an
       unreasonably thin sliver of its own. */
    var bbMinX = Infinity, bbMinY = Infinity, bbMaxX = -Infinity, bbMaxY = -Infinity;
    function growBBox(x0, y0, x1, y1) {
      if (x0 < bbMinX) bbMinX = x0;
      if (y0 < bbMinY) bbMinY = y0;
      if (x1 > bbMaxX) bbMaxX = x1;
      if (y1 > bbMaxY) bbMaxY = y1;
    }
    for (var bv = 1; bv <= n; bv++) {
      var bp = pts[bv];
      growBBox(bp.x - params.r, bp.y - params.r, bp.x + params.r, bp.y + params.r);
    }
    if (showWeightLabels) {
      edges.forEach(function (e) {
        var ba = pts[e[0]], bb = pts[e[1]];
        var bMidX = (ba.x + bb.x) / 2, bMidY = (ba.y + bb.y) / 2;
        var bDx = bb.x - ba.x, bDy = bb.y - ba.y;
        var bLen = Math.sqrt(bDx * bDx + bDy * bDy) || 1;
        var bPerpX = -bDy / bLen, bPerpY = bDx / bLen;
        var bOffset = Math.max(11, params.fontSize * 0.9);
        var bLx = bMidX + bPerpX * bOffset, bLy = bMidY + bPerpY * bOffset;
        var bWStr = formatReal(e[2]);
        var bFSize = Math.max(params.fontSize - 1, 9);
        var bBoxW = Math.max(bWStr.length * bFSize * 0.62 + 8, bFSize + 8);
        var bBoxH = bFSize + 6;
        growBBox(bLx - bBoxW / 2, bLy - bBoxH / 2, bLx + bBoxW / 2, bLy + bBoxH / 2);
      });
    }
    if (!isFinite(bbMinX)) { bbMinX = 0; bbMinY = 0; bbMaxX = size; bbMaxY = size; }

    var VIEW_PAD = 26;
    var MIN_VIEW_DIM = 170;
    var MAX_VIEW_RATIO = 2.6;
    var MIN_VIEW_RATIO = 1 / MAX_VIEW_RATIO;
    var viewW = Math.max(bbMaxX - bbMinX + VIEW_PAD * 2, MIN_VIEW_DIM);
    var viewH = Math.max(bbMaxY - bbMinY + VIEW_PAD * 2, MIN_VIEW_DIM);
    var viewRatio = viewW / viewH;
    if (viewRatio > MAX_VIEW_RATIO) viewH = viewW / MAX_VIEW_RATIO;
    else if (viewRatio < MIN_VIEW_RATIO) viewW = viewH * MIN_VIEW_RATIO;
    var bbCx = (bbMinX + bbMaxX) / 2, bbCy = (bbMinY + bbMaxY) / 2;
    var viewX = bbCx - viewW / 2, viewY = bbCy - viewH / 2;

    graphContainer.style.aspectRatio = (viewW / viewH).toFixed(4) + ' / 1';

    var styles = getComputedStyle(document.documentElement);
    function cssVar(name, fallback) {
      var v = styles.getPropertyValue(name);
      return v ? v.trim() : fallback;
    }
    var colorEdge = cssVar('--color-text-muted', '#5c6a62');
    var colorNodeFill = cssVar('--color-surface', '#ffffff');
    var colorNodeStroke = cssVar('--color-green', '#078055');
    var colorLabel = cssVar('--color-text', '#1b2420');
    var colorHighlight = cssVar('--color-pink', '#b84472');
    var colorHighlightSoft = cssVar('--color-pink-soft', '#fbebf1');

    var svg = svgEl('svg', {
      viewBox: viewX.toFixed(2) + ' ' + viewY.toFixed(2) + ' ' + viewW.toFixed(2) + ' ' + viewH.toFixed(2),
      width: '100%',
      height: '100%',
      preserveAspectRatio: 'xMidYMid meet',
      class: 'graph-svg',
      focusable: 'false'
    });

    var titleEl = svgEl('title');
    titleEl.textContent = 'Graph with ' + n + ' vert' + (n === 1 ? 'ex' : 'ices') + ' and ' + edges.length +
      ' edge' + (edges.length === 1 ? '' : 's') + ', reconstructed from the Laplacian matrix.' +
      (weighted ? ' Edge weights are shown.' : '');
    svg.appendChild(titleEl);
    svg.appendChild(svgEl('rect', {
      x: viewX.toFixed(2), y: viewY.toFixed(2), width: viewW.toFixed(2), height: viewH.toFixed(2),
      fill: colorNodeFill, rx: 14
    }));

    var edgeLayer = svgEl('g', { class: 'graph-edges' });
    var nodeLayer = svgEl('g', { class: 'graph-nodes' });
    var adjacency = {};
    var edgeWeightGroups = {};
    var edgeLineByKey = {};
    for (var v = 1; v <= n; v++) adjacency[v] = [];

    edges.forEach(function (e) {
      var a = pts[e[0]], b = pts[e[1]];
      var edgeKey = e[0] + '-' + e[1];
      var line = svgEl('line', {
        x1: a.x.toFixed(2), y1: a.y.toFixed(2),
        x2: b.x.toFixed(2), y2: b.y.toFixed(2),
        stroke: colorEdge,
        'stroke-width': params.strokeWidth,
        'stroke-opacity': baseOpacity,
        'stroke-linecap': 'round',
        class: 'graph-edge',
        'data-edge': edgeKey
      });
      var edgeTitle = svgEl('title');
      edgeTitle.textContent = weighted
        ? ('Edge ' + e[0] + '-' + e[1] + '\nWeight: ' + formatReal(e[2]))
        : ('Edge ' + e[0] + '-' + e[1]);
      line.appendChild(edgeTitle);
      edgeLayer.appendChild(line);
      adjacency[e[0]].push(line);
      adjacency[e[1]].push(line);
      edgeLineByKey[edgeKey] = line;

      if (showWeightLabels) {
        var midX = (a.x + b.x) / 2;
        var midY = (a.y + b.y) / 2;
        var edgeDx = b.x - a.x, edgeDy = b.y - a.y;
        var edgeLen = Math.sqrt(edgeDx * edgeDx + edgeDy * edgeDy) || 1;
        var perpX = -edgeDy / edgeLen, perpY = edgeDx / edgeLen;
        var labelOffset = Math.max(11, params.fontSize * 0.9);
        var mx = midX + perpX * labelOffset;
        var my = midY + perpY * labelOffset;
        var wStr = formatReal(e[2]);
        var fSize = Math.max(params.fontSize - 1, 9);
        var boxW = Math.max(wStr.length * fSize * 0.62 + 8, fSize + 8);
        var boxH = fSize + 6;
        var weightGroup = svgEl('g', { class: 'graph-edge-weight', 'data-edge': edgeKey });
        var weightRect = svgEl('rect', {
          x: (mx - boxW / 2).toFixed(2), y: (my - boxH / 2).toFixed(2),
          width: boxW.toFixed(2), height: boxH.toFixed(2),
          rx: 4, fill: colorNodeFill, 'fill-opacity': 0.92,
          stroke: colorEdge, 'stroke-width': 1
        });
        var weightText = svgEl('text', {
          x: mx.toFixed(2), y: my.toFixed(2),
          'text-anchor': 'middle', 'dominant-baseline': 'central',
          'font-size': fSize, 'font-weight': 600,
          fill: colorLabel, class: 'graph-edge-weight-text'
        });
        weightText.textContent = wStr;
        var weightTitle = svgEl('title');
        weightTitle.textContent = 'Edge ' + e[0] + '-' + e[1] + '\nWeight: ' + wStr;
        weightGroup.appendChild(weightTitle);
        weightGroup.appendChild(weightRect);
        weightGroup.appendChild(weightText);
        edgeLayer.appendChild(weightGroup);
        edgeWeightGroups[edgeKey] = { rect: weightRect, text: weightText };
      }
    });

    var vertexGroups = {};
    for (var i = 1; i <= n; i++) {
      var p = pts[i];
      var g = svgEl('g', {
        class: 'graph-vertex',
        'data-vertex': i,
        tabindex: '0',
        role: 'button',
        'aria-label': 'Vertex ' + i + ', degree ' + (degrees[i - 1] || 0)
      });

      var circle = svgEl('circle', {
        cx: p.x.toFixed(2), cy: p.y.toFixed(2), r: params.r,
        fill: colorNodeFill, stroke: colorNodeStroke, 'stroke-width': 2,
        class: 'graph-vertex-circle'
      });
      var vTitle = svgEl('title');
      vTitle.textContent = 'Vertex ' + i + ' (degree ' + (degrees[i - 1] || 0) + ')';
      circle.appendChild(vTitle);
      g.appendChild(circle);

      if (params.showLabels) {
        var text = svgEl('text', {
          x: p.x.toFixed(2), y: p.y.toFixed(2),
          'text-anchor': 'middle', 'dominant-baseline': 'central',
          'font-size': params.fontSize,
          fill: colorLabel, class: 'graph-vertex-label'
        });
        text.textContent = String(i);
        g.appendChild(text);
      }

      nodeLayer.appendChild(g);
      vertexGroups[i] = g;
    }

    svg.appendChild(edgeLayer);
    svg.appendChild(nodeLayer);

    function setEdgeHighlight(line, on) {
      line.setAttribute('stroke', on ? colorHighlight : colorEdge);
      line.setAttribute('stroke-width', on ? (params.strokeWidth + 1.4) : params.strokeWidth);
      line.setAttribute('stroke-opacity', on ? 1 : baseOpacity);
      var key = line.getAttribute('data-edge');
      var wg = key && edgeWeightGroups[key];
      if (wg) {
        wg.rect.setAttribute('stroke', on ? colorHighlight : colorEdge);
        wg.rect.setAttribute('stroke-width', on ? 1.6 : 1);
        wg.text.setAttribute('fill', on ? colorHighlight : colorLabel);
      }
    }
    function setVertexHighlight(vid, on) {
      var g = vertexGroups[vid];
      if (!g) return;
      var circle = g.querySelector('circle');
      circle.setAttribute('fill', on ? colorHighlightSoft : colorNodeFill);
      circle.setAttribute('stroke', on ? colorHighlight : colorNodeStroke);
      circle.setAttribute('stroke-width', on ? 3 : 2);
      (adjacency[vid] || []).forEach(function (line) { setEdgeHighlight(line, on); });
    }
    function closestClass(target, className) {
      return (target && target.closest) ? target.closest('.' + className) : null;
    }

    svg.addEventListener('mouseover', function (evt) {
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) { setVertexHighlight(parseInt(vg.getAttribute('data-vertex'), 10), true); return; }
      var el = closestClass(evt.target, 'graph-edge');
      if (el) { setEdgeHighlight(el, true); return; }
      var wl = closestClass(evt.target, 'graph-edge-weight');
      if (wl) {
        var wlLine = edgeLineByKey[wl.getAttribute('data-edge')];
        if (wlLine) setEdgeHighlight(wlLine, true);
      }
    });
    svg.addEventListener('mouseout', function (evt) {
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) { setVertexHighlight(parseInt(vg.getAttribute('data-vertex'), 10), false); return; }
      var el = closestClass(evt.target, 'graph-edge');
      if (el) { setEdgeHighlight(el, false); return; }
      var wl = closestClass(evt.target, 'graph-edge-weight');
      if (wl) {
        var wlLine = edgeLineByKey[wl.getAttribute('data-edge')];
        if (wlLine) setEdgeHighlight(wlLine, false);
      }
    });
    svg.addEventListener('focusin', function (evt) {
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) setVertexHighlight(parseInt(vg.getAttribute('data-vertex'), 10), true);
    });
    svg.addEventListener('focusout', function (evt) {
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) setVertexHighlight(parseInt(vg.getAttribute('data-vertex'), 10), false);
    });

    function activateVertex(vid) {
      if (weighted) {
        graphVertexInfo.textContent = 'Vertex ' + vid + ' — Degree: ' + (degrees[vid - 1] || 0) +
          ', Weighted Degree: ' + formatReal((weightedDegrees && weightedDegrees[vid - 1]) || 0) + '.';
      } else {
        graphVertexInfo.textContent = 'Vertex ' + vid + ' — degree ' + (degrees[vid - 1] || 0) + '.';
      }
    }
    svg.addEventListener('click', function (evt) {
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) activateVertex(parseInt(vg.getAttribute('data-vertex'), 10));
    });
    svg.addEventListener('keydown', function (evt) {
      if (evt.key !== 'Enter' && evt.key !== ' ') return;
      var vg = closestClass(evt.target, 'graph-vertex');
      if (vg) { evt.preventDefault(); activateVertex(parseInt(vg.getAttribute('data-vertex'), 10)); }
    });

    return svg;
  }

  function renderGraphSection(matrix, n) {
    var analysis = laplacianAnalysis(matrix, n);
    if (!analysis.ok) {
      graphWrap.hidden = true;
      graphUnavailable.hidden = false;
      graphUnavailable.textContent = 'Graph visualization is unavailable: ' + analysis.reason;
      graphWeightedBadge.hidden = true;
      lastGraph = null;
      return;
    }

    graphUnavailable.hidden = true;
    graphWrap.hidden = false;
    graphVertexInfo.textContent = '';

    var edges = analysis.edges;
    var weighted = analysis.weighted;
    var degrees = computeDegrees(n, edges);
    var weightedDegrees = weighted ? computeWeightedDegrees(n, edges) : null;
    var svg = buildGraphSVG(n, edges, degrees, weighted, weightedDegrees);

    graphContainer.innerHTML = '';
    graphContainer.appendChild(svg);
    graphContainer.setAttribute('aria-label', weighted
      ? 'Weighted graph reconstructed from the Laplacian matrix, with edge weights shown'
      : 'Graph reconstructed from the Laplacian matrix');

    graphWeightedBadge.hidden = !weighted;
    graphHint.textContent = weighted
      ? 'Hover a vertex or edge to highlight it. Click a vertex to view its degree and weighted degree.'
      : 'Hover a vertex or edge to highlight it. Click a vertex to see its degree.';

    var shapeLabel = classifyGraph(n, edges, degrees);
    var typeLabel = weighted ? (shapeLabel ? ('Weighted ' + shapeLabel) : 'Weighted Graph') : shapeLabel;
    var infoHtml = '<strong>Vertices:</strong> ' + n + ' &nbsp;•&nbsp; <strong>Edges:</strong> ' + edges.length;
    if (typeLabel) infoHtml += '<br><strong>Graph Type:</strong> ' + escapeHtml(typeLabel);
    graphInfo.innerHTML = infoHtml;

    var params = graphVisualParams(n);
    var densityMessages = [];
    if (n > 60 || edges.length > 400) {
      densityMessages.push('This graph is large and dense. Vertex numbers are shown on hover or focus instead of as permanent labels, and individual edges may be hard to distinguish visually.');
    }
    if (weighted && !params.showLabels) {
      densityMessages.push('Edge weights are shown on hover instead of as permanent labels due to the graph size.');
    }
    if (densityMessages.length) {
      graphDensityNote.hidden = false;
      graphDensityNote.textContent = densityMessages.join(' ');
    } else {
      graphDensityNote.hidden = true;
    }

    lastGraph = { n: n, edges: edges, type: typeLabel, weighted: weighted };
  }

  function buildGraphCopyText() {
    if (!lastGraph) return '';
    var lines = [];
    lines.push('Vertices: ' + lastGraph.n);
    lines.push('Edges: ' + lastGraph.edges.length);
    if (lastGraph.type) lines.push('Graph Type: ' + lastGraph.type);
    lines.push('');
    lines.push('Edges:');
    lastGraph.edges.forEach(function (e) {
      lines.push(lastGraph.weighted ? (e[0] + '-' + e[1] + ': weight ' + formatReal(e[2])) : (e[0] + '-' + e[1]));
    });
    return lines.join('\n') + '\n';
  }

  function copyTextToClipboard(text, btn, failureMessage) {
    if (!text) return;
    var restoreLabel = btn.textContent;
    function onSuccess() {
      btn.textContent = '✓ Copied';
      btn.classList.add('btn--success');
      window.setTimeout(function () {
        btn.textContent = restoreLabel;
        btn.classList.remove('btn--success');
      }, 1500);
    }
    function onFailure() {
      showError(failureMessage);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(onSuccess, onFailure);
    } else {
      try {
        var textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.style.position = 'fixed';
        textarea.style.opacity = '0';
        document.body.appendChild(textarea);
        textarea.select();
        document.execCommand('copy');
        document.body.removeChild(textarea);
        onSuccess();
      } catch (e) {
        onFailure();
      }
    }
  }

  function copyGraphInfo() {
    if (!lastGraph) return;
    copyTextToClipboard(buildGraphCopyText(), copyGraphBtn, 'Could not copy automatically. Please select and copy the graph information manually.');
  }
  copyGraphBtn.addEventListener('click', copyGraphInfo);

  function downloadGraph() {
    var svgSource = graphContainer.querySelector('svg');
    if (!svgSource || !lastGraph) return;
    var clone = svgSource.cloneNode(true);
    clone.setAttribute('xmlns', SVG_NS);
    var serializer = new XMLSerializer();
    var source = '<?xml version="1.0" standalone="no"?>\r\n' + serializer.serializeToString(clone);
    var blob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = 'laplacian-graph.svg';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }
  downloadGraphBtn.addEventListener('click', downloadGraph);

  /* ---------------- Calculate ---------------- */
  function runCalculation() {
    showError('');
    var matrix = readMatrix();
    if (!matrix) return;

    calculateBtn.disabled = true;
    calculatingIndicator.hidden = false;

    window.setTimeout(function () {
      try {
        var groups = computeEigen(matrix, currentN);
        renderResults(matrix, currentN, groups);
        renderGraphSection(matrix, currentN);
      } catch (err) {
        showError('Something went wrong while calculating. Please check your matrix values and try again.');
      } finally {
        calculateBtn.disabled = false;
        calculatingIndicator.hidden = true;
      }
    }, 30);
  }

  calculateBtn.addEventListener('click', runCalculation);

  /* ---------------- Clear / Reset ---------------- */
  function clearMatrixValues() {
    if (!currentN) return;
    for (var i = 0; i < currentN; i++) {
      for (var j = 0; j < currentN; j++) {
        getCell(i, j).value = '';
      }
    }
    showError('');
    clearResults();
    var first = getCell(0, 0);
    if (first) first.focus();
  }

  clearBtn.addEventListener('click', clearMatrixValues);

  /* ---------------- Load example ---------------- */
  var EXAMPLE = [[2, 1], [0, 3]];

  function loadExample() {
    sizeInput.value = '2';
    createMatrix();
    buildMatrixGrid(2, EXAMPLE);
    showError('');
    clearResults();

    var cells = matrixContainer.querySelectorAll('.matrix-cell');
    cells.forEach(function (cell) { cell.classList.add('matrix-cell--filled'); });
    window.setTimeout(function () {
      cells.forEach(function (cell) { cell.classList.remove('matrix-cell--filled'); });
    }, 500);
  }

  exampleBtn.addEventListener('click', loadExample);

  /* ---------------- Copy results ---------------- */
  function buildCopyText() {
    if (!lastResults) return '';
    var lines = [];
    var n = lastResults.n;

    lines.push('Input Matrix:');
    lastResults.matrix.forEach(function (row) {
      lines.push('[' + row.map(formatReal).join(' ') + ']');
    });
    lines.push('');

    lines.push('Eigenvalues:');
    lastResults.groups.forEach(function (g) {
      var line = 'λ = ' + formatComplex(g.value) + ', Algebraic Multiplicity = ' + g.multiplicity;
      if (g.vectors.length < g.multiplicity) {
        line += ', Geometric Multiplicity = ' + g.vectors.length;
      }
      lines.push(line);
    });
    lines.push('');

    lines.push('Eigenvectors:');
    lastResults.groups.forEach(function (g) {
      lines.push('');
      lines.push('λ = ' + formatComplex(g.value));
      g.vectors.forEach(function (v) {
        v.forEach(function (component) {
          lines.push('[' + formatComplex(component) + ']');
        });
        lines.push('');
      });
    });

    return lines.join('\n').replace(/\n{3,}/g, '\n\n').trim() + '\n';
  }

  function copyResults() {
    var text = buildCopyText();
    if (!text) return;

    var restoreLabel = copyBtn.textContent;
    function onSuccess() {
      copyBtn.textContent = '✓ Copied';
      copyBtn.classList.add('btn--success');
      window.setTimeout(function () {
        copyBtn.textContent = restoreLabel;
        copyBtn.classList.remove('btn--success');
      }, 1500);
    }
    function onFailure() {
      showError('Could not copy automatically. Please select and copy the results manually.');
    }

    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(onSuccess, onFailure);
    } else {
      try {
        var textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.style.position = 'fixed';
        textarea.style.opacity = '0';
        document.body.appendChild(textarea);
        textarea.select();
        document.execCommand('copy');
        document.body.removeChild(textarea);
        onSuccess();
      } catch (e) {
        onFailure();
      }
    }
  }

  copyBtn.addEventListener('click', copyResults);

  /* ---------------- Init ---------------- */
  matrixSection.hidden = true;
  resultsSection.hidden = true;
  calculatingIndicator.hidden = true;
  errorBox.hidden = true;
  largeMatrixNotice.hidden = true;
})();
