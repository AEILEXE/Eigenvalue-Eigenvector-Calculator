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

  var currentN = 0;
  var lastResults = null; // holds { matrix, n, groups } for copy functionality

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
