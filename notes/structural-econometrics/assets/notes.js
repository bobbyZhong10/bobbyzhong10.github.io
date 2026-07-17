/*
 * From Choices to Counterfactuals — shared script
 * Features: theme (light/dark) apply & toggle, KaTeX auto-render (local assets),
 *           code tab switching, sidebar current-section highlight, prev/next chapter nav.
 * Include (before </body>):
 *   <script defer src="assets/katex/katex.min.js"></script>
 *   <script defer src="assets/katex/auto-render.min.js"></script>
 *   <script defer src="assets/notes.js"></script>
 * Note: each page also has a tiny inline <head> script that applies the theme
 *       before first paint to avoid a flash.
 */

/* ---------- Theme: apply the saved preference as early as possible ---------- */
(function () {
  try {
    var saved = localStorage.getItem('theme');
    if (saved === 'dark' || saved === 'light') {
      document.documentElement.setAttribute('data-theme', saved);
    }
  } catch (e) { /* if localStorage is unavailable, follow the system */ }
})();

/* Chapter order (used by the prev/next navigation) */
var CHAPTER_ORDER = [
  { file: 'index.html',                   label: 'Cover' },
  { file: '01_overview_roadmap.html',     label: '01 Overview & Decision Framework' },
  { file: '02_why_structural.html',       label: '02 Why Structural Models' },
  { file: '03_rum_logit.html',            label: '03 RUM and Logit' },
  { file: '04_gev_mixed_logit.html',      label: '04 Nested and Mixed Logit' },
  { file: '05_simulation_estimation.html',label: '05 Simulation and Numerical Estimation' },
  { file: '06_endogeneity_blp.html',      label: '06 Endogeneity and BLP' },
  { file: '07_dynamic_choice.html',       label: '07 Dynamic Discrete Choice' },
  { file: '08_applications.html',         label: '08 Platform Research Applications' },
  { file: '09_workflow_reference.html',   label: '09 Workflow, Validation, and Reference' }
];

function effectiveTheme() {
  var attr = document.documentElement.getAttribute('data-theme');
  if (attr === 'dark' || attr === 'light') return attr;
  return (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches)
    ? 'dark' : 'light';
}

function buildThemeToggle() {
  var btn = document.createElement('button');
  btn.className = 'theme-toggle';
  btn.setAttribute('type', 'button');
  btn.setAttribute('aria-label', 'Toggle dark / light mode');
  var setIcon = function () { btn.textContent = effectiveTheme() === 'dark' ? '☀' : '☾'; };
  setIcon();
  btn.addEventListener('click', function () {
    var next = effectiveTheme() === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('theme', next); } catch (e) { /* ignore */ }
    setIcon();
  });
  document.body.appendChild(btn);
}

function buildChapterNav() {
  var main = document.querySelector('.main-content');
  if (!main) return;
  var footer = main.querySelector('.page-footer');
  var path = location.pathname.split('/').pop() || 'index.html';
  var idx = -1;
  for (var i = 0; i < CHAPTER_ORDER.length; i++) {
    if (CHAPTER_ORDER[i].file === path) { idx = i; break; }
  }
  if (idx < 1) return; /* cover or unrecognized page: do not inject */

  var nav = document.createElement('nav');
  nav.className = 'chapter-nav';
  nav.setAttribute('aria-label', 'Chapter navigation');

  var prev = CHAPTER_ORDER[idx - 1];
  var next = CHAPTER_ORDER[idx + 1];

  if (prev) {
    var a = document.createElement('a');
    a.className = 'nav-prev';
    a.href = prev.file;
    a.innerHTML = '<div class="nav-dir">← Previous</div><div class="nav-title"></div>';
    a.querySelector('.nav-title').textContent = prev.label;
    nav.appendChild(a);
  }
  if (next) {
    var b = document.createElement('a');
    b.className = 'nav-next';
    b.href = next.file;
    b.innerHTML = '<div class="nav-dir">Next →</div><div class="nav-title"></div>';
    b.querySelector('.nav-title').textContent = next.label;
    nav.appendChild(b);
  }
  if (!nav.children.length) return;

  if (footer) footer.parentNode.insertBefore(nav, footer);
  else main.appendChild(nav);
}

function syncResponsiveTables(mediaQuery) {
  if (!mediaQuery.matches) {
    document.querySelectorAll('.table-scroll[data-mobile-table-scroll]').forEach(function (wrapper) {
      var table = wrapper.querySelector(':scope > table');
      if (table) wrapper.parentNode.insertBefore(table, wrapper);
      wrapper.remove();
    });
    return;
  }

  document.querySelectorAll('.main-content table').forEach(function (table) {
    if (table.parentElement && table.parentElement.classList.contains('table-scroll')) return;
    var wrapper = document.createElement('div');
    wrapper.className = 'table-scroll';
    wrapper.setAttribute('data-mobile-table-scroll', '');
    table.parentNode.insertBefore(wrapper, table);
    wrapper.appendChild(table);
  });
}

function labelRenderedMath() {
  document.querySelectorAll('.katex').forEach(function (formula) {
    var source = formula.querySelector('annotation[encoding="application/x-tex"]');
    if (!source || formula.hasAttribute('aria-label')) return;
    formula.setAttribute('role', 'math');
    formula.setAttribute('aria-label', source.textContent.trim());
  });
}

document.addEventListener('DOMContentLoaded', function () {
  /* ---------- KaTeX auto-render ---------- */
  if (typeof renderMathInElement === 'function') {
    renderMathInElement(document.body, {
      delimiters: [
        { left: '$$', right: '$$', display: true },
        { left: '$', right: '$', display: false }
      ],
      throwOnError: false
    });
  }


  /* Retain an accessible formula label when mobile CSS suppresses KaTeX's
     duplicate MathML layout layer to prevent viewport overflow. */
  labelRenderedMath();

  /* ---------- Theme toggle button ---------- */
  buildThemeToggle();

  /* ---------- Keep wide tables inside the reading column on mobile ---------- */
  var tableMediaQuery = window.matchMedia('(max-width: 768px)');
  var updateResponsiveTables = function () { syncResponsiveTables(tableMediaQuery); };
  updateResponsiveTables();
  if (typeof tableMediaQuery.addEventListener === 'function') {
    tableMediaQuery.addEventListener('change', updateResponsiveTables);
  } else if (typeof tableMediaQuery.addListener === 'function') {
    tableMediaQuery.addListener(updateResponsiveTables);
  }

  /* ---------- Prev / next chapter navigation ---------- */
  buildChapterNav();

  /* ---------- Code tabs ---------- */
  document.querySelectorAll('.code-tabs').forEach(function (tabs) {
    var btns = tabs.querySelectorAll('.tab-btn');
    var panels = tabs.querySelectorAll('.tab-panel');
    btns.forEach(function (btn, i) {
      btn.addEventListener('click', function () {
        btns.forEach(function (b) { b.classList.remove('active'); });
        panels.forEach(function (p) { p.classList.remove('active'); });
        btn.classList.add('active');
        if (panels[i]) panels[i].classList.add('active');
      });
    });
  });

  /* ---------- Sidebar current-section highlight ---------- */
  var tocLinks = Array.prototype.slice.call(
    document.querySelectorAll('.sidebar .toc a[href^="#"]')
  );
  if (tocLinks.length && 'IntersectionObserver' in window) {
    var byId = {};
    tocLinks.forEach(function (a) {
      byId[decodeURIComponent(a.getAttribute('href').slice(1))] = a;
    });
    var current = null;
    var observer = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          var link = byId[entry.target.id];
          if (link) {
            if (current) current.classList.remove('active');
            link.classList.add('active');
            current = link;
          }
        }
      });
    }, { rootMargin: '0px 0px -75% 0px' });
    Object.keys(byId).forEach(function (id) {
      var el = document.getElementById(id);
      if (el) observer.observe(el);
    });
  }
});
