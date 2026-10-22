(function () {
  'use strict';

  var CHAPTER_ORDER = [
    { file: 'index.html', label: 'Cover' },
    { file: 'introduction.html', label: 'Introduction · How to Use' },
    { file: '01_consumer_producer_welfare.html', label: '01 Consumer, Producer & Welfare' },
    { file: '02_game_theory.html', label: '02 Game Theory' },
    { file: '03_information_mechanism_design.html', label: '03 Information & Mechanism Design' },
    { file: '04_pricing_imperfect_competition.html', label: '04 Pricing & Imperfect Competition' },
    { file: '05_m_estimation_nonlinear_regression.html', label: '05 M-Estimation & Nonlinear Regression' },
    { file: '06_maximum_likelihood.html', label: '06 Maximum Likelihood' },
    { file: '07_gmm_minimum_distance.html', label: '07 GMM & Minimum Distance' },
    { file: '08_inference.html', label: '08 Inference: Bootstrap & Clustering' },
    { file: '09_iv_system_estimation.html', label: '09 Instrumental Variables & System Estimation' },
    { file: '10_linear_panel_fe_re.html', label: '10 Linear Panels: FE, RE, FD & CRE' },
    { file: '11_dynamic_panel.html', label: '11 Dynamic Panels' },
    { file: '12_binary_multinomial_response.html', label: '12 Binary & Multinomial Response Models' },
    { file: '13_quantile_regression.html', label: '13 Quantile Regression' },
    { file: '14_nonparametrics.html', label: '14 Nonparametrics & Semiparametrics' },
    { file: '15_potential_outcomes.html', label: '15 Potential Outcomes & Selection' },
    { file: '16_research_design_architecture.html', label: '16 Research Design' },
    { file: '17_measurement_digital_traces.html', label: '17 Measurement & Digital Traces' },
    { file: '18_iv_late_mte.html', label: '18 IV: LATE & MTE' },
    { file: '19_did_event_studies.html', label: '19 DiD & Event Studies' },
    { file: '20_did_beyond_binary.html', label: '20 DiD II: Non-Binary Treatments' },
    { file: '21_regression_discontinuity.html', label: '21 Regression Discontinuity' },
    { file: '22_synthetic_control_factor_models.html', label: '22 Synthetic Control & Panel Factor Models' },
    { file: '23_sensitivity_partial_identification.html', label: '23 Sensitivity & Partial Identification' },
    { file: '24_causal_ml_dml.html', label: '24 Causal ML & DML' },
    { file: '25_digital_experimentation.html', label: '25 Digital Experimentation' },
    { file: '26_interference_marketplace_experiments.html', label: '26 Interference & Marketplace Experiments' },
    { file: '27_discrete_choice.html', label: '27 Discrete Choice & Mixed Logit' },
    { file: '28_demand_blp.html', label: '28 Demand Estimation & BLP' },
    { file: '29_supply_markups_mergers.html', label: '29 Supply, Markups & Mergers' },
    { file: '30_entry_partial_id.html', label: '30 Entry & Partial Identification' },
    { file: '31_dynamic_structural.html', label: '31 Dynamic Structural Models' },
    { file: '32_platforms_networks.html', label: '32 Platforms & Network Effects' },
    { file: '33_search_advertising.html', label: '33 Search, Advertising & Experience Goods' },
    { file: '34_ai_algorithms_data.html', label: '34 AI, Adoption & Data Economics' },
    { file: '35_quick_reference.html', label: '35 Quick Reference & Tools' }
  ];

  function chapterLink(chapter, direction, label) {
    var link = document.createElement('a');
    link.className = 'nav-' + direction;
    link.href = chapter.file;

    var directionText = document.createElement('div');
    directionText.className = 'nav-dir';
    directionText.textContent = label;

    var title = document.createElement('div');
    title.className = 'nav-title';
    title.textContent = chapter.label;

    link.appendChild(directionText);
    link.appendChild(title);
    return link;
  }

  function buildChapterNav() {
    var main = document.querySelector('.main-content');
    if (!main) return;

    var currentFile = location.pathname.split('/').pop() || 'index.html';
    var currentIndex = CHAPTER_ORDER.findIndex(function (chapter) {
      return chapter.file === currentFile;
    });
    if (currentIndex < 1) return;

    var nav = document.createElement('nav');
    nav.className = 'chapter-nav';
    nav.setAttribute('aria-label', 'Chapter navigation');

    if (currentIndex > 0) {
      nav.appendChild(chapterLink(CHAPTER_ORDER[currentIndex - 1], 'prev', 'Previous'));
    }
    if (currentIndex < CHAPTER_ORDER.length - 1) {
      nav.appendChild(chapterLink(CHAPTER_ORDER[currentIndex + 1], 'next', 'Next'));
    }

    var footer = main.querySelector('.page-footer');
    if (footer) footer.parentNode.insertBefore(nav, footer);
    else main.appendChild(nav);
  }

  function buildNotesTopbar() {
    if (document.querySelector('.notes-topbar')) return;
    var topbar = document.createElement('div');
    topbar.className = 'notes-topbar';
    topbar.innerHTML = '<a class="notes-home" href="https://bobbyzhong10.github.io/">← Yunzhou (Bobby) Zhong</a><div class="notes-topbar__right"><a href="https://bobbyzhong10.github.io/notes/">Research Notes</a><button id="theme-toggle" class="notes-theme-toggle" type="button" aria-label="Toggle theme"><span></span></button></div>';
    document.body.insertBefore(topbar, document.body.firstChild);
  }

  function bindThemeToggle() {
    var button = document.getElementById('theme-toggle');
    if (!button || button.dataset.bound) return;
    button.dataset.bound = 'true';
    button.addEventListener('click', function () {
      var next = document.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      try { localStorage.setItem('theme', next); } catch (e) { /* ignore */ }
    });
  }

  function setupMobileSidebar() {
    var sidebar = document.querySelector('.sidebar');
    if (!sidebar || sidebar.querySelector('.notes-nav-toggle')) return;

    var button = document.createElement('button');
    button.className = 'notes-nav-toggle';
    button.type = 'button';
    button.textContent = 'Contents & chapters';
    button.setAttribute('aria-expanded', 'false');
    var anchor = sidebar.querySelector('.sidebar-subtitle') || sidebar.querySelector('.sidebar-title');
    anchor.insertAdjacentElement('afterend', button);

    var media = window.matchMedia('(max-width: 860px)');
    function setCompactMode() {
      if (media.matches) {
        sidebar.classList.add('is-collapsed');
        button.setAttribute('aria-expanded', 'false');
      } else {
        sidebar.classList.remove('is-collapsed');
        button.setAttribute('aria-expanded', 'true');
      }
    }
    button.addEventListener('click', function () {
      var collapsed = sidebar.classList.toggle('is-collapsed');
      button.setAttribute('aria-expanded', String(!collapsed));
    });
    setCompactMode();
    if (media.addEventListener) media.addEventListener('change', setCompactMode);
    else media.addListener(setCompactMode);
  }

  document.addEventListener('DOMContentLoaded', function () {
    buildNotesTopbar();
    bindThemeToggle();
    setupMobileSidebar();
    buildChapterNav();
  });
}());
