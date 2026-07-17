(function () {
  'use strict';

  var CHAPTER_ORDER = [
    { file: 'index.html', label: 'Cover' },
    { file: 'introduction.html', label: 'Introduction · How to Use' },
    { file: '01_consumer_producer_welfare.html', label: '01 Consumer, Producer & Welfare' },
    { file: '02_game_theory.html', label: '02 Game Theory' },
    { file: '03_information_mechanism_design.html', label: '03 Information & Mechanism Design' },
    { file: '04_pricing_imperfect_competition.html', label: '04 Pricing & Imperfect Competition' },
    { file: '05_estimation_mle_gmm.html', label: '05 Estimation: MLE & GMM' },
    { file: '06_inference.html', label: '06 Inference: Bootstrap & Clustering' },
    { file: '07_panel_dynamic_panel.html', label: '07 Panel & Dynamic Panel' },
    { file: '08_nonparametrics.html', label: '08 Nonparametrics & Semiparametrics' },
    { file: '09_potential_outcomes.html', label: '09 Potential Outcomes & Selection' },
    { file: '10_research_design_architecture.html', label: '10 Research Design' },
    { file: '11_measurement_digital_traces.html', label: '11 Measurement & Digital Traces' },
    { file: '12_iv_late_mte.html', label: '12 IV: LATE & MTE' },
    { file: '13_did_event_studies.html', label: '13 DiD & Event Studies' },
    { file: '14_rdd_synthetic_control.html', label: '14 RDD & Synthetic Control' },
    { file: '15_sensitivity_partial_identification.html', label: '15 Sensitivity & Partial Identification' },
    { file: '16_causal_ml_dml.html', label: '16 Causal ML & DML' },
    { file: '17_digital_experimentation.html', label: '17 Digital Experimentation' },
    { file: '18_interference_marketplace_experiments.html', label: '18 Interference & Marketplace Experiments' },
    { file: '19_discrete_choice.html', label: '19 Discrete Choice & Mixed Logit' },
    { file: '20_demand_blp.html', label: '20 Demand Estimation & BLP' },
    { file: '21_supply_markups_mergers.html', label: '21 Supply, Markups & Mergers' },
    { file: '22_entry_partial_id.html', label: '22 Entry & Partial Identification' },
    { file: '23_dynamic_structural.html', label: '23 Dynamic Structural Models' },
    { file: '24_platforms_networks.html', label: '24 Platforms & Network Effects' },
    { file: '25_search_advertising.html', label: '25 Search, Advertising & Information Goods' },
    { file: '26_ai_algorithms_data.html', label: '26 AI, Adoption & Data Economics' },
    { file: '27_quick_reference.html', label: '27 Quick Reference & Tools' }
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

  document.addEventListener('DOMContentLoaded', buildChapterNav);
}());
