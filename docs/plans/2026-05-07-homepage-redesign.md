# Homepage Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the left sidebar with a top hero on the home page, remove the sidebar from inner pages, refine the masthead/cards/typography toward a Vercel-docs feel, while preserving all content and the existing light/dark Penn-blue/Penn-red palette.

**Architecture:** Jekyll site (Minimal Mistakes / academicpages base). The bulk of the change is CSS in `assets/css/main.scss`, plus a new `_includes/home-hero.html` partial wired into `_layouts/default.html`. The Susy grid `.page` rule and the sidebar include get adjusted so single-column full-width layout works on every page. No content (markdown bodies) changes except renaming `# CV` → `# Curriculum Vitae` and the Education page H2.

**Tech Stack:** Jekyll, Liquid, Sass/SCSS, vanilla JS for theme + masthead scroll shadow. Local server: `bundle exec jekyll serve` (or the existing `run_server.sh`).

**Reference design:** `docs/plans/2026-05-07-homepage-redesign-design.md`

---

## Verification approach

This is a static-site visual change, not unit-tested code. Each task ends with a verification step that runs the local server and visually checks the result in both themes at the relevant breakpoints. The local server runs in the background:

```bash
bundle exec jekyll serve --host 127.0.0.1 --port 4000
```

For each visual task:
1. Save edits.
2. Wait for Jekyll to regenerate (watch terminal for "regenerated").
3. Refresh `http://127.0.0.1:4000` in browser.
4. Toggle theme button, verify both light + dark.
5. At minimum check viewport widths: 1440 (desktop), 1024 (laptop), 768 (tablet), 414 (phone).

If a task is text-only (markdown / yaml), the verification is just "page loads and the new label appears."

Commit after each task or each tightly related cluster of tasks.

---

## Phase 0: Setup

### Task 0.1: Start the local server

**Step 1:** Start the Jekyll dev server in the background.

```bash
bundle exec jekyll serve --host 127.0.0.1 --port 4000
```

**Step 2:** Open `http://127.0.0.1:4000` in the browser. Confirm the current site renders, both themes work, and you can navigate between pages.

**Step 3:** No commit (no code changed).

---

## Phase 1: Low-risk text changes

### Task 1.1: Rename "CV" to "Curriculum Vitae" in nav and on the page

**Files:**
- Modify: `_data/navigation.yml`
- Modify: `_pages/cv.md`

**Step 1:** In `_data/navigation.yml`, change the last menu item.

Current:
```yaml
  - title: "CV"
    url: "/cv/"
```

New:
```yaml
  - title: "Curriculum Vitae"
    url: "/cv/"
```

**Step 2:** In `_pages/cv.md`, change the H1.

Current:
```markdown
# CV
```

New:
```markdown
# Curriculum Vitae
```

URL frontmatter (`permalink: /cv/`) is unchanged.

**Step 3:** Verify in browser. Refresh, click "Curriculum Vitae" in the nav. The CV page H1 reads "Curriculum Vitae" and the URL is still `/cv/`.

**Step 4:** Commit.

```bash
git add _data/navigation.yml _pages/cv.md
git commit -m "Rename CV nav label and page heading to Curriculum Vitae"
```

### Task 1.2: Rename Education inner H2 to "Academic Background"

**Files:**
- Modify: `_pages/education.md`

**Step 1:** In `_pages/education.md`, the first content card H2 changes.

Current:
```markdown
<section class="content-card" markdown="1">
## Education
```

New:
```markdown
<section class="content-card" markdown="1">
## Academic Background
```

The page H1 (`# Education`) and the second H2 (`## PhD Level Courses`) stay untouched.

**Step 2:** Verify in browser. The Education page now shows H1 "Education", first card H2 "Academic Background", second card H2 "PhD Level Courses".

**Step 3:** Commit.

```bash
git add _pages/education.md
git commit -m "Rename Education inner heading to Academic Background"
```

---

## Phase 2: Masthead refinements

### Task 2.1: Drop uppercase + tune spacing on nav links

**Files:**
- Modify: `assets/css/main.scss` (around line 349–356)

**Step 1:** Locate the `.greedy-nav .visible-links li a` rule. Change three properties.

Current:
```scss
.greedy-nav .visible-links li a {
    margin: 0 0.7rem;
    padding: 0.55rem 0;
    font-size: 0.9rem;
    font-weight: 600;
    letter-spacing: 0.02em;
    text-transform: uppercase;
}
```

New:
```scss
.greedy-nav .visible-links li a {
    margin: 0 0.7rem;
    padding: 0.5rem 0;
    font-size: 0.92rem;
    font-weight: 550;
    letter-spacing: 0.005em;
}
```

(The `text-transform: uppercase;` line is removed entirely.)

**Step 2:** Also bump the masthead inner padding slightly (line 344–347).

Current:
```scss
.masthead__inner-wrap {
    padding-top: 0.55rem;
    padding-bottom: 0.55rem;
}
```

New:
```scss
.masthead__inner-wrap {
    padding-top: 0.5rem;
    padding-bottom: 0.5rem;
}
```

**Step 3:** Verify in browser. Menu now reads "Home / Research and Projects / Education / Awards and Services / Teaching / Curriculum Vitae" in title case. Hover state still shows the underline animation. Test both themes and 1024px / 1440px viewports.

**Step 4:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Drop uppercase and tighten spacing in masthead nav"
```

### Task 2.2: Refine active-state visual

**Files:**
- Modify: `assets/css/main.scss` (around line 358–382)

**Step 1:** Replace the heavy red rounded-rect indicator with a thin gradient underline plus a soft glow.

Current block (lines 358–382 area):
```scss
.masthead__menu-item--active a {
    font-weight: 700;
}

.masthead__menu-item--active {
    position: relative;
}

.masthead__menu-item--active::after {
    content: "";
    position: absolute;
    left: 0.55rem;
    right: 0.55rem;
    bottom: 0.42rem;
    height: 0.3rem;
    border-radius: 999px;
    background: color-mix(in srgb, var(--site-accent-soft) 75%, transparent);
    z-index: -1;
}

.masthead__menu-item--active a:before {
    transform: scaleX(1) !important;
    height: 2px !important;
    background: var(--site-accent) !important;
}
```

New block:
```scss
.masthead__menu-item--active a {
    font-weight: 650;
}

.masthead__menu-item--active {
    position: relative;
}

.masthead__menu-item--active a:before {
    transform: scaleX(1) !important;
    height: 2px !important;
    border-radius: 999px;
    background: linear-gradient(
        90deg,
        var(--site-link),
        var(--site-accent)
    ) !important;
    box-shadow: 0 0 8px color-mix(in srgb, var(--site-link) 35%, transparent);
}
```

(The `::after` red pill is removed; the existing `:before` underline is repurposed with a gradient + glow.)

**Step 2:** Also reduce the hover underline thickness for non-active links. Find the `.greedy-nav .visible-links a` styles in `_sass/_navigation.scss` (lines ~226–249). The `&:before` pseudo currently has `height: 4px`. Override it in `assets/css/main.scss` so we don't touch the vendored sass:

Add this block in `assets/css/main.scss` near the active-state rules:

```scss
.greedy-nav .visible-links a:before {
    height: 2px !important;
    border-radius: 999px;
    background: linear-gradient(90deg, var(--site-link), var(--site-accent)) !important;
}
```

**Step 3:** Verify in browser. Click each nav item. The active item shows a thin gradient underline with a faint glow; no red box anymore. Hover on inactive items shows the underline expanding from scale-0 to scale-1, also gradient and 2px thick. Both themes.

**Step 4:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Replace heavy red active state with gradient underline + glow"
```

### Task 2.3: Update theme toggle to ☼ / ☾ glyphs

**Files:**
- Modify: `_includes/scripts.html` (the inline theme script)
- Modify: `_includes/masthead.html` (initial icon)
- Modify: `assets/css/main.scss` (icon styling for the new glyph)

**Step 1:** In `_includes/masthead.html`, change the initial icon span content.

Current (line 21):
```html
<span id="theme-toggle-icon" aria-hidden="true">☾</span>
```

New:
```html
<span id="theme-toggle-icon" aria-hidden="true">☼</span>
```

(Initial render is light mode, so the icon shows ☼ — meaning "currently sun / light".)

**Step 2:** In `_includes/scripts.html`, the JS currently does `icon.textContent = '';`. Replace that line so the icon glyph reflects the active theme.

Current `setTheme` function:
```javascript
function setTheme(theme) {
  root.setAttribute('data-theme', theme);
  localStorage.setItem('theme', theme);
  btn.setAttribute('data-theme', theme);
  icon.textContent = '';
  btn.setAttribute('aria-label', theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode');
}
```

New `setTheme` function:
```javascript
function setTheme(theme) {
  root.setAttribute('data-theme', theme);
  localStorage.setItem('theme', theme);
  btn.setAttribute('data-theme', theme);
  icon.textContent = theme === 'dark' ? '☾' : '☼';
  btn.setAttribute('aria-label', theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode');
}
```

**Step 3:** In `assets/css/main.scss`, the existing `#theme-toggle-icon::before` rule (around line 310–317) renders a small solid dot inside the indicator. With a real glyph now in `textContent`, the `::before` dot would conflict. Delete the `#theme-toggle-icon::before { ... }` block entirely.

Also update the icon size so the glyph fits nicely. Current rule (around line 296–308):

```scss
#theme-toggle-icon {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 1.05rem;
    height: 1.05rem;
    border-radius: 999px;
    background: linear-gradient(180deg, #ffffff, #e8edf7);
    color: #011f5b;
    box-shadow: 0 2px 8px rgba(1, 31, 91, 0.18);
    transform: translateX(0);
    transition: transform 0.2s ease, background 0.2s ease, color 0.2s ease;
}
```

Add `font-size: 0.78rem; line-height: 1; font-weight: 600;` and remove the `align-items` value if needed. Keep transform behavior. New version:

```scss
#theme-toggle-icon {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 1.15rem;
    height: 1.15rem;
    border-radius: 999px;
    background: linear-gradient(180deg, #ffffff, #e8edf7);
    color: #011f5b;
    box-shadow: 0 2px 8px rgba(1, 31, 91, 0.18);
    font-size: 0.78rem;
    line-height: 1;
    font-weight: 600;
    transform: translateX(0);
    transition: transform 0.2s ease, background 0.2s ease, color 0.2s ease;
}
```

The dark variant (around line 319–324) needs its width matched too: `transform: translateX(0.98rem);` (replaces `1.08rem`) so the dark indicator sits flush with the right edge of the slightly larger pill.

```scss
html[data-theme="dark"] #theme-toggle-icon {
    transform: translateX(0.98rem);
    background: linear-gradient(180deg, #1f2d40, #101926);
    color: #f3d7d7;
    box-shadow: 0 2px 10px rgba(0, 0, 0, 0.34);
}
```

Also widen the toggle pill itself so the glyph + slide range fit. Update `.theme-toggle-btn` (around line 257–275): change `width: 2.7rem;` to `width: 2.85rem;`.

**Step 4:** Verify in browser. Refresh in light mode → icon shows ☼ on the left, with light gradient. Click toggle → icon slides to the right, glyph becomes ☾ with dark gradient. Click again → returns. Test in both desktop and mobile widths.

**Step 5:** Commit.

```bash
git add _includes/masthead.html _includes/scripts.html assets/css/main.scss
git commit -m "Switch theme toggle glyph to ☼/☾ characters"
```

### Task 2.4: Add scroll-state shadow to masthead

**Files:**
- Modify: `_includes/scripts.html` (or `assets/js/_main.js` — choose `scripts.html` since it already has inline JS)
- Modify: `assets/css/main.scss`

**Step 1:** In `assets/css/main.scss`, add a `.is-scrolled` modifier on `.masthead`. Find the existing `.masthead` rule (around line 203–208) and after it add:

```scss
.masthead {
    transition: box-shadow 0.18s ease, background-color 0.2s ease;
}

.masthead.is-scrolled {
    box-shadow:
        0 1px 0 color-mix(in srgb, var(--site-border) 75%, transparent),
        0 6px 24px color-mix(in srgb, var(--site-link) 12%, transparent);
}
```

**Step 2:** In `_includes/scripts.html`, append a new IIFE that toggles the class. Add it after the existing theme IIFE (before the analytics include).

```html
<script>
  (function() {
    var masthead = document.querySelector('.masthead');
    if (!masthead) return;
    function update() {
      if (window.scrollY > 4) {
        masthead.classList.add('is-scrolled');
      } else {
        masthead.classList.remove('is-scrolled');
      }
    }
    update();
    window.addEventListener('scroll', update, { passive: true });
  })();
</script>
```

**Step 3:** Verify. Load page → masthead has no shadow at top. Scroll a bit → soft shadow appears, fades back when scrolled to top. Both themes.

**Step 4:** Commit.

```bash
git add _includes/scripts.html assets/css/main.scss
git commit -m "Add scroll-state shadow to masthead for docs feel"
```

---

## Phase 3: Hero block + sidebar removal

### Task 3.1: Create the home hero include

**Files:**
- Create: `_includes/home-hero.html`

**Step 1:** Create the new file with the hero markup. The data comes from `site.author` (already in `_config.yml`).

```html
<section class="home-hero" aria-label="About">
  <div class="home-hero__inner">
    <div class="home-hero__avatar">
      <img src="{{ site.author.avatar | relative_url }}" alt="{{ site.author.name }}">
    </div>
    <div class="home-hero__content">
      <h1 class="home-hero__name">{{ site.author.name }}</h1>
      <ul class="home-hero__meta">
        {% if site.author.location %}
          <li><i class="fa fa-fw fa-map-marker" aria-hidden="true"></i> {{ site.author.location }}</li>
        {% endif %}
        {% if site.author.employer %}
          <li><i class="fa fa-fw fa-university" aria-hidden="true"></i> {{ site.author.employer }}</li>
        {% endif %}
      </ul>
      <ul class="home-hero__links" aria-label="Profiles and contact">
        {% if site.author.email %}
          <li><a href="mailto:{{ site.author.email }}" aria-label="Email">
            <i class="fas fa-fw fa-envelope" aria-hidden="true"></i>
          </a></li>
        {% endif %}
        {% if site.author.linkedin %}
          <li><a href="https://www.linkedin.com/in/{{ site.author.linkedin }}" aria-label="LinkedIn" target="_blank" rel="noopener">
            <i class="fab fa-fw fa-linkedin" aria-hidden="true"></i>
          </a></li>
        {% endif %}
        {% if site.author.googlescholar %}
          <li><a href="{{ site.author.googlescholar }}" aria-label="Google Scholar" target="_blank" rel="noopener">
            <i class="fas fa-fw fa-graduation-cap" aria-hidden="true"></i>
          </a></li>
        {% endif %}
        {% if site.author.orcid %}
          <li><a href="{{ site.author.orcid }}" aria-label="ORCID" target="_blank" rel="noopener">
            <i class="ai ai-orcid-square ai-fw" aria-hidden="true"></i>
          </a></li>
        {% endif %}
      </ul>
    </div>
  </div>
  <div class="home-hero__divider" aria-hidden="true"></div>
</section>
```

**Step 2:** No verification yet — the include is not rendered anywhere.

**Step 3:** Commit.

```bash
git add _includes/home-hero.html
git commit -m "Add home-hero include scaffolding"
```

### Task 3.2: Wire the hero into the home page only

**Files:**
- Modify: `_layouts/default.html`

**Step 1:** Open `_layouts/default.html`. After the masthead include and just before `<div id="main">`, add a Liquid conditional that renders the hero only when the current page is the home (`/`).

Current relevant block:
```html
{% include browser-upgrade.html %}
{% include masthead.html %}

<div id="main" role="main">
  {% include sidebar.html %}
```

New block:
```html
{% include browser-upgrade.html %}
{% include masthead.html %}

{% if current_url == "/" %}
  {% include home-hero.html %}
{% endif %}

<div id="main" role="main">
  {% include sidebar.html %}
```

(`current_url` is already assigned earlier in the layout.)

**Step 2:** Remove the now-redundant `# Home` H1 from `_pages/about.md`.

Current:
```markdown
---
permalink: /
title: ""
excerpt: ""
author_profile: true
redirect_from:
  - /about/
  - /about.html
---

# Home

<section class="home-intro-card" markdown="1">
```

New:
```markdown
---
permalink: /
title: ""
excerpt: ""
author_profile: true
redirect_from:
  - /about/
  - /about.html
---

<section class="home-intro-card" markdown="1">
```

**Step 3:** Verify. Reload the home page — the hero block now appears (currently unstyled or partially styled by inherited rules) above the intro paragraph. The "# Home" heading is gone. Inner pages are unaffected.

**Step 4:** Commit.

```bash
git add _layouts/default.html _pages/about.md
git commit -m "Render hero on home page and drop redundant Home heading"
```

### Task 3.3: Style the hero

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Append a new `Hero` section to `main.scss` (place it before the masthead-overrides for cohesion, or at end of the file with the other home rules — end is fine). Use the variables already defined.

```scss
/* ==========================================================================
   HOME HERO
   ========================================================================== */

.home-hero {
    @include container;
    margin: 1.6rem auto 0.4rem;
    padding: 0 1em;

    @include breakpoint($x-large) {
        max-width: $x-large;
    }
}

.home-hero__inner {
    display: flex;
    flex-direction: row;
    align-items: center;
    gap: 1.5rem;
}

.home-hero__avatar {
    flex: 0 0 auto;
}

.home-hero__avatar img {
    display: block;
    width: 130px;
    height: 130px;
    object-fit: cover;
    border-radius: 18px;
    border: 2px solid color-mix(in srgb, var(--site-link) 38%, var(--glass-border-strong));
    box-shadow:
        0 14px 32px color-mix(in srgb, var(--site-link) 14%, transparent),
        0 1px 0 var(--glass-highlight) inset;
    background: var(--glass-bg-strong);
}

.home-hero__content {
    min-width: 0;
}

.home-hero__name {
    margin: 0 0 0.45rem 0;
    padding: 0;
    border: 0 !important;
    font-size: 1.8rem;
    font-weight: 700;
    line-height: 1.15;
    color: var(--site-text);
}

.home-hero__name::after {
    display: none !important;
}

.home-hero__meta {
    list-style: none;
    margin: 0 0 0.7rem 0;
    padding: 0;
    display: flex;
    flex-wrap: wrap;
    gap: 0.4rem 1rem;
    font-family: "JetBrains Mono", "SF Mono", "Cascadia Code", ui-monospace, monospace;
    font-size: 0.84rem;
    color: var(--site-muted);
}

.home-hero__meta li {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    white-space: nowrap;
}

.home-hero__meta i {
    color: var(--site-link);
}

.home-hero__links {
    list-style: none;
    margin: 0;
    padding: 0;
    display: flex;
    align-items: center;
    gap: 0.55rem;
}

.home-hero__links a {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 2.05rem;
    height: 2.05rem;
    border-radius: 10px;
    border: 1px solid color-mix(in srgb, var(--site-link) 22%, var(--glass-border));
    background: var(--glass-bg-strong);
    color: var(--site-link) !important;
    text-decoration: none;
    transition: transform 0.18s ease, border-color 0.18s ease, color 0.18s ease, background-color 0.18s ease;
}

.home-hero__links a:hover {
    color: var(--site-link-hover) !important;
    border-color: color-mix(in srgb, var(--site-link) 50%, var(--glass-border-strong));
    transform: translateY(-1px);
}

.home-hero__links i,
.home-hero__links .fa,
.home-hero__links .fab,
.home-hero__links .fas,
.home-hero__links .ai {
    color: inherit !important;
    font-size: 0.95rem;
}

.home-hero__divider {
    margin: 1.4rem 0 0;
    height: 2px;
    width: 4rem;
    border-radius: 999px;
    background: linear-gradient(90deg, var(--site-link), var(--site-accent));
    opacity: 0.78;
}

@include breakpoint(max-width $small) {
    .home-hero__inner {
        flex-direction: column;
        align-items: center;
        text-align: center;
        gap: 1.1rem;
    }

    .home-hero__avatar img {
        width: 96px;
        height: 96px;
        border-radius: 14px;
    }

    .home-hero__name {
        font-size: 1.6rem;
    }

    .home-hero__meta {
        justify-content: center;
        font-size: 0.8rem;
    }

    .home-hero__links {
        justify-content: center;
    }

    .home-hero__divider {
        margin-left: auto;
        margin-right: auto;
    }
}
```

**Step 2:** Verify in browser at 1440 → 1024 → 768 → 414 widths. Hero looks balanced: photo on left, name/meta/icons on right. Below the hero a short gradient divider sits before the intro paragraph. Toggle dark mode — hero photo ring, meta text, and icon backgrounds all read clearly. Hover icon backgrounds → lift + border darken.

**Step 3:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Style the home hero block with avatar, meta, and icon links"
```

### Task 3.4: Remove the sidebar from all pages and let main go full-width

**Files:**
- Modify: `_layouts/default.html`
- Modify: `assets/css/main.scss`

**Step 1:** In `_layouts/default.html`, delete the sidebar include line.

Current:
```html
<div id="main" role="main">
  {% include sidebar.html %}

  <article class="page" itemscope itemtype="http://schema.org/CreativeWork">
```

New:
```html
<div id="main" role="main">
  <article class="page" itemscope itemtype="http://schema.org/CreativeWork">
```

**Step 2:** In `assets/css/main.scss`, override the Susy `.page` width so it takes the full container. Append:

```scss
/* ==========================================================================
   FULL-WIDTH MAIN COLUMN
   ========================================================================== */

.page {
    float: none !important;
    width: 100% !important;
    margin: 0 auto !important;
    padding: 0 !important;
    max-width: 880px;
}

#main {
    max-width: 980px;
}

@media (min-width: 1440px) {
    #main {
        max-width: 1040px;
    }

    .page {
        max-width: 920px;
    }
}
```

(The existing `@media (min-width: 1440px) #main { max-width: 1360px; }` rule near the bottom of `main.scss` should be deleted to avoid conflict — it was tuned for the two-column layout. Locate it and remove.)

**Step 3:** Verify each page:
- Home: hero + intro + PhD callout + News + Site Visitors all stacked, no left empty column.
- Research / Education / Awards / Teaching / Curriculum Vitae: H1 + content cards, centered, no left empty column.
- 1024 / 1440 / 768 / 414 widths.
- Light + dark.

**Step 4:** Commit.

```bash
git add _layouts/default.html assets/css/main.scss
git commit -m "Remove sidebar and let main column run full width"
```

### Task 3.5: Clean up the now-orphaned profile_box rules and includes

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** The `_includes/sidebar.html` and `_includes/author-profile.html` partials are no longer rendered anywhere. Leave them on disk (no behavior cost, easier to revert) but delete the CSS that targeted `.profile_box`, `.author__urls li`, and `.sidebar` overrides we wrote in `main.scss`. Search for `profile_box`, `.author__urls`, and `.sidebar` and remove the bespoke rules added on top of the base theme. Specifically remove these blocks (line numbers approximate):

- `.profile_box { ... }` (~line 406–420)
- `.profile_box::before { ... }` (~line 422–432)
- `.profile_box::after { ... }` (~line 434–441)
- `.author__avatar img { ... }` (~line 443–446)
- `.author__urls li { ... }` (~line 448–455)
- `.author__urls li > i { ... }` (~line 457–461)
- `.author__urls li > a { ... }` (~line 463–467)
- `.author__urls li > a i { ... }` (~line 469–474)
- `.sidebar .author__urls li { ... }` (~line 476–478)
- `.sidebar .author__name { ... }` (~line 480–486)
- `.sidebar .author__meta-line { ... }` (~line 488–493)
- `.sidebar .author__meta-line--employer { ... }` (~line 495–498)
- `.author__urls .author__meta-line i { ... }` (~line 500–502)
- `.author__urls a i, .author__urls a .fa, ... { ... }` (~line 504–510)
- `.author__urls a:hover i, ... { ... }` (~line 512–518)
- The `@include breakpoint($large) { .profile_box .author__... }` block (~line 520–538)
- The `.author__bio { color: ... }` rule (~line 677–680)
- The `.author__content, .author__urls-wrapper, .author__urls { background: transparent; }` block (~line 682–686)
- The `.author__urls { background: transparent; }` standalone block (~line 688–690)
- The `.author__urls i, ... { color: var(--site-text) !important; }` block (~line 331–342)

**Step 2:** Also remove the bespoke `body, .greedy-nav, .page, .sidebar, .author__urls-wrapper { ... }` rule (~line 167–174) and replace with just `body, .greedy-nav, .page { ... }`.

**Step 3:** Verify all pages still render correctly (no visual regressions; the rules removed only targeted now-hidden elements). Confirm light + dark.

**Step 4:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Drop orphaned profile_box and sidebar styles"
```

---

## Phase 4: Card and typography refinement

### Task 4.1: Soften card glass strength

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Reduce blur and shadow strength globally on cards. Find `:root` (~line 111–135) and update the relevant variables.

Current:
```scss
--site-shadow: 0 18px 45px rgba(1, 31, 91, 0.08);
--glass-shadow: 0 20px 46px rgba(1, 31, 91, 0.11), 0 1px 0 rgba(255, 255, 255, 0.72) inset;
--glass-blur: 13px;
```

New:
```scss
--site-shadow: 0 14px 36px rgba(1, 31, 91, 0.07);
--glass-shadow: 0 14px 34px rgba(1, 31, 91, 0.09), 0 1px 0 rgba(255, 255, 255, 0.72) inset;
--glass-blur: 11px;
```

Also update the dark variant (~line 137–159):

Current:
```scss
--site-shadow: 0 24px 50px rgba(0, 0, 0, 0.32);
--glass-shadow: 0 22px 50px rgba(0, 0, 0, 0.44), 0 1px 0 rgba(220, 232, 255, 0.16) inset;
--glass-blur: 15px;
```

New:
```scss
--site-shadow: 0 18px 38px rgba(0, 0, 0, 0.28);
--glass-shadow: 0 16px 36px rgba(0, 0, 0, 0.36), 0 1px 0 rgba(220, 232, 255, 0.16) inset;
--glass-blur: 12px;
```

**Step 2:** Reduce the gap between successive content cards. Find `.content-card` (~line 765–777). Change `margin: 1.25rem 0 0;` to `margin: 1.1rem 0 0;`. Find the news-list grid gap rule (`.news-list { display: grid; gap: 0.85rem; }`) and leave at `0.85rem` (already tight).

**Step 3:** Verify in browser. Cards still feel dimensional but less heavy. Check dark mode — shadows should not feel muddy.

**Step 4:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Soften card glass blur and shadow for docs feel"
```

### Task 4.2: Tune section heading H2 size and tracking

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Find `.page__content h2` (~line 389–395). Change font-size from `1.25rem` to `1.18rem`. Add `letter-spacing: 0.005em;`.

```scss
.page__content h2 {
    margin-top: 1.85rem;
    margin-bottom: 0.7rem;
    padding-bottom: 0.28rem;
    border-bottom: 1px solid var(--site-border);
    font-size: 1.18rem;
    letter-spacing: 0.005em;
}
```

**Step 2:** Verify on Research, Education, Awards, Teaching pages. H2s feel slightly smaller and more refined.

**Step 3:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Refine section heading H2 size and letter-spacing"
```

### Task 4.3: Switch news date pill to monospace, drop uppercase

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Find `.news-item__date` (~line 647–660). Change typography.

Current:
```scss
.news-item__date {
    display: inline-flex;
    align-items: center;
    width: fit-content;
    padding: 0.18rem 0.5rem;
    border: 1px solid color-mix(in srgb, var(--site-link) 22%, transparent);
    border-radius: 999px;
    background: color-mix(in srgb, var(--site-link) 10%, var(--site-surface));
    color: color-mix(in srgb, var(--site-link) 88%, var(--site-text));
    font-size: 0.82rem;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}
```

New:
```scss
.news-item__date {
    display: inline-flex;
    align-items: center;
    width: fit-content;
    padding: 0.18rem 0.55rem;
    border: 1px solid color-mix(in srgb, var(--site-link) 22%, transparent);
    border-radius: 999px;
    background: color-mix(in srgb, var(--site-link) 10%, var(--site-surface));
    color: color-mix(in srgb, var(--site-link) 88%, var(--site-text));
    font-family: "JetBrains Mono", "SF Mono", "Cascadia Code", ui-monospace, monospace;
    font-size: 0.78rem;
    font-weight: 600;
    letter-spacing: 0;
}
```

**Step 2:** Verify. News list shows dates like "Apr 2026", "Mar 2026" in monospace, no uppercase, well-aligned, both themes.

**Step 3:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Style news date pill with monospace digits, no uppercase"
```

### Task 4.4: CV "Last updated" treatment

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Find `.cv-updated` (~line 561–566). Add monospace.

Current:
```scss
.cv-updated {
    margin-top: 0.55rem;
    color: var(--site-text);
    opacity: 0.86;
    font-size: 0.9rem;
}
```

New:
```scss
.cv-updated {
    margin-top: 0.55rem;
    color: var(--site-muted);
    font-family: "JetBrains Mono", "SF Mono", "Cascadia Code", ui-monospace, monospace;
    font-size: 0.82rem;
}
```

**Step 2:** Verify on the Curriculum Vitae page.

**Step 3:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Style CV last-updated line with muted monospace"
```

### Task 4.5: Polish the PhD callout border to coordinate with hero ring

**Files:**
- Modify: `assets/css/main.scss`

**Step 1:** Find `.phd-callout` (~line 692–703). Increase the border tone slightly to match the hero photo ring.

Current border-color line:
```scss
border: 1px solid color-mix(in srgb, var(--site-link) 20%, var(--glass-border-strong));
```

New:
```scss
border: 1px solid color-mix(in srgb, var(--site-link) 30%, var(--glass-border-strong));
border-left: 3px solid color-mix(in srgb, var(--site-link) 70%, var(--site-accent));
```

(The existing `::before` gradient bar can stay — visually layered with the new left rule, it gives the callout a clearer "noteworthy" mark.)

**Step 2:** Verify on the home page. The PhD callout reads clearly as a callout, with a left accent rule plus the existing top gradient.

**Step 3:** Commit.

```bash
git add assets/css/main.scss
git commit -m "Refine PhD callout border to coordinate with hero ring"
```

---

## Phase 5: Responsive verification and dark-mode pass

### Task 5.1: Verify each breakpoint on each page

**Step 1:** With the dev server running, walk through:

| Page                  | 1440 | 1280 | 1024 | 768 | 414 | 375 |
|-----------------------|------|------|------|-----|-----|-----|
| Home                  |      |      |      |     |     |     |
| Research and Projects |      |      |      |     |     |     |
| Education             |      |      |      |     |     |     |
| Awards and Services   |      |      |      |     |     |     |
| Teaching              |      |      |      |     |     |     |
| Curriculum Vitae      |      |      |      |     |     |     |

For each cell: load page, check no horizontal scroll, no overlapping text, hero stacks correctly on mobile, masthead nav wraps gracefully (the existing `greedy-nav` handles this — verify the overflow menu opens on the smallest widths).

**Step 2:** If any breakpoint breaks, add a targeted `@include breakpoint(...)` rule in `assets/css/main.scss` next to the relevant component. Likely candidates:
- Hero on tablet (768–1023): may need `gap` reduction.
- Long names/links in the masthead at 768–900px: may need smaller font-size override.

**Step 3:** Commit any fixes.

```bash
git add assets/css/main.scss
git commit -m "Tighten responsive behavior at <breakpoint>"
```

### Task 5.2: Dark mode contrast pass

**Step 1:** Toggle to dark mode. Walk every page again. Check:
- Hero photo ring is visible against dark background.
- Hero meta monospace text passes legibility (`var(--site-muted)` is light in dark mode — should be fine).
- Hero icon link backgrounds (`--glass-bg-strong`) stay distinguishable from page background.
- News date pill text contrast.
- Active nav underline glow is visible but not blown out.
- Card top gradient hairline is visible.

**Step 2:** Adjust any rule that fails by editing the relevant CSS variable (in `:root` or `html[data-theme="dark"]`) rather than hard-coding. Commit any changes.

```bash
git add assets/css/main.scss
git commit -m "Dark-mode contrast adjustments"
```

### Task 5.3: Cross-browser sanity check

**Step 1:** Open the local site in Safari (macOS) and Firefox. Specifically check:
- Backdrop-filter renders on Safari (`-webkit-backdrop-filter` is already in the stylesheet — verify no obvious gaps).
- Theme toggle slide animation runs smoothly.
- Masthead scroll-shadow toggles correctly.

**Step 2:** No commit unless adjustments needed.

### Task 5.4: Final visual review

**Step 1:** Walk through every page once more in the user's preferred browser (Chrome by default), both themes, at 1440 and 414. Confirm the design matches the intent: simple, elegant, restrained, with subtle tech accents (monospace dates, gradient underline, scroll shadow, hero gradient divider).

**Step 2:** Stop the dev server.

**Step 3:** No commit. Implementation is complete; the worktree branch has a clean linear history from Task 1.1 onward.

---

## Closeout

When all phases are green, the worktree branch should be ready to merge. Use the superpowers:finishing-a-development-branch skill to decide how to integrate (squash, fast-forward, or PR).
