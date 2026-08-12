# Personal Academic Homepage — Layout & Visual Refinement

Date: 2026-05-07
Owner: Bobby Zhong
Status: Approved (design phase complete)

## Goal

Refine the layout and visual treatment of the personal academic homepage. Keep
content unchanged. Preserve the existing light/dark themes and the
glass-morphism + Penn-blue/Penn-red palette. Move closer to a Vercel-docs feel:
simple, elegant, restrained, with a subtle tech accent.

## Non-goals

- Do not change page content (intro paragraph, news items, papers, etc.).
- Do not introduce new pages or remove existing pages.
- Do not change URL paths (`/cv/` stays `/cv/` even though the label changes).
- Do not add scroll progress bars or scanning-line decorations (keep restraint).

## Scope

Five blocks of changes, all approved with the user.

### 1. Masthead (top navigation)

- Remove `text-transform: uppercase` from the visible nav links.
- Reduce font size slightly (`0.9rem` → `0.92rem`), font weight (`600` → `550`),
  and letter-spacing (`0.02em` → `0.005em`).
- Replace the heavy red rounded-rect active state with a thin 2px gradient
  (Penn blue → Penn red) underline plus a soft glow.
- Reduce hover underline thickness from 4px to 2px; keep the scaleX animation.
- Tighten masthead vertical padding (`0.55rem` → `0.5rem`).
- Theme toggle: keep the pill slider, but render explicit `☼` / `☾` characters
  inside the indicator instead of an abstract dot.
- Add a subtle scroll-state shadow to the masthead: no shadow at top of page,
  small shadow once `window.scrollY > 4`.

### 2. Home page restructure

New top-down order on the home page:

1. Masthead.
2. Hero block (replaces the left sidebar).
3. Intro paragraph (content unchanged).
4. PhD-seeking callout (content unchanged, mild style polish).
5. News section (content unchanged).
6. Site Visitors section (content unchanged, lighter card).

Hero block details:
- Photo: square with 16px corner radius, ~130×130px on desktop, thin colored
  ring border, soft shadow. Square + rounded corners chosen to match the
  existing rounded-rectangle card language; circular would feel disconnected.
- Photo can sit on either left or right; final placement decided in
  implementation based on visual balance.
- Name: H1-sized, sans, normal title case.
- Meta line(s): location (Philadelphia, PA) and employer (University of
  Pennsylvania) with their existing icons; the location/employer line may be
  rendered in a subtle monospace if it fits visually.
- Social row: icon-only links (Email, LinkedIn, Google Scholar, ORCID) — no
  text labels, just the icons that click through to the destinations.
- Hero is NOT wrapped in a glass card. It reads as page identity, not as
  content. A short gradient divider (Penn blue → Penn red, ~4rem wide) sits
  between the hero and the intro paragraph as a visual anchor.

Removed from the home page:
- The current `# Home` H1.
- The left sidebar / `profile_box` (its data is absorbed into the hero).

Preserved from the home page:
- Intro paragraph (content unchanged).
- PhD callout (mild polish only — left border accent now coordinates with the
  hero ring color).
- News list — but date pills lose `text-transform: uppercase` and the
  `0.08em` letter-spacing, switching to monospace digits for an even,
  consistent look.
- Visitor map — kept in place but the wrapping card is lighter.

NOT added (explicitly rejected during brainstorming):
- `01 ───` numbered prefix on section headings.
- Role line ("Master Student · UPenn MBDS") under the name.
- Text labels on social icons.

### 3. Inner pages

All inner pages (Research, Education, Awards, Teaching, Curriculum Vitae)
remove the left sidebar and become single-column, full-width.

Common treatment:
- Drop the left sidebar / `profile_box` entirely.
- Center the main content with a comfortable reading max-width
  (~800–880px on desktop, scaling up to ~960px on very wide screens).
- Add a small breathing space above the H1.

Per-page adjustments:
- **Research and Projects**: keep the MIS Quarterly research-quote card and the
  Research Interests card. List items in Working Papers / Conference
  Proceedings / Presentations get a subtle visual layer for venue / year using
  muted color or a small chip.
- **Education**: rename the inner H2 from "Education" (which duplicates the
  page H1) to "Academic Background". The "PhD Level Courses" H2 stays as is.
- **Awards and Services**: list items get a subtle year chip on the right (or
  muted year text) for visual rhythm.
- **Teaching**: no structural change beyond removing the sidebar.
- **Curriculum Vitae**:
  - Navigation label changes from `CV` to `Curriculum Vitae`.
  - Page H1 changes from `# CV` to `# Curriculum Vitae`.
  - URL path remains `/cv/` (do not break inbound links).
  - The PDF button is preserved; the "Last updated" line uses a small
    monospace muted treatment.

### 4. Cards, typography, color

Cards:
- Reduce glass blur slightly (`13px` → `11px`) and soften shadows one notch.
  Goal is a more docs-like, less heavy frosted look.
- Keep the top gradient hairline (`::before`), the inner 1px highlight border,
  and the 18px corner radius.
- Hover keeps the existing `translateY(-1px)` lift but with a softer shadow.
- Vertical spacing between successive cards reduced from `1.25rem` to
  `1.1rem`.

Section headings (H2):
- Keep the existing 2.6rem-wide gradient underline anchor.
- Slightly smaller font size (`1.25rem` → `1.18rem`), tighter letter-spacing
  (`0.01em` → `0.005em`).

Typography:
- Body sans stack unchanged.
- Add a monospace stack for accents only:
  `"JetBrains Mono", "SF Mono", "Cascadia Code", ui-monospace, monospace`.
- Monospace usage limited to: News date pills, the CV "Last updated" line,
  optionally the hero meta line.

Color:
- Penn blue `#011f5b` and Penn red `#990000` preserved unchanged.
- Dark-mode palette preserved unchanged.
- Slightly reduce gradient saturation in a few places where Penn red was
  appearing too often.

Explicitly NOT added: scroll progress bar, scanning-line decoration, global
hover sweep on every card. The Research-page sweep stays where it already is.

### 5. Responsive

Breakpoints and behavior:
- `< 480px`: hero stacks — photo centered (~96×96), name/meta below centered,
  icons centered in a row. Card padding `1rem`.
- `480–767px`: hero is horizontal but compact — photo ~110×110 on the left,
  text on the right, content cards use standard padding.
- `768–1023px`: hero photo ~130×130, two-line meta, content centered.
- `1024–1439px`: standard desktop, content max-width ~880px.
- `≥ 1440px`: content max-width ~960px, hero photo can grow to ~150×150,
  card padding slightly larger.

Verification:
- Test viewport widths: 375, 414, 768, 1024, 1366, 1440, 1920.
- Test both `data-theme="light"` and `data-theme="dark"`.
- Verify `backdrop-filter` rendering on Safari (uses `-webkit-` prefix already)
  in addition to Chrome / Firefox.
- Verify the new active-state nav underline and the hero photo ring have
  acceptable contrast in both themes.
- Verify News date pills and CV "Last updated" monospace text remain legible
  in both themes.

## Files expected to change

- `_data/navigation.yml` — relabel `CV` → `Curriculum Vitae`.
- `_pages/cv.md` — H1 to `# Curriculum Vitae`.
- `_pages/about.md` — remove `# Home`, insert hero block markup, keep all
  other content.
- `_pages/education.md` — rename inner H2 to `## Academic Background`.
- `_layouts/default.html` — remove the sidebar include for inner pages (or
  conditionally hide), since the sidebar is no longer used.
- `_includes/sidebar.html` — superseded by the hero (kept on disk but not
  rendered on any page, or removed entirely).
- `_includes/masthead.html` — minor markup if needed; mostly CSS work.
- `assets/css/main.scss` — the bulk of the visual changes.
- `assets/js/_main.js` (or a small inline script in `scripts.html`) — add the
  scroll-state shadow toggle on the masthead.

## Risks

- The site uses the Minimal Mistakes / academicpages base theme; removing the
  sidebar must not break the existing two-column grid (`_sass/_sidebar.scss`,
  Susy grid in `_layouts/default.html`).
- The CV URL path must remain `/cv/` regardless of the new label.
- Backdrop-filter cards under reduced-motion or coarse-pointer devices already
  have overrides; those overrides must not be removed inadvertently.
- News and CV monospace must not regress on Windows/Linux where the listed
  fonts may not be installed (the fallback `ui-monospace, monospace` covers
  this; visually verify).
