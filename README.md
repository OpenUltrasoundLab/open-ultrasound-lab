# Open Ultrasound Lab

A free, open, interactive course on the mathematics, physics, hardware and signal processing of ultrasound imaging: 33 modules in five parts, live in-browser figures, a formula sheet, interview questions with model answers, and runnable MATLAB/Octave scripts. By M. Nguyen, with Claude (Anthropic). Partly reviewed; it will contain errors (see the disclaimer and errata pages on the site).

Site URL: `COURSE_URL` in `site.config.json` (placeholders there are substituted at build time).

## Build

```bash
python3 build_site.py            # writes public/ (MathJax from the jsDelivr CDN) -- the deploy version
python3 build_site.py --local    # same, MathJax served from public/vendor/ (offline testing)
```

Inputs read at build time (never modified): `shell.html`, `m1.html` ... `m33.html`, `rail*.html`, `home*.html`, `matlab/*.m`, `audit_summary.json`, `errata.json`, `site.config.json`, `CITATION.cff`, `build_sheet.py`. Edit `site.config.json` (`GH_USER`, `REPO_URL`, `AUTHOR_URL`, `COURSE_URL`, `VERSION`, `DATE`, `DOI`) before publishing.

## Test

```bash
pip install playwright pyyaml && playwright install chromium
python3 build_site.py --local
python3 check_site.py --fast     # headless render check of every page: JS errors, overflow, blank canvases, MathJax errors, raw TeX
python3 check_links.py           # every internal href/src and #fragment resolves
python3 mtest.py matlab/*.m      # all MATLAB scripts in GNU Octave (headless)
```

CI: `.github/workflows/deploy.yml` builds and deploys `public/` to GitHub Pages on push to `main`; `test.yml` runs the render check, link check and Octave tests on pull requests.

## Contribute

Report errors through the "Report an error" link on any module page or open an issue (module, equation number, what is wrong, a source). For content changes read `CONTRIBUTING.md` (authoring rules and quality gates).

## License

Text and figures: CC BY 4.0 (`LICENSE-CONTENT`). Code (page scripts, MATLAB/Octave, build tools): MIT (`LICENSE`). Cite via `CITATION.cff` or the site's citation page.
