# Contributing

Thank you for helping. The most valuable contribution is a precise error report: module, equation or section number, what is wrong, and a source or short derivation. Use the "Report an error" link on a module page or open an issue.

## Content changes

- Follow these rules: module structure, accuracy rules (one-way vs two-way, dB conventions, FWHM vs sigma), provenance badges on every figure (exact formula, simulation, approximate model, schematic), cite only works you have verified, MathJax markup for equations, base-MATLAB-only scripts that also run in Octave 8.
- Modules live in `mN.html` (fragments assembled by `build_site.py`); MATLAB scripts in `matlab/`. Open an issue before large changes.
- Do not commit `public/`, `build/` or `shots/` (they are generated).

## Quality gates (all must pass; CI runs them on pull requests)

1. `python3 build_site.py --local` completes with no unresolved placeholders.
2. `python3 check_site.py --fast`: every page renders at 1280 px light and 400 px dark with no JS errors, horizontal overflow, blank canvases, `mjx-merror` or raw TeX.
3. `python3 check_links.py`: every internal link and fragment resolves.
4. `python3 mtest.py matlab/*.m`: every script runs in GNU Octave without error, in under about 30 s, and prints the numbers the page text quotes.
5. Numbers you add are checked (Python/Octave/node), and anything unverifiable is labeled as such and added to the errata list.

By contributing you agree that text and figures are released under CC BY 4.0, MATLAB scripts and in-page demo code under MIT, and changes to the site design and build tools (LICENSE-SITE) under terms that let M. Nguyen distribute them under any license.
