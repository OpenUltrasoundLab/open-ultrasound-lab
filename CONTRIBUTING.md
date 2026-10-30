# Contributing

Thank you for helping. Contributions of every kind are welcome: error reports, technical corrections, better derivations, data, demo code and MATLAB/Octave scripts.

## How to contribute

- **Report an error**: use the "Report an error" link on any module page, or open an issue. Please give the module, the section or equation number, what is wrong, and a source or a short derivation.
- **Suggest content**: open an issue describing the idea, the sources, and where it fits. For larger additions, please open an issue before writing.
- **Pull requests** are welcome for code, scripts, data and technical fixes. Keep them focused on one change.

## Editorial integration

To keep the course consistent in structure, accuracy rules and style, the maintainer (M. Nguyen) integrates and edits all published text. Your technical contribution may be reworded, restructured or merged with existing material before it appears on the site; it will be credited in the errata and review record.

## Quality gates for pull requests

1. `python3 build_site.py --local` completes.
2. `python3 check_site.py --fast` passes (no JS errors, overflow, blank canvases or TeX errors).
3. `python3 check_links.py` reports no broken links.
4. `python3 mtest.py matlab/<script>.m` runs in GNU Octave without error and prints the numbers the page quotes.

By contributing you agree that text and figures are released under CC BY 4.0, MATLAB scripts and in-page demo code under MIT, and changes to the site design and build tools (LICENSE-SITE) under terms that let M. Nguyen distribute them under any license.
