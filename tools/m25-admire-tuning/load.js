// Load the page's own engine (M16E from m25.html) and Lab.rng/Lab.fft from shell.html into Node.
const fs = require('fs');
const SITE = './';
function grab(src, startMarker, endMarker) {
  const i = src.indexOf(startMarker); if (i < 0) throw new Error('marker not found: ' + startMarker);
  const j = src.indexOf(endMarker, i); if (j < 0) throw new Error('end marker not found: ' + endMarker);
  return src.slice(i, j + endMarker.length);
}
const shell = fs.readFileSync(SITE + '/shell.html', 'utf8');
const Lab = {};
eval(grab(shell, 'Lab.rng = function (seed) {', '\n  };'));
eval(grab(shell, 'Lab.fft = function (re, im, inverse) {', '\n  };'));
function loadEngine(file) {
  const page = fs.readFileSync(file || SITE + '/m25.html', 'utf8');
  const code = grab(page, 'const M16E = (function () {', '\n})();');
  return new Function('Lab', 'performance', code + '\nreturn M16E;')(Lab, { now: () => Number(process.hrtime.bigint()) / 1e6 });
}
module.exports = { Lab, loadEngine };
// Engine + ADMIRE: from the page if it contains ADM, else from adm_engine.js (development)
function loadAll(file) {
  const page = fs.readFileSync(file || SITE + '/m25.html', 'utf8');
  const code = grab(page, 'const M16E = (function () {', '\n})();');
  let adm = page.includes('const ADM = (function (E) {') ? grab(page, 'const ADM = (function (E) {', '\n})(M16E);') : fs.readFileSync(__dirname + '/adm_engine.js', 'utf8');
  return new Function('Lab', 'performance', code + '\n' + adm + '\nreturn { E: M16E, ADM };')(Lab, { now: () => Number(process.hrtime.bigint()) / 1e6 });
}
module.exports.loadAll = loadAll;
