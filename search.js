/* Open Ultrasound Lab: client-side search. No dependencies. Index: search-index.json (built by build_site.py). */
(function () {
  'use strict';
  var me = document.currentScript, ROOT = (me && me.getAttribute('data-root')) || '';
  var TYPES = { mod: 'Module', sec: 'Section', demo: 'Demo', eq: 'Equation', rule: 'Rule', cheat: 'Cheat sheet', iv: 'Interview', chk: 'Self-check', gl: 'Glossary' };
  var BASE = { mod: 100, sec: 60, demo: 50, eq: 35, rule: 30, cheat: 25, iv: 22, chk: 20, gl: 70 };
  var GREEK = { 'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta', 'ε': 'epsilon', 'ζ': 'zeta', 'η': 'eta', 'θ': 'theta', 'κ': 'kappa', 'λ': 'lambda', 'μ': 'mu', 'ν': 'nu',
    'ξ': 'xi', 'π': 'pi', 'ρ': 'rho', 'σ': 'sigma', 'ς': 'sigma', 'τ': 'tau', 'φ': 'phi', 'χ': 'chi', 'ψ': 'psi', 'ω': 'omega', 'Δ': 'delta', 'Γ': 'gamma', 'Σ': 'sigma', 'Ω': 'omega',
    'Φ': 'phi', 'Ψ': 'psi', 'Λ': 'lambda', 'Θ': 'theta', 'ϑ': 'theta', 'ϕ': 'phi', 'ϵ': 'epsilon', '×': 'x', '−': '-', '–': '-', '’': "'" };
  var data = null, loading = null, failed = false, items = [], active = -1, box = null, input = null, timer = 0;

  /* Greek letters become their Latin names; accents are stripped; case is folded. */
  function norm(s) {
    var out = '';
    for (var i = 0; i < s.length; i++) {
      var c = s.charAt(i), g = GREEK[c];
      if (g) { out += g; continue; }
      var d = c.normalize ? c.normalize('NFD').replace(/[̀-ͯ]/g, '') : c;
      out += d.toLowerCase();
    }
    return out;
  }
  /* normalised string plus a map from each normalised char back to the original index (for highlighting) */
  function normMap(s) {
    var out = '', map = [];
    for (var i = 0; i < s.length; i++) {
      var n = norm(s.charAt(i));
      for (var k = 0; k < n.length; k++) { out += n.charAt(k); map.push(i); }
    }
    return { s: out, map: map };
  }
  function esc(s) { return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;'); }

  function load() {
    if (data || loading) return loading;
    loading = fetch(ROOT + 'search-index.json').then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); }).then(function (j) {
      data = j;
      data.n = j.e.map(function (e) { return { t: norm(e[1]), x: e[4] ? norm(e[4]) : '' }; });
      data.mt = {}; Object.keys(j.mods).forEach(function (k) { data.mt[k] = norm(j.mods[k][0]); });
    }).catch(function () { failed = true; });
    return loading;
  }

  function search(q) {
    var terms = norm(q).split(/[^a-z0-9.'+\-]+/).filter(Boolean);
    if (!terms.length) return { terms: terms, groups: [] };
    var phrase = terms.join(' '), hits = [];
    for (var i = 0; i < data.e.length; i++) {
      var e = data.e[i], nn = data.n[i], hay = nn.t + ' \u0001 ' + nn.x, ok = true, sc = BASE[e[0]] || 10;
      for (var k = 0; k < terms.length; k++) {
        var p = hay.indexOf(terms[k]);
        if (p < 0) { ok = false; break; }
        if (p < nn.t.length) {
          sc += 4;
          if (p === 0 || /[^a-z0-9]/.test(nn.t.charAt(p - 1))) sc += 6;  /* word start */
          if (p < 40) sc += 3;
        } else sc -= 4;  /* only matched in TeX / where-list text */
      }
      if (!ok) continue;
      if (terms.length > 1 && nn.t.indexOf(phrase) >= 0) sc += 25;
      if (e[0] === 'mod' || e[0] === 'sec') { var tt = e[0] === 'mod' ? data.mt[e[2]] : nn.t; if (terms.every(function (t) { return tt.indexOf(t) >= 0; })) sc += 30; }
      sc -= Math.min(nn.t.length, 400) / 80;
      hits.push({ e: e, sc: sc });
    }
    hits.sort(function (a, b) { return b.sc - a.sc; });
    var order = [], by = {};
    hits.slice(0, 120).forEach(function (h) {
      var n = h.e[2];
      if (!by[n]) { by[n] = { n: n, best: h.sc, list: [] }; order.push(by[n]); }
      if (by[n].list.length < 5) by[n].list.push(h);
    });
    order.sort(function (a, b) { return b.best - a.best; });
    return { terms: terms, groups: order.slice(0, 8) };
  }

  function snippet(text, terms) {
    var nm = normMap(text), first = -1, ranges = [];
    terms.forEach(function (t) {
      var from = 0, p, c = 0;
      while (c++ < 6 && (p = nm.s.indexOf(t, from)) >= 0) {
        if (first < 0 || p < first) first = p;
        ranges.push([nm.map[p], nm.map[p + t.length - 1] + 1]); from = p + t.length;
      }
    });
    var start = 0, end = text.length, W = 150;
    if (first >= 0 && text.length > W) {
      start = Math.max(0, nm.map[first] - 40); end = Math.min(text.length, start + W);
      var sp = text.indexOf(' ', start); if (start > 0 && sp >= 0 && sp < nm.map[first]) start = sp + 1;
    } else if (text.length > W) end = W;
    ranges = ranges.filter(function (r) { return r[0] >= start && r[1] <= end; }).sort(function (a, b) { return a[0] - b[0]; });
    var out = '', pos = start;
    ranges.forEach(function (r) { if (r[0] < pos) return; out += esc(text.slice(pos, r[0])) + '<mark>' + esc(text.slice(r[0], r[1])) + '</mark>'; pos = r[1]; });
    out += esc(text.slice(pos, end));
    return (start > 0 ? '… ' : '') + out + (end < text.length ? ' …' : '');
  }

  var list = document.createElement('div');
  list.id = 'os-list'; list.setAttribute('role', 'listbox'); list.setAttribute('aria-label', 'Search results'); list.hidden = true;
  function mount() { if (!list.parentNode) document.body.appendChild(list); }

  function place() {
    if (!input || list.hidden) return;
    var r = input.getBoundingClientRect(), vw = document.documentElement.clientWidth;
    var w = Math.min(Math.max(r.width, 420), vw - 16), left = Math.max(8, Math.min(r.left, vw - w - 8));
    list.style.left = left + 'px'; list.style.width = w + 'px'; list.style.top = (r.bottom + 4) + 'px';
    list.style.maxHeight = Math.max(160, window.innerHeight - r.bottom - 16) + 'px';
  }
  function setOpen(o) {
    list.hidden = !o;
    document.querySelectorAll('.os-q').forEach(function (i) { i.setAttribute('aria-expanded', o && i === input ? 'true' : 'false'); if (!o) i.removeAttribute('aria-activedescendant'); });
    if (o) place(); else active = -1;
  }
  function setActive(i) {
    var opts = list.querySelectorAll('a[role=option]');
    if (!opts.length) return;
    active = (i + opts.length) % opts.length;
    opts.forEach(function (o, k) { o.setAttribute('aria-selected', k === active ? 'true' : 'false'); });
    input.setAttribute('aria-activedescendant', opts[active].id);
    opts[active].scrollIntoView({ block: 'nearest' });
  }

  function render() {
    mount();
    var q = input.value.trim();
    box.classList.toggle('has-q', !!q);
    if (!q) { setOpen(false); return; }
    if (failed) { list.innerHTML = '<div class="os-none">Search index unavailable (serve the site over http, not file://).</div>'; setOpen(true); return; }
    if (!data) { list.innerHTML = '<div class="os-none">Loading index…</div>'; setOpen(true); return; }
    var res = search(q), html = '', id = 0;
    res.groups.forEach(function (g) {
      var m = data.mods[g.n], nn = ('0' + g.n).slice(-2);
      html += '<div class="os-grp" role="presentation"><b>' + nn + '</b>' + esc(m[0]) + '</div>';
      g.list.forEach(function (h) {
        var e = h.e, href = ROOT + m[1] + (e[3] ? '#' + e[3] : '');
        html += '<a role="option" id="os-opt-' + (id++) + '" aria-selected="false" tabindex="-1" href="' + esc(href) + '"><span class="os-t t-' + e[0] + '">' + TYPES[e[0]] +
          '</span><span class="os-s">' + snippet(e[1], res.terms) + '</span></a>';
      });
    });
    if (!html) html = '<div class="os-none">No results for “' + esc(q) + '”.</div>';
    list.innerHTML = html;
    setOpen(true);
    if (id) setActive(0);
  }

  function activate(inp) { input = inp; box = inp.closest('.os-box'); }
  function init(inp) {
    inp.addEventListener('focus', function () { activate(inp); load().then(function () { if (input === inp) render(); }); if (inp.value.trim()) render(); });
    inp.addEventListener('input', function () { activate(inp); render(); if (!data) load().then(function () { if (input === inp) render(); }); });
    inp.addEventListener('keydown', function (ev) {
      var k = ev.key;
      if (k === 'ArrowDown' || k === 'ArrowUp') {
        if (list.hidden && inp.value.trim()) render();
        if (!list.hidden) { ev.preventDefault(); setActive(active + (k === 'ArrowDown' ? 1 : -1)); }
      } else if (k === 'Enter') {
        var a = list.querySelectorAll('a[role=option]')[active < 0 ? 0 : active];
        if (!list.hidden && a) { ev.preventDefault(); location.href = a.href; }
      } else if (k === 'Escape') {
        if (!list.hidden) { ev.preventDefault(); setOpen(false); } else if (inp.value) inp.value = ''; else inp.blur();
      }
    });
    inp.addEventListener('blur', function () { clearTimeout(timer); timer = setTimeout(function () { if (document.activeElement !== inp && !list.contains(document.activeElement)) setOpen(false); }, 120); });
  }
  function ready() {
    var inputs = document.querySelectorAll('.os-q');
    if (!inputs.length) return;
    inputs.forEach(init);
    list.addEventListener('mousedown', function (e) { e.preventDefault(); });
    window.addEventListener('resize', place);
    window.addEventListener('scroll', place, { passive: true });
    document.addEventListener('keydown', function (e) {
      if (e.key !== '/' || e.ctrlKey || e.metaKey || e.altKey) return;
      var t = e.target, tag = t && t.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || (t && t.isContentEditable)) return;
      var vis = Array.prototype.filter.call(document.querySelectorAll('.os-q'), function (i) { return i.offsetParent !== null; });
      if (!vis.length) return;
      e.preventDefault(); vis[0].focus(); vis[0].select();
    });
    document.addEventListener('click', function (e) { if (!list.hidden && !list.contains(e.target) && !e.target.closest('.os-box')) setOpen(false); });
  }
  /* arriving on a #anchor that is a closed <details> (interview / self-check question): open it */
  function openTarget() {
    var id = decodeURIComponent((location.hash || '').slice(1)), el = id && document.getElementById(id);
    if (el && el.tagName === 'DETAILS') { el.open = true; el.scrollIntoView(); }
  }
  window.addEventListener('hashchange', openTarget);
  window.addEventListener('load', openTarget);
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', ready); else ready();
})();
