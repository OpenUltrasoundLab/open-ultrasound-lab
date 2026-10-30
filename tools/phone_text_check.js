// Flags canvas text drawn outside its canvas at phone width. Usage: node tools/phone_text_check.js public [width] > report.json
const { chromium } = require('playwright'); const fs = require('fs'); const path = require('path');
const root = path.resolve(process.argv[2] || 'public'); const W = +(process.argv[3] || 390);
const pages = []; (function walk(d) { for (const f of fs.readdirSync(d)) { const p = path.join(d, f); if (fs.statSync(p).isDirectory()) { if (!/vendor|matlab/.test(f)) walk(p); } else if (f === 'index.html' && /modules/.test(p)) pages.push(p); } })(root);
const INIT = `(() => { const P = CanvasRenderingContext2D.prototype, of = P.fillText; window.__ovf = [];
  P.fillText = function (t, x, y, mw) { try { const c = this.canvas, m = this.getTransform(), w0 = this.measureText(t).width, w = (mw !== undefined && mw < w0) ? mw : w0, sq = w / w0, al = this.textAlign;
    let l = al === 'center' ? x - w / 2 : (al === 'right' || al === 'end') ? x - w : x; let r = l + w;
    const sx = m.a || 1, d = window.devicePixelRatio || 1; const L = (l * sx + m.e) / d, R = (r * sx + m.e) / d, cw = c.clientWidth;
    if (cw > 0 && (L < -2 || R > cw + 2 || sq < 0.8)) { const id = c.id || (c.closest('[id]') || {}).id || '?'; window.__ovf.push({ id, t: String(t).slice(0, 60), L: Math.round(L), R: Math.round(R), cw, sq: Math.round(sq * 100) / 100 }); } } catch (e) {}
    return of.apply(this, arguments); }; })();`;
(async () => { const b = await chromium.launch(); const out = {};
  for (const pg of pages) { const p = await b.newPage({ viewport: { width: W, height: 800 }, deviceScaleFactor: 2 }); await p.addInitScript(INIT);
    await p.goto('file://' + pg); await p.waitForTimeout(2500);
    const h = await p.evaluate(() => document.body.scrollHeight); for (let y = 0; y < h; y += 700) { await p.evaluate(v => window.scrollTo(0, v), y); await p.waitForTimeout(60); }
    await p.waitForTimeout(800); const v = await p.evaluate(() => window.__ovf);
    const seen = new Set(); out[path.basename(path.dirname(pg))] = v.filter(o => { const k = o.id + '|' + o.t; if (seen.has(k)) return false; seen.add(k); return true; });
    await p.close(); }
  await b.close(); console.log(JSON.stringify(out, null, 1)); })();
