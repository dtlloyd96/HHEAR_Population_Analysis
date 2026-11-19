// Science house style for supplementary text, applied by build_science.js to
// everything it prints (Methods, Supplementary Text, captions, table titles).
//
//   sci(text)   mechanical style rules (idempotent):
//     fig. S3 / figs. S3 and S4 (supplementary), Fig. 3B / Figs. 2 and 3 (main),
//     table S7 (supplementary), Table 1 (main), capitalized at sentence start;
//     italic P for P values (*P* markup); US spelling.
//   runs(text)  TextRuns with **bold** and *italic* markup.
//
// Run directly to restyle a source file in place:
//   node science_style.js <file.md | captions.json>
const { TextRun } = require('docx');

const SPELLING = {
  millilitre: 'milliliter', millilitres: 'milliliters', standardisation: 'standardization',
  standardised: 'standardized', standardise: 'standardize', harmonised: 'harmonized',
  harmonise: 'harmonize', harmonisation: 'harmonization', summarised: 'summarized',
  summarise: 'summarize', normalised: 'normalized', normalise: 'normalize', centred: 'centered',
  colour: 'color', colours: 'colors', coloured: 'colored', analysed: 'analyzed', analyse: 'analyze',
  modelled: 'modeled', modelling: 'modeling', grey: 'gray', behaviour: 'behavior',
  organisation: 'organization', recognised: 'recognized', characterise: 'characterize',
  characterised: 'characterized', utilise: 'utilize', minimise: 'minimize', labelled: 'labeled',
  favour: 'favor', programme: 'program', randomised: 'randomized', emphasise: 'emphasize',
  prioritise: 'prioritize', categorise: 'categorize', categorised: 'categorized',
};

function sci(t) {
  let s = t;
  // supplementary figures: plural first, then singular
  s = s.replace(/\b(?:Supplementary |Supplemental )?(?:Figures|Figs\.?)\s+S(\d+)/gi, 'figs. S$1');
  s = s.replace(/\b(?:Supplementary |Supplemental )?(?:Figure|Fig\.?)\s+S(\d+)/gi, 'fig. S$1');
  // main-text figures
  s = s.replace(/\b(?:Figures|Figs\.?)\s+(\d)/g, 'Figs. $1').replace(/\b(?:Figure|Fig\.?|fig\.)\s+(\d)/g, 'Fig. $1');
  // supplementary tables (Supplementary Table 7, Supplemental table 2, Table S8)
  s = s.replace(/\b(?:Supplementary|Supplemental) tables? S?(\d+|\[[^\]]*\])/gi, (m, n) => (/ables/.test(m) ? 'tables S' : 'table S') + n);
  s = s.replace(/\bTables S(\d+|\[)/g, 'tables S$1').replace(/\bTable S(\d+|\[)/g, 'table S$1');
  // sentence start: capitalize
  s = s.replace(/(^|[.!?]\s+|\n)(figs?\.|tables?) S/g, (m, pre, w) => pre + w[0].toUpperCase() + w.slice(1) + ' S');
  // P values
  s = s.replace(/(^|[^*\w])p(\s*[<=>≤≥])/g, '$1*P*$2');
  s = s.replace(/\bp-values?\b/g, (m) => '*P* value' + (m.endsWith('s') ? 's' : ''));
  s = s.replace(/\bpermutation p\b/g, 'permutation *P*');
  // spelling (keep initial capital)
  s = s.replace(/\b[A-Za-z]+\b/g, (w) => {
    const r = SPELLING[w.toLowerCase()];
    if (!r) return w;
    return w[0] === w[0].toUpperCase() ? r[0].toUpperCase() + r.slice(1) : r;
  });
  return s;
}

function runs(text, base = {}) {
  return text.split(/(\*\*[^*]+\*\*|\*[^*\s][^*]*\*)/).filter(Boolean).map((p) => {
    if (p.startsWith('**')) return new TextRun({ ...base, text: p.slice(2, -2), bold: true });
    if (p.startsWith('*') && p.endsWith('*') && p.length > 2) return new TextRun({ ...base, text: p.slice(1, -1), italics: true });
    return new TextRun({ ...base, text: p });
  });
}

module.exports = { sci, runs };

if (require.main === module) {
  const fs = require('fs');
  const file = process.argv[2];
  const old = fs.readFileSync(file, 'utf8');
  let out;
  if (file.endsWith('.json')) {
    const o = JSON.parse(old);
    for (const k of Object.keys(o)) o[k] = sci(o[k]);
    out = JSON.stringify(o, null, 1) + '\n';
  } else {
    // leave the reference list and <!-- notes --> alone
    const cut = old.indexOf('\n## References and notes');
    const head = cut === -1 ? old : old.slice(0, cut);
    out = head.split(/(<!--[\s\S]*?-->)/).map((x) => (x.startsWith('<!--') ? x : sci(x))).join('') +
          (cut === -1 ? '' : old.slice(cut));
  }
  const a = old.split('\n'), b = out.split('\n');
  a.forEach((l, i) => { if (l !== b[i]) console.log(`${i + 1}\n  - ${l.slice(0, 160)}\n  + ${b[i].slice(0, 160)}`); });
  fs.writeFileSync(file, out);
}
