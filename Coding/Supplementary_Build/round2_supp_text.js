// round2_supp_text.js -- Supplementary Text S8 onward, read from
// Output/Paper1_Submission3/Text/Supplementary_Text_Round2_DRAFT.md.
//
// "## S8 Title" starts a section; paragraphs are separated by blank lines; an
// HTML comment that follows a paragraph is returned as that paragraph's note
// (shown in the annotated build, dropped elsewhere). Anything before the first
// "## " heading is ignored.
//
//   readRound2Text(root) -> [{ num: 'S8', title: 'Overlap of ...', paras: [{ text, notes: [] }] }]
const fs = require('fs');
const path = require('path');

function readRound2Text(root) {
  const f = path.join(root, 'Output', 'Paper1_Submission3', 'Text', 'Supplementary_Text_Round2_DRAFT.md');
  if (!fs.existsSync(f)) return [];
  const src = fs.readFileSync(f, 'utf8');
  const sections = [];
  for (const chunk of src.split(/^## /m).slice(1)) {
    const nl = chunk.indexOf('\n');
    const head = chunk.slice(0, nl).trim();
    const m = head.match(/^(S\d+)\s+(.*)$/);
    if (!m) continue;
    const sec = { num: m[1], title: m[2], paras: [] };
    // split into paragraphs and comments, in order
    const parts = chunk.slice(nl + 1).split(/(<!--[\s\S]*?-->)/);
    for (const part of parts) {
      if (part.startsWith('<!--')) {
        const note = part.slice(4, -3).replace(/\s+/g, ' ').trim();
        if (sec.paras.length && note) sec.paras[sec.paras.length - 1].notes.push(note);
        continue;
      }
      for (const p of part.split(/\n\s*\n/)) {
        const text = p.replace(/\s*\n\s*/g, ' ').trim();
        if (text) sec.paras.push({ text, notes: [] });
      }
    }
    sections.push(sec);
  }
  return sections;
}

// *italic* -> plain, for builds whose paragraph helper only understands **bold**
const dropItalics = (t) => t.replace(/(^|[^*])\*([^*\s][^*]*?)\*(?!\*)/g, '$1$2');

module.exports = { readRound2Text, dropItalics };
