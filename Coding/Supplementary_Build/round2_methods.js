// Supplementary Methods, read from the author-editable Markdown file
// Output/Paper1_Submission3/Text/Methods_MERGED_DRAFT.md (published Methods as
// corrected, with the Round 2 analyses in place, in Science style). Used by
// build.js, build_clean.js and build_science.js.
//
// Everything above the first '---' is skipped. '## ' -> heading 2, '### ' ->
// heading 3, blank-line-separated paragraphs, **bold** and *italic* kept (Science
// reference numbers are written (*47*)). A paragraph that is a whole <!-- comment -->
// goes to onNote (annotated build) or is dropped.
//
// readMethods(root, opts)     the Methods, up to '## References and notes'
//   opts.onNote(text)         render editorial notes (annotated build only)
//   opts.style                {h2, h3, body} named paragraph styles (Science template)
//   opts.transform(text)      applied to every heading and paragraph (sci())
//   opts.withReferences       also include the reference list (plain builds)
// readReferences(root, opts)  the reference list alone (Science build, end of document)
const fs = require('fs');
const path = require('path');
const { Paragraph, TextRun, HeadingLevel } = require('docx');
const { runs } = require('./science_style.js');

const FILE = 'Output/Paper1_Submission3/Text/Methods_MERGED_DRAFT.md';
const REFS = '## References and notes';

function parse(root, part, { onNote, style, transform = (t) => t } = {}) {
  const md = fs.readFileSync(path.join(root, FILE), 'utf8');
  let lines = md.split('\n');
  lines = lines.slice(lines.findIndex((l) => l.trim() === '---') + 1);
  const cut = lines.findIndex((l) => l.trim() === REFS);
  if (part === 'methods' && cut !== -1) lines = lines.slice(0, cut);
  if (part === 'refs') lines = cut === -1 ? [] : lines.slice(cut + 1);

  const out = [];
  let para = [];
  const flush = () => {
    const t = para.join(' ').trim();
    para = [];
    if (!t) return;
    const c = t.match(/^<!--\s*([\s\S]*?)\s*-->$/);
    if (c) { if (onNote) out.push(onNote(c[1])); return; }
    const text = part === 'refs' ? t : transform(t);
    out.push(style ? new Paragraph({ style: part === 'refs' ? style.ref || style.body : style.body, children: runs(text) })
                   : new Paragraph({ children: runs(text), spacing: { after: 160, line: 300 } }));
  };
  for (const raw of lines) {
    const l = raw.trim();
    if (!l) { flush(); continue; }
    const h = l.match(/^(#{2,3})\s+(.*)$/);
    if (h) {
      flush();
      const text = transform(h[2]);
      if (style) out.push(new Paragraph({ style: h[1] === '##' ? style.h2 : style.h3,
                                          children: [new TextRun({ text, italics: h[1] === '###' })] }));
      else out.push(new Paragraph({ text, heading: h[1] === '##' ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3,
                                    spacing: { before: 240, after: 140 } }));
    } else para.push(l);
  }
  flush();
  return out;
}

function readMethods(root, opts = {}) {
  const m = parse(root, 'methods', opts);
  if (!opts.withReferences) return m;
  const r = parse(root, 'refs', opts);
  if (!r.length) return m;
  const h = opts.style ? new Paragraph({ style: opts.style.h2, children: [new TextRun('References and notes')] })
                       : new Paragraph({ text: 'References and notes', heading: HeadingLevel.HEADING_2 });
  return [...m, h, ...r];
}
const readReferences = (root, opts = {}) => parse(root, 'refs', opts);

module.exports = { readMethods, readReferences };
