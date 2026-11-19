// Supplementary materials in the Science template
// (HHEAR/Manuscript/Paper1/Submit3/Science_Supplementary_Materials_Word_template.docx).
//
// Same content as the clean build (build_clean.js): merged Methods, Supplementary
// Text, Figs. S1-S22 and Tables S1-S8, with no editorial notes. Layout follows the
// template: its styles.xml is loaded as-is (Times New Roman 12 pt; SM Heading,
// SM Subheading, SM Text, SM caption), US Letter with 1-inch margins, page number
// in the footer, cover page as in the template, each figure on its own page
// above "Fig. S#." and its caption, each table on its own page below its caption.
// The template's instruction page and the Movie/Audio/Data entries are left out.
//
//   HHEAR_ROOT=/path/to/HHEAR_Analysis node build_science.js
// Writes Output/Supplementary/HHEAR_Supplementary_Science.docx.
const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, ImageRun, AlignmentType, PageBreak, Footer,
  PageNumber, Table, TableRow, TableCell, WidthType, ShadingType,
} = require('docx');
const { parseCSVLine } = require('./csvparse.js');
const { readMethods, readReferences } = require('./round2_methods.js');
const { sci, runs } = require('./science_style.js');

const HERE = __dirname;
const ROOT = process.env.HHEAR_ROOT || path.resolve(HERE, '..', '..');
const FIGS = process.env.HHEAR_SUPP_FIGS || path.join(ROOT, 'Output/Supplementary/figs');
const HN = path.join(ROOT, 'Output/Hero_Numbers');
const OUT = process.env.HHEAR_SUPP_OUT || path.join(ROOT, 'Output/Supplementary');
const TEMPLATE = process.env.HHEAR_SCIENCE_TEMPLATE || path.resolve(ROOT,
  '../HHEAR/Manuscript/Paper1/Submit3/Science_Supplementary_Materials_Word_template.docx');

// ---- the template's styles.xml ------------------------------------------------
// Read straight out of the .docx (a zip) so the template file is the only source.
function zipEntry(file, name) {
  const buf = fs.readFileSync(file);
  const zlib = require('zlib');
  let eocd = buf.length - 22;
  while (buf.readUInt32LE(eocd) !== 0x06054b50) eocd--;
  let p = buf.readUInt32LE(eocd + 16);
  for (let i = 0, n = buf.readUInt16LE(eocd + 10); i < n; i++) {
    const nameLen = buf.readUInt16LE(p + 28), extra = buf.readUInt16LE(p + 30), comment = buf.readUInt16LE(p + 32);
    const entry = buf.toString('utf8', p + 46, p + 46 + nameLen);
    if (entry === name) {
      const off = buf.readUInt32LE(p + 42), method = buf.readUInt16LE(off + 8), size = buf.readUInt32LE(p + 20);
      const start = off + 30 + buf.readUInt16LE(off + 26) + buf.readUInt16LE(off + 28);
      const data = buf.subarray(start, start + size);
      return (method === 8 ? zlib.inflateRawSync(data) : data).toString('utf8');
    }
    p += 46 + nameLen + extra + comment;
  }
  throw new Error(`${name} not found in ${file}`);
}
if (!fs.existsSync(TEMPLATE)) throw new Error('Science template not found: ' + TEMPLATE + ' (set HHEAR_SCIENCE_TEMPLATE)');
const STYLES = zipEntry(TEMPLATE, 'word/styles.xml');

// ---- paragraph helpers in the template's styles --------------------------------
// Ledger text marks current numbers with **; Science style has no such emphasis.
const plain = (t) => t.replace(/\*\*/g, '');
const P = (style, text, opts = {}) => new Paragraph({ style, children: runs(sci(plain(text))), ...opts });
const smHeading = (t, pageBreak = false) => new Paragraph({ style: 'SMHeading',
  children: [...(pageBreak ? [new PageBreak()] : []), new TextRun({ text: t })] });
const smSub = (t) => P('SMSubheading', t);
const smText = (t) => P('SMText', t);
const pageBreak = () => new Paragraph({ children: [new PageBreak()] });
// Caption as Science sets it: the first sentence (the title) and the panel
// labels -- (A), (A to D), (I and J) -- in bold.
function smCaption(text) {
  const t = sci(plain(text));
  const m = t.match(/^(.+?\.)(\s+[\s\S]*)?$/);
  const rest = (m ? m[2] || '' : t).replace(/\(([A-M](?: (?:to|and) [A-M])?)\)/g, '**($1)**');
  return new Paragraph({ style: 'SMcaption', children: [...(m ? runs('**' + m[1] + '**') : []), ...runs(rest)] });
}
const centred = (text, size, extra = {}) => new Paragraph({ alignment: AlignmentType.CENTER,
  children: [new TextRun({ text, size, ...extra })] });
const blank = () => new Paragraph({ children: [] });

// Text block 6.5 x 9 in at 96 px/in; leave room under the figure for its caption.
const MAXW = 624, MAXH = 700;
function figure(file) {
  const buf = fs.readFileSync(path.join(FIGS, file + '.png'));
  const pw = buf.readUInt32BE(16), ph = buf.readUInt32BE(20);
  let w = MAXW, h = Math.round(MAXW * (ph / pw));
  if (h > MAXH) { h = MAXH; w = Math.round(MAXH * (pw / ph)); }
  return new Paragraph({ alignment: AlignmentType.CENTER,
    children: [new ImageRun({ data: buf, type: 'png', transformation: { width: w, height: h } })] });
}

// ---- content -----------------------------------------------------------------
const TITLE = 'Population Structure of the Chemical Exposome Across Global Cohorts';
const AUTHORS = 'Dillon T. Lloyd, Braden T. Tierney, Shakson Isaac, Dimitri Abrahamsson, Kwaku Poku Asante, ' +
  'Siddhi Desai, Peng Gao, Biban Gill, Rima Habre, Darby Jack, Seyram Kaali, Kurunthachalam Kannan, Alison Lee, ' +
  'Donghai Liang, Gary W. Miller, Alison Motsinger-Reif, John L. Pearce, Krystal J. Godri Pollitt, Marina Sirota, ' +
  'Jeanette A. Stingone, Susan L. Teitelbaum, Chirag J. Patel*';
const EMAIL = 'chirag_patel@hms.harvard.edu';

const caps = JSON.parse(fs.readFileSync(path.join(HERE, 'captions.json'), 'utf8'));
const N_FIG = 22;

const TABLES = [
  ['Cohort demographics by child-based or maternal-based cohort.', 'SuppTable1_regenerated.csv', [2400, 1900, 1900, 1700, 1700], ','],
  ['Chemical abbreviation, full name, CAS ID and category.', 'SuppTable2_FINAL.csv', [1500, 4200, 1800, 2100], ','],
  ['Participants by study.', 'SuppTable_Participants_FINAL.tsv', [3000, 2600, 3000], '\t'],
  ['Study location, HHEAR center, analysis year and type.', 'SuppTable_StudyMeta_FINAL.tsv', [1900, 2400, 2100, 1600, 1600], '\t'],
  ['Demographic features available by study.', 'SuppTable_DemoFeatures_FINAL.tsv', [1500, 7100], '\t'],
  ['Chemicals measured by study.', 'SuppTable_AllChems_FINAL.tsv', [1700, 6900], '\t'],
  ['Cross-cohort model specifications.', 'SuppTable_Models_CrossCohort.tsv', [1500, 2200, 4200, 2100], '\t'],
  ['Within-cohort model specifications.', 'SuppTable_Models_WithinCohort.tsv', [1400, 2100, 4400, 2100], '\t'],
  ['Cohorts and participants by continent for each chemical class. Analytes on ≥2 continents are those measured in cohorts on at least two continents. TEDDY, which has sites in the United States and Europe, is counted in both regions, so regional cohort counts can sum to more than the class total; two ZIP participants with no recorded site are included in the class total but not in a region.', 'SuppTable_ClassByContinent.tsv', [1300, 750, 950, 700, 2000, 900, 2300, 900], '\t'],
  ['Questionnaire items describing exposure sources, by cohort. Items with at least 30 responses and at least two distinct values, grouped into nine source domains and paired with the chemical classes each source would plausibly affect. DISCOVER diet and product-use items (unlabeled response codes) and MADRES rice intake (reported by 237 of 241 participants) are excluded. Every item, with its response distribution, is listed in data S1.', 'SuppTable_SourceInventory.tsv', [1100, 600, 3900, 2700], '\t'],
  ['Regional data layers linked to HHEAR locations. Gridded layers were averaged over each location\'s area, weighted by population (GHS-POP R2023A); national layers were assigned to every location in the country. Layers matched to sampling years were weighted by the number of samples collected in each year, and years outside a layer\'s coverage took the nearest available year.', 'SuppTable_RegionalLayers.tsv', [1700, 2500, 1100, 1300, 2800], '\t'],
  ['HHEAR locations, the areas used to link them to regional data, and layer values at each location\'s sampling years. Values are population-weighted area means matched to sampling years and weighted by samples per year (NO2, 2015 only; groundwater arsenic, static). *Sampling windows taken from published cohort descriptions (ATLAA 2014 to 2018; ELEMENT 1997 to 2004, 2008 to 2012 and 2015 to 2018; MEBCO 2008 to 2011; SWSD 2017 to 2019). PM2.5 begins in 1998, so NYUWHS (1985 to 1990) and the earliest ELEMENT samples take the 1998 value. Dashes, not available (national data sources do not report Puerto Rico separately).', 'SuppTable_LocationLinkage.tsv', [1150, 900, 1350, 700, 600, 550, 500, 650, 650, 600, 550, 550, 600], '\t'],
];
const TEXT_W = 9360;   // 6.5 in in twentieths of a point
function table(file, widths, sep) {
  const rows = fs.readFileSync(file, 'utf8').trim().split('\n')
    .map((l) => sep === '\t' ? l.split('\t') : parseCSVLine(l));
  const bad = rows.findIndex((r) => r.length !== widths.length);
  if (bad !== -1) throw new Error(`${path.basename(file)}: row ${bad} has ${rows[bad].length} fields, expected ${widths.length}`);
  const total = widths.reduce((a, b) => a + b, 0);
  const w = widths.map((x) => Math.floor(x * TEXT_W / total));   // fit the 1-inch margins
  return new Table({
    columnWidths: w, width: { size: TEXT_W, type: WidthType.DXA },
    rows: rows.map((cells, ri) => new TableRow({ tableHeader: ri === 0,
      children: w.map((cw, ci) => new TableCell({ width: { size: cw, type: WidthType.DXA },
        shading: ri === 0 ? { type: ShadingType.CLEAR, fill: 'E7E6E6' } : undefined,
        children: [new Paragraph({ children: [new TextRun({ text: (cells[ci] || '').replace(/^"|"$/g, ''),
          bold: ri === 0, size: 16 })] })] })) })) });
}

const children = [];

// cover page, as laid out in the template
children.push(centred('Supplementary Materials for', 36), blank(),
  centred(TITLE, 28, { bold: true }), blank(),
  centred(AUTHORS, 24), blank(),
  centred('*Corresponding author. Email: ' + EMAIL, 24), blank(),
  new Paragraph({ children: [new TextRun({ text: 'This PDF file includes:', bold: true })] }), blank());
const REF_NUMS = fs.readFileSync(path.join(ROOT, 'Output/Paper1_Submission3/Text/Methods_MERGED_DRAFT.md'), 'utf8')
  .split('## References and notes')[1]?.match(/^\d+(?=\. )/gm)?.map(Number) || [];
const coverItems = ['Materials and Methods', 'Supplementary Text', `Figs. S1 to S${N_FIG}`, `Tables S1 to S${TABLES.length}`, 'Data S1'];
if (REF_NUMS.length) coverItems.push(`References (${Math.min(...REF_NUMS)} to ${Math.max(...REF_NUMS)})`);
for (const item of coverItems) {
  children.push(new Paragraph({ indent: { left: 720 }, children: [new TextRun(item)] }));
}

// Materials and Methods
children.push(smHeading('Materials and Methods', true));
children.push(...readMethods(ROOT, { transform: sci, style: { h2: 'SMSubheading', h3: 'SMSubheading', body: 'SMText' } }));

// Supplementary Text (the corrected-sentence ledger, editorial comments stripped)
const blocks = [];
for (const line of fs.readFileSync(path.join(HN, 'supp_corrected_sentences.md'), 'utf8').split('\n')) {
  if (!line.trim()) continue;
  if (line.startsWith('## ')) blocks.push({ heading: line.slice(3).trim(), lines: [] });
  else if (blocks.length) blocks[blocks.length - 1].lines.push(line.trim());
}
const clean = (l) => l.replace(/\s*<!--.*?-->\s*$/, '').trim();
const TITLES = {
  'S1 ': 'Text S1. HHEAR combined cohort', 'S2 ': 'Text S2. Limit of detection analysis',
  'S3 ': 'Text S3. Location associated with exposure levels',
  'S4 ': 'Text S4. Secular time associated with exposure levels',
  'S5 ': 'Text S5. Demographics associated with exposure levels',
  'S6 ': 'Text S6. Quantile power analysis', 'S7 ': 'Text S7. Exposure correlation and the cumulative score',
};
children.push(smHeading('Supplementary Text', true));
// Text S1 to S13, all from Output/Paper1_Submission3/Text/Supplementary_Text_Round2_DRAFT.md
for (const s of require('./round2_supp_text.js').readRound2Text(ROOT)) {
  children.push(smSub(`Text ${s.num}. ${s.title}`));
  for (const p of s.paras) children.push(smText(p.text));
}

// Figures: page break, figure, "Fig. S#.", caption
for (let i = 1; i <= N_FIG; i++) {
  const label = 'S' + i;
  const text = caps[label];
  if (!text) throw new Error('no caption for Fig. ' + label);
  children.push(pageBreak(), figure(label), smHeading(`Fig. ${label}.`), smCaption(text));
}

// Tables: page break, "Table S#.", caption, table
TABLES.forEach(([title, file, widths, sep], i) => {
  children.push(pageBreak(), smHeading(`Table S${i + 1}.`), smCaption(title), blank(),
                table(path.join(HN, file), widths, sep));
});

// Data S1 is a separate file (Output/Paper1_Submission3/Supplementary_Document/DataS1_Source_Items.csv)
children.push(pageBreak(), smHeading('Data S1. (separate file)'),
  smCaption('Questionnaire items describing exposure sources. One row per retained item (373 items, 23 cohorts): cohort, source table, variable, source domain, number of responses, number of distinct values, response distribution and the chemical classes paired with it (table S10).'));

// References cited only in the supplement, numbered on from the main text
if (REF_NUMS.length) {
  children.push(smHeading('References and Notes', true));
  children.push(...readReferences(ROOT, { style: { body: 'SMcaption' } }));
}

const doc = new Document({
  creator: 'HHEAR analysis pipeline', title: 'Supplementary Materials: ' + TITLE,
  externalStyles: STYLES,
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 },
      margin: { top: 1440, right: 1440, bottom: 1440, left: 1440, header: 720, footer: 720 },
      pageNumbers: { start: 1 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ children: [PageNumber.CURRENT] })] })] }) },
    children }],
});
Packer.toBuffer(doc).then((b) => {
  const out = path.join(OUT, 'HHEAR_Supplementary_Science.docx');
  fs.writeFileSync(out, b);
  console.log('written:', out, (b.length / 1024 / 1024).toFixed(1), 'MB');
});
