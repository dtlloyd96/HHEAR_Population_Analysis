// Proper CSV line parser. The regex this replaces
//   l.match(/("([^"]|"")*"|[^,]*)/g)
// alternates real fields with empty strings, so taking the first N cells for an
// N-column table silently dropped every second column -- CAS and Category
// vanished from Table S2, and Table S1 rendered with blank columns interleaved.
function parseCSVLine(line) {
  const out = [];
  let cur = '', inQ = false;
  for (let i = 0; i < line.length; i++) {
    const c = line[i];
    if (inQ) {
      if (c === '"') {
        if (line[i + 1] === '"') { cur += '"'; i++; } else { inQ = false; }
      } else cur += c;
    } else if (c === '"') inQ = true;
    else if (c === ',') { out.push(cur); cur = ''; }
    else cur += c;
  }
  out.push(cur);
  return out;
}
module.exports = { parseCSVLine };
