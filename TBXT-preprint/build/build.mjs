import { readFile, writeFile } from 'node:fs/promises';
import { dirname, join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const [manuscript, scratch, modules] = process.argv.slice(2).map(arg => resolve(arg));
const { marked } = await import(pathToFileURL(join(modules, 'marked/lib/marked.esm.js')));
const source = (await readFile(manuscript, 'utf8')).replace(/^\uFEFF/, '');
const tokens = marked.lexer(source, { gfm: true });
const escape = text => text.replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;');
const caption = /^(?:\*\*)?(Figure|Table)\s+[\dS.]+[.:]?\b/i;
let inReferences = false;
let html = '';
for (const [index, token] of tokens.entries()) {
  if (token.type === 'space') continue;
  if (token.type === 'heading') {
    inReferences = /^References\s*$/i.test(token.text);
  }
  const image = token.type === 'paragraph' && token.tokens.length === 1 && token.tokens[0].type === 'image' ? token.tokens[0] : null;
  const prev = tokens.slice(0, index).findLast(entry => entry.type !== 'space');
  const next = tokens.slice(index + 1).find(entry => entry.type !== 'space');
  if (image) {
    const file = pathToFileURL(resolve(dirname(manuscript), image.href));
    const separate = next?.type === 'paragraph' && /^([*]{2})?Figure\s+[\dS.]+[.:]/i.test(next.text);
    const rendered = image.tokens?.length ? marked.Parser.parseInline(image.tokens) : escape(image.text);
    const size = /(?:^|\/)(?:figure|fig)[-_ ]?([1-9])/i.exec(image.href)?.[1];
    html += `<div class="figure${separate ? ' attached' : ''}"${size ? ` data-figure="${size}"` : ''}><img src="${file}" alt="${escape(image.text || 'See figure caption')}">`;
    if (image.text && !separate) html += `<div class="caption">${rendered}</div>`;
    html += '</div>';
    continue;
  }
  if (token.type === 'paragraph') {
    const isCaption = caption.test(token.text);
    const style = isCaption ? (next?.type === 'table' ? 'table-caption' : 'caption') : inReferences ? 'reference' : prev?.type === 'heading' && prev.depth === 1 ? 'metadata' : '';
    html += `<p class="${style}">${marked.Parser.parseInline(token.tokens)}</p>`;
    continue;
  }
  if (token.type === 'table') {
    const weights = token.header.map((header, i) => Math.max(header.text.length, ...token.rows.map(row => row[i]?.text?.length ?? 0)));
    const proportions = weights.map(weight => Math.min(80, Math.max(14, weight)) ** .62);
    const widths = proportions.map(value => 100 * value / proportions.reduce((a, b) => a + b, 0));
    const cols = '<colgroup>' + widths.map(width => `<col style="width:${width.toFixed(2)}%">`).join('') + '</colgroup>';
    html += marked.Parser.parse([token]).replace('<table>', '<table>' + cols);
    continue;
  }
  html += marked.Parser.parse([token]);
}
const css = await readFile(new URL('./publication.css', import.meta.url), 'utf8');
const title = tokens.find(token => token.type === 'heading' && token.depth === 1)?.text || 'Preprint';
await writeFile(join(scratch, 'article.html'), `<!doctype html><html lang="en"><head><meta charset="utf-8"><title>${escape(title)}</title><style>${css}</style></head><body><article>${html}</article></body></html>`, 'utf8');
await writeFile(join(scratch, 'article-tokens.json'), JSON.stringify(tokens), 'utf8');
