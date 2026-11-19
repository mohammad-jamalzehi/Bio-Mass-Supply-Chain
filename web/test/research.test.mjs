import assert from 'node:assert/strict'
import { existsSync, readFileSync } from 'node:fs'
import { test } from 'node:test'
import React from 'react'
import { renderToStaticMarkup } from 'react-dom/server'
import Markdown from 'react-markdown'
import math from 'remark-math'
import katex from 'rehype-katex'

const read = (path) => readFileSync(new URL(path, import.meta.url), 'utf8')
const equations = JSON.parse(read('../src/content/equations.json'))
const results = JSON.parse(read('../src/content/results.json'))
test('all 13 equations map to existing Python locations and render accessible mathematics', () => {
  assert.deepEqual(equations.map(({ id }) => id), Array.from({ length: 13 }, (_, i) => i + 1))
  for (const { latex, source, title } of equations) {
    const [path, line] = source.split('#L')
    assert.ok(read('../../' + path).split(/\r?\n/)[Number(line) - 1]?.trim(), title)
    const rendered = renderToStaticMarkup(React.createElement(Markdown,
      { remarkPlugins: [math], rehypePlugins: [katex] }, '$1+1$\n\n$$\n' + latex + '\n$$'))
    assert.match(rendered, /katex-display/)
    assert.match(rendered, /<math/)
    assert.doesNotMatch(rendered, /katex-error/)
  }
})
test('historical numbers and synthetic provenance are distinct', () => {
  const data = read('../src/content/data.md')
  for (const value of ['12,624', '1,555.2', '1,036.8', '0.094', '0.056', 'Synthetic inputs']) assert.ok(data.includes(value))
  assert.deepEqual(results.fixtures.map((row) => row.totals.cost), [33, 5, 8.4, 10])
  for (const row of results.fixtures) {
    assert.equal(row.data_kind, 'synthetic_test_fixture')
    assert.equal(row.solver.status, 0)
    assert.equal(row.solver.mip_gap, 0)
    assert.equal(row.solver.best_cost_bound, row.totals.cost)
    assert.ok(row.days.length)
  }
  assert.equal(results.annual.filter((row) => row.days_requested === 365).length, 3)
  for (const row of results.annual.filter((row) => row.days_requested === 365)) {
    assert.equal(row.campus_inputs_validated, false)
    assert.equal(row.validation, 'audited_feasible_incumbent')
    assert.equal(row.solver.status_code, 1)
  }
})
test('paper links use the Pages base and preserve author order', () => {
  const app = read('../src/main.tsx')
  assert.match(app, /github\.com\/mohammad-jamalzehi\/Bio-Mass-Supply-Chain/)
  assert.match(app, /import\.meta\.env\.BASE_URL/)
  assert.match(app, /'Mohamad Jamalzehi', 'Kosar Ghazali'/)
  assert.match(read('../index.html'), /name="author" content="Mohamad Jamalzehi, Kosar Ghazali"/)
  for (const name of ['biomass-supply-chain.pdf', 'manuscript.tex', 'references.bib']) {
    assert.ok(existsSync(new URL('../public/research/' + name, import.meta.url)), name)
  }
  assert.equal(read('../public/research/manuscript.tex'), read('../../paper/manuscript.tex'))
  assert.equal(read('../public/research/references.bib'), read('../../paper/references.bib'))
})
test('Vite is configured for the repository Pages path', () => {
  assert.match(read('../vite.config.ts'), /base: '\/Bio-Mass-Supply-Chain\/'/)
})
